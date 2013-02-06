#include "config.h"
#include "codon_smithwaterman_gotoh.hpp"

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include <iostream>
#include <string>
#include <vector>


using namespace std;

struct Options
{
    std::string ref_fasta_path;
    std::string qry_fasta_path;
    std::string output_prefix;
    seqan::CharString command_line;
    bool is_fastq;
    bool quiet;
    codonalign::ScoringScheme<int, seqan::Blosum62> scoring_scheme;
};

/// Parse arguments from command line, storing the parsed values in \c options.
/// \returns Exit code.
int parse_args(const int argc, const char** argv, Options& options)
{
    using namespace seqan;
    using codonalign::Macse454Default;
    ArgumentParser parser("codon-align");
    addArgument(parser, ArgParseArgument(seqan::ArgParseArgument::STRING, "ref_fasta"));
    addArgument(parser, ArgParseArgument(seqan::ArgParseArgument::STRING, "qry_fastx"));
    addArgument(parser, ArgParseArgument(seqan::ArgParseArgument::STRING, "output_prefix"));

    addOption(parser, ArgParseOption("", "fastq", "Input is FASTQ."));
    addOption(parser, ArgParseOption("q", "quiet", "Align quietly."));
    addOption(parser, ArgParseOption("fo", "gap-open", "", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "gap-open", Macse454Default.gapopen);
    addOption(parser, ArgParseOption("ge", "gap-extend", "", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "gap-extend", Macse454Default.gapextend);
    addOption(parser, ArgParseOption("st", "stop-codon", "", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "stop-codon", Macse454Default.stop);
    addOption(parser, ArgParseOption("fs", "frameshift", "", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "frameshift", Macse454Default.frameshift);

    addUsageLine(parser, "[options] <reference_fasta> <query_fasta> <output_prefix>");
    setDate(parser, __DATE__);
    setVersion(parser, CODON_ALIGN_VERSION);

    const seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if(res != ArgumentParser::PARSE_OK) {
        return res;
    }

    getArgumentValue(options.ref_fasta_path, parser, 0);
    getArgumentValue(options.qry_fasta_path, parser, 1);
    getArgumentValue(options.output_prefix, parser, 2);
    getOptionValue(options.is_fastq, parser, "fastq");
    getOptionValue(options.quiet, parser, "quiet");

    // Scoring
    codonalign::ScoringScheme<int, seqan::Blosum62> &scheme = options.scoring_scheme;
    getOptionValue(scheme.gapopen, parser, "gap-open");
    getOptionValue(scheme.gapextend, parser, "gap-extend");
    getOptionValue(scheme.stop, parser, "stop-codon");
    getOptionValue(scheme.frameshift, parser, "frameshift");

    return 0;
}

template<typename T, typename TName>
void print_alignment(const seqan::Align<T>& align, std::ostream& out, const TName& ref_name, const TName& read_name) {
    using namespace seqan;
    for(unsigned short i = 0; i < 2; ++i) {
        if(i == 0)
            out << '>' << ref_name << std::endl;
        else
            out << '>' << read_name << std::endl;
        auto r = row(align, i);
        for(auto b = begin(r), e = end(r); b != e; ++b) {
            if(isGap(b))
                out << gapValue<char>();
            else
                out << *b;
        }
        out << std::endl;;
    }
}

template<typename T>
void align_cigar(const seqan::Align<T>& alignment, seqan::String<seqan::CigarElement<>>& result)
{
    using namespace seqan;
    auto clip_begin = clippedBeginPosition(row(alignment, 1));
    auto query_length = length(source(row(alignment, 1)));
    if(clip_begin > 0)
        appendValue(result, seqan::CigarElement<>('S', clip_begin));
    // Generate cigar
    auto ref_it = begin(row(alignment, 0)),
         qry_it = begin(row(alignment, 1)),
         ref_end = end(row(alignment, 0));
    char status = 0;
    long count = 0, qry_count = 0;

    for(; ref_it != ref_end; ++ref_it, ++qry_it) {
        if(isGap(ref_it) && isGap(qry_it)) continue;
        char new_status = 'M';
        // Padding
        if(isGap(ref_it)) {
            new_status = 'I';
        }
        else if (isGap(qry_it))
            new_status = 'D';
        if(new_status == status) count++;
        else {
            if(status != 0) {
                appendValue(result, seqan::CigarElement<>(status, count));
                if(status != 'D') qry_count += count;
            }
            status = new_status;
            count = 1;
        }
    }
    assert(qry_it == end(row(alignment, 1)));
    appendValue(result, CigarElement<>(status, count));
    if(status != 'D') qry_count += count;

    if(query_length - clip_begin - qry_count > 0)
        appendValue(result, CigarElement<>('S', query_length - clip_begin - qry_count));

    size_t total_count = 0;
    for(auto i = begin(result), e = end(result); i != e; i++) {
        char op = i->operation;
        if(op != 'D')
            total_count += i->count;
    }

    assert(total_count == query_length);
}

template<typename T>
int calculate_nm(const seqan::Align<T>& alignment)
{
    int result = 0;
    auto ref_it = begin(row(alignment, 0)),
         qry_it = begin(row(alignment, 1)),
         ref_end = end(row(alignment, 0));

    for(; ref_it != ref_end; ++ref_it, ++qry_it) {
        if(isGap(ref_it) && isGap(qry_it)) continue;
        if(isGap(ref_it) || isGap(qry_it) || *qry_it != *ref_it)
            result++;
    }
    return result;
}

template<typename Alphabet>
int perform_alignment(const Options& options)
{
    using namespace seqan;

    SequenceStream ref_stream(options.ref_fasta_path.c_str());
    if (!isGood(ref_stream)) {
        std::cerr << "Could not open" << options.ref_fasta_path << '\n';
        return 1;
    }

    seqan::CharString ref_id;
    seqan::IupacString ref_seq;

    if (readRecord(ref_id, ref_seq, ref_stream) != 0) {
        std::cerr << "Error reading: " << options.ref_fasta_path << '\n';
    }

    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;

    TNameStore      nameStore;
    appendValue(nameStore, ref_id);
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);
    BamHeader header;
    BamHeaderRecord firstRecord;
    firstRecord.type = BAM_HEADER_FIRST;
    appendValue(firstRecord.tags, BamHeaderRecord::TTag("VN", "1.0"));
    appendValue(header.records, firstRecord);

    BamHeaderRecord seq_record;
    seq_record.type = BAM_HEADER_REFERENCE;
    appendValue(seq_record.tags, BamHeaderRecord::TTag("SN", ref_id));
    appendValue(seq_record.tags, BamHeaderRecord::TTag("LN", std::to_string(length(ref_seq))));
    appendValue(header.records, seq_record);
    appendValue(header.sequenceInfos, BamHeader::TSequenceInfo(ref_id, seqan::length(ref_seq)));
    BamHeaderRecord program_record;
    program_record.type = seqan::BAM_HEADER_PROGRAM;
    appendValue(program_record.tags, BamHeaderRecord::TTag("ID", 1));
    appendValue(program_record.tags, BamHeaderRecord::TTag("PN", "codon-align"));
    appendValue(program_record.tags, BamHeaderRecord::TTag("VN", CODON_ALIGN_VERSION));
    appendValue(program_record.tags, BamHeaderRecord::TTag("CL", options.command_line));
    appendValue(header.records, program_record);

    std::fstream qry_fasta(options.qry_fasta_path, std::ios_base::in | std::ios_base::binary);
    if (!qry_fasta.good()) {
        std::cerr << "Could not open" << options.qry_fasta_path << '\n';
        return 1;
    }
    RecordReader<std::fstream, SinglePass<> > qry_reader(qry_fasta);

    seqan::CharString qry_name;
    IupacString qry_seq;
    seqan::String<char> qry_qualities;

    std::fstream out_stream(options.output_prefix + ".sam", std::ios::binary | std::ios::out);
    if(write2(out_stream, header, context, seqan::Sam()) != 0) {
        std::cerr << "Could not write header\n";
        return 1;
    }

    std::ofstream out_nt(options.output_prefix + "_nt.fasta");
    std::ofstream out_aa(options.output_prefix + "_aa.fasta");

    auto read_record = [&qry_name, &qry_seq, &qry_reader, &qry_qualities]() -> int
    {
        if(HasQualities<Alphabet>::VALUE)
            return readRecord(qry_name, qry_seq, qry_qualities, qry_reader, Fastq());

        int result = readRecord(qry_name, qry_seq, qry_reader, Fasta());
        qry_qualities = "*";
        return result;
    };

    size_t c = 0;
    while(read_record() == 0) {
        StringSet<CharString> parts;
        strSplit(parts, qry_name, ' ', true, 1);
        codonalign::CodonAlignment<int> a = codonalign::codon_align_sw(ref_seq, qry_seq, options.scoring_scheme);
        BamAlignmentRecord record;

        record.qName = parts[0];
        record.flag = 0;
        record.rId = 0;
        record.pos = beginPosition(row(a.dna_alignment, 0));
        record.seq = qry_seq;
        record.qual = qry_qualities;
        record.mapQ = 30; // TODO: Assign
        align_cigar(a.dna_alignment, record.cigar);
        BamTagsDict d(record.tags);
        setTagValue(d, "AS", a.max_score);
        setTagValue(d, "NM", calculate_nm(a.dna_alignment));

        if(write2(out_stream, record, context, seqan::Sam()) != 0) {
            std::cerr << "Could not write BAM record\n";
            return 1;
        }

        if(c++ % 10 == 0 && !options.quiet)
            std::cerr << std::setw(10) << c << " " << qry_name << '\r';

        print_alignment(a.dna_alignment, out_nt, ref_id, qry_name);
        print_alignment(a.aa_alignment, out_aa, ref_id, qry_name);
    }

    return 0;

}

int main(const int argc, const char** argv)
{
    Options options;
    int result = parse_args(argc, argv, options);
    if(result) return result;

    for(int i = 0; i < argc; i++) {
        if(i > 0)
            appendValue(options.command_line, ' ');
        append(options.command_line, argv[i]);
    }

    // Parse
    if(options.is_fastq) {
        return perform_alignment<seqan::Dna5Q>(options);
    } else {
        return perform_alignment<seqan::Iupac>(options);
    }
}
