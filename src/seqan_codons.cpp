#include "seqan_codons.hpp"
#include <array>
#include <algorithm>
#include <cassert>
#include <iterator>
#include <iostream>

typedef unsigned char byte;
using seqan::Iupac;
using seqan::IupacString;
using seqan::AminoAcid;
using seqan::Peptide;

namespace codonalign
{

static std::array<AminoAcid,32768> codon_table; //32K table for fasta codon decoding
// codons are encoded as triplets of 5-bit-encoded nucleotides
// (so any codon can be encoded/decoded as a unique 15-bit value)

static char codon_data[] = { //long list of 3+1 characters (codon+translation)
    'A', 'A', 'A', 'K', 'A', 'A', 'C', 'N', 'A', 'A', 'G', 'K', 'A', 'A', 'R', 'K', 'A', 'A', 'T', 'N',
    'A', 'A', 'Y', 'N', 'A', 'C', 'A', 'T', 'A', 'C', 'B', 'T', 'A', 'C', 'C', 'T', 'A', 'C', 'D', 'T',
    'A', 'C', 'G', 'T', 'A', 'C', 'H', 'T', 'A', 'C', 'K', 'T', 'A', 'C', 'M', 'T', 'A', 'C', 'N', 'T',
    'A', 'C', 'R', 'T', 'A', 'C', 'S', 'T', 'A', 'C', 'T', 'T', 'A', 'C', 'V', 'T', 'A', 'C', 'W', 'T',
    'A', 'C', 'Y', 'T', 'A', 'G', 'A', 'R', 'A', 'G', 'C', 'S', 'A', 'G', 'G', 'R', 'A', 'G', 'R', 'R',
    'A', 'G', 'T', 'S', 'A', 'G', 'Y', 'S', 'A', 'T', 'A', 'I', 'A', 'T', 'C', 'I', 'A', 'T', 'G', 'M',
    'A', 'T', 'H', 'I', 'A', 'T', 'M', 'I', 'A', 'T', 'T', 'I', 'A', 'T', 'W', 'I', 'A', 'T', 'Y', 'I',
    'C', 'A', 'A', 'Q', 'C', 'A', 'C', 'H', 'C', 'A', 'G', 'Q', 'C', 'A', 'R', 'Q', 'C', 'A', 'T', 'H',
    'C', 'A', 'Y', 'H', 'C', 'C', 'A', 'P', 'C', 'C', 'B', 'P', 'C', 'C', 'C', 'P', 'C', 'C', 'D', 'P',
    'C', 'C', 'G', 'P', 'C', 'C', 'H', 'P', 'C', 'C', 'K', 'P', 'C', 'C', 'M', 'P', 'C', 'C', 'N', 'P',
    'C', 'C', 'R', 'P', 'C', 'C', 'S', 'P', 'C', 'C', 'T', 'P', 'C', 'C', 'V', 'P', 'C', 'C', 'W', 'P',
    'C', 'C', 'Y', 'P', 'C', 'G', 'A', 'R', 'C', 'G', 'B', 'R', 'C', 'G', 'C', 'R', 'C', 'G', 'D', 'R',
    'C', 'G', 'G', 'R', 'C', 'G', 'H', 'R', 'C', 'G', 'K', 'R', 'C', 'G', 'M', 'R', 'C', 'G', 'N', 'R',
    'C', 'G', 'R', 'R', 'C', 'G', 'S', 'R', 'C', 'G', 'T', 'R', 'C', 'G', 'V', 'R', 'C', 'G', 'W', 'R',
    'C', 'G', 'Y', 'R', 'C', 'T', 'A', 'L', 'C', 'T', 'B', 'L', 'C', 'T', 'C', 'L', 'C', 'T', 'D', 'L',
    'C', 'T', 'G', 'L', 'C', 'T', 'H', 'L', 'C', 'T', 'K', 'L', 'C', 'T', 'M', 'L', 'C', 'T', 'N', 'L',
    'C', 'T', 'R', 'L', 'C', 'T', 'S', 'L', 'C', 'T', 'T', 'L', 'C', 'T', 'V', 'L', 'C', 'T', 'W', 'L',
    'C', 'T', 'Y', 'L', 'G', 'A', 'A', 'E', 'G', 'A', 'C', 'D', 'G', 'A', 'G', 'E', 'G', 'A', 'R', 'E',
    'G', 'A', 'T', 'D', 'G', 'A', 'Y', 'D', 'G', 'C', 'A', 'A', 'G', 'C', 'B', 'A', 'G', 'C', 'C', 'A',
    'G', 'C', 'D', 'A', 'G', 'C', 'G', 'A', 'G', 'C', 'H', 'A', 'G', 'C', 'K', 'A', 'G', 'C', 'M', 'A',
    'G', 'C', 'N', 'A', 'G', 'C', 'R', 'A', 'G', 'C', 'S', 'A', 'G', 'C', 'T', 'A', 'G', 'C', 'V', 'A',
    'G', 'C', 'W', 'A', 'G', 'C', 'Y', 'A', 'G', 'G', 'A', 'G', 'G', 'G', 'B', 'G', 'G', 'G', 'C', 'G',
    'G', 'G', 'D', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'H', 'G', 'G', 'G', 'K', 'G', 'G', 'G', 'M', 'G',
    'G', 'G', 'N', 'G', 'G', 'G', 'R', 'G', 'G', 'G', 'S', 'G', 'G', 'G', 'T', 'G', 'G', 'G', 'V', 'G',
    'G', 'G', 'W', 'G', 'G', 'G', 'Y', 'G', 'G', 'T', 'A', 'V', 'G', 'T', 'B', 'V', 'G', 'T', 'C', 'V',
    'G', 'T', 'D', 'V', 'G', 'T', 'G', 'V', 'G', 'T', 'H', 'V', 'G', 'T', 'K', 'V', 'G', 'T', 'M', 'V',
    'G', 'T', 'N', 'V', 'G', 'T', 'R', 'V', 'G', 'T', 'S', 'V', 'G', 'T', 'T', 'V', 'G', 'T', 'V', 'V',
    'G', 'T', 'W', 'V', 'G', 'T', 'Y', 'V', 'M', 'G', 'A', 'R', 'M', 'G', 'G', 'R', 'M', 'G', 'R', 'R',
    'N', 'N', 'N', 'X', 'R', 'A', 'Y', 'B', 'S', 'A', 'R', 'Z', 'T', 'A', 'A', '.', 'T', 'A', 'C', 'Y',
    'T', 'A', 'G', '.', 'T', 'A', 'R', '.', 'T', 'A', 'T', 'Y', 'T', 'A', 'Y', 'Y', 'T', 'C', 'A', 'S',
    'T', 'C', 'B', 'S', 'T', 'C', 'C', 'S', 'T', 'C', 'D', 'S', 'T', 'C', 'G', 'S', 'T', 'C', 'H', 'S',
    'T', 'C', 'K', 'S', 'T', 'C', 'M', 'S', 'T', 'C', 'N', 'S', 'T', 'C', 'R', 'S', 'T', 'C', 'S', 'S',
    'T', 'C', 'T', 'S', 'T', 'C', 'V', 'S', 'T', 'C', 'W', 'S', 'T', 'C', 'Y', 'S', 'T', 'G', 'A', '.',
    'T', 'G', 'C', 'C', 'T', 'G', 'G', 'W', 'T', 'G', 'T', 'C', 'T', 'G', 'Y', 'C', 'T', 'R', 'A', '.',
    'T', 'T', 'A', 'L', 'T', 'T', 'C', 'F', 'T', 'T', 'G', 'L', 'T', 'T', 'R', 'L', 'T', 'T', 'T', 'F',
    'T', 'T', 'Y', 'F', 'X', 'X', 'X', 'X', 'Y', 'T', 'A', 'L', 'Y', 'T', 'G', 'L', 'Y', 'T', 'R', 'L'
};


static bool is_codon_table_ready = codon_table_init();

unsigned short pack_codon(const Iupac n1, const Iupac n2, const Iupac n3)
{
    const int bpv = seqan::BitsPerValue<Iupac>::VALUE;
    static_assert(bpv*3 <= 8*sizeof(unsigned short), "unsigned short is not big enough to hold codon of T.");
    unsigned short result = (seqan::ordValue(n1) << (bpv*2)) |
                            (seqan::ordValue(n2) << bpv) |
                            seqan::ordValue(n3);
    if(result >= codon_table.size()) {
        std::cerr << "Invalid codon! " << n1 << n2 << n3 << std::endl;
        throw std::runtime_error("Invalid codon");
    }
    return result;
}

bool codon_table_init()
{
    std::fill(std::begin(codon_table), std::end(codon_table), 'X');
    int cdsize = sizeof(codon_data);
    for(int i = 0; i < cdsize; i += 4) {
        unsigned short aacode = pack_codon(codon_data[i], codon_data[i + 1], codon_data[i + 2]);
        assert(codon_table[aacode] == 'X');
        codon_table[aacode] = codon_data[i + 3];
    }
    return true;
}

AminoAcid translate_codon(const Iupac& n1, const Iupac& n2, const Iupac& n3)
{
    return codon_table[pack_codon(n1, n2, n3)];
}

//simple 1st frame forward translation of a given DNA string
// allocate and returns the translation string
Peptide translate_dna(const IupacString& dnastr)
{
    const size_t dnalen = seqan::length(dnastr);
    const size_t aalen = dnalen / 3;
    Peptide r;
    seqan::resize(r, aalen, 'X');
    for(size_t ai = 0; ai < aalen; ai++) {
        const size_t i = ai * 3;
        r[ai] = codon_table[pack_codon(dnastr[i], dnastr[i + 1], dnastr[i + 2])];
    }
    return r;
}

} // namespace codonalign
