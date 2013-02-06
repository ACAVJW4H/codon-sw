#ifndef CODONS_H
#define CODONS_H
#include <ctype.h>
#include <stdexcept>
#include <string>

#include <seqan/sequence.h>

namespace codonalign
{

unsigned short pack_codon(const seqan::Iupac n1, const seqan::Iupac n2, const seqan::Iupac n3);

seqan::AminoAcid translate_codon(const seqan::Iupac& n1, const seqan::Iupac& n2, const seqan::Iupac& n3);

//simple 1st frame forward translation of a given DNA string
seqan::Peptide translate_dna(const seqan::IupacString& dnastr);

bool codon_table_init();

} // namespace codonalign

#endif
