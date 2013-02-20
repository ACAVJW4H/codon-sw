# WARNING

This is alpha-quality software.

# Codon-Smith-Waterman

Pairwise alignments of coding sequences, using a modification of the algorithm described in:

> V. Ranwez, S. Harispe, F. Delsuc, and E. J. P. Douzery 
> “MACSE: Multiple Alignment of Coding SEquences accounting for frameshifts and stop codons.” 
> PloS one, vol. 6, no. 9, p. e22594, Jan. 2011.

[link](http://dx.doi.org/10.1371/journal.pone.0022594)

# Building

Requires [cmake](http://www.cmake.org), and a compiler supporting C++11 (tested on g++ 4.6.3, clang++ 3.2).

To build:

    make release

The resulting executable will be placed in `_build/Release/codon-sw`

# Usage

    codon-sw <ref.fasta> <query.fast[aq]> <output_prefix>

Generates three files:

* `<output_prefix>.sam` - [SAM](http://samtools.sourceforge.net/)-format nucleotide alignment
* `<output_prefix>_nt.fasta` - Pairwise alignments.
* `<output_prefix>_aa.fasta` - Translation of `output_prefix_nt.fasta`

## Options

To modify scoring parameters, see options in `codon-sw -h`.

# License

BSD 3-clause.

The excellent [SeqAn](http://www.seqan.de/) library is bundled with the source
- see `deps/seqan/LICENSE` for license information.
