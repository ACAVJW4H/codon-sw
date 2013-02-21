# Codon-Smith-Waterman
**This is alpha-quality software.**

Pairwise alignments of coding sequences, using a modification of the algorithm described in the
[MACSE paper](http://dx.doi.org/10.1371/journal.pone.0022594):

> V. Ranwez, S. Harispe, F. Delsuc, and E. J. P. Douzery 
> “MACSE: Multiple Alignment of Coding SEquences accounting for frameshifts and stop codons.” 
> PloS one, vol. 6, no. 9, p. e22594, Jan. 2011.

# Building

Requires [cmake](http://www.cmake.org), `libbz2-dev`, and a compiler supporting C++11 (tested on g++ 4.6.3, clang++ 3.2).

To build:

    make

The resulting executable will be placed in `_build/Release/codon-sw`

# Usage

    codon-sw <ref.fasta> <query.fast[aq]> <output_prefix>

Generates three files:

* `<output_prefix>.sam` - [SAM](http://samtools.sourceforge.net/)-format nucleotide alignment
* `<output_prefix>_nt.fasta` - Pairwise alignments in FASTA format, alternating between reference and query sequences.
* `<output_prefix>_aa.fasta` - Translation of `output_prefix_nt.fasta`

## Options

To modify scoring parameters, see options in `codon-sw -h`.

# License

BSD 3-clause.

The excellent [SeqAn](http://www.seqan.de/) library is bundled with the source
in the `deps/seqan` subdirectory. See `deps/seqan/LICENSE` for license information.
