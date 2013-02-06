#ifndef CODON_SMITH_WATERMAN_GOTOH_H
#define CODON_SMITH_WATERMAN_GOTOH_H

#include "seqan_codons.hpp"
#include "matrix.hpp"

#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/sequence.h>

#include <limits>
#include <vector>
#include <utility>
#include <functional>
#include <set>

namespace codonalign
{


template<typename TScore, typename TMatrix>
struct ScoringScheme
{
    /// <b>amino acid</b> Gap open penalty
    TScore gapopen;
    /// <b>amino acid</b> Gap extension penalty
    TScore gapextend;
    /// Cost of a frame shift
    TScore frameshift;
    /// Cost of a stop codon
    TScore stop;
    /// Substitution matrix
    TMatrix substitution;
};

//matrix (BLOSUM 62), gap opening (−7), gap extension (−1), frameshift (−30), and stop codon (−100) except for 454 reads for which lower penalties were assigned to frameshift (−10) and stop codon (−60)
const ScoringScheme<int,seqan::Blosum62> MacseDefault { -7, -1, -30, -100, seqan::Blosum62() };
const ScoringScheme<int,seqan::Blosum62> Macse454Default { -7, -1, -10, -60, seqan::Blosum62() };
//const ScoringScheme<int,seqan::Blosum62> Macse454Default { -10, -1, -10, -10, seqan::Blosum62() };


struct Traceback
{
    char x_offset;
    char y_offset;

    inline bool stop() const { return x_offset == 0 && y_offset == 0; };
};

// Tracebacks - from 0 to 3 nucleotides in either direction
const Traceback T00{0,0},
                T01{0,1},
                T02{0,2},
                T03{0,3},
                T10{1,0},
                T11{1,1},
                T12{1,2},
                T13{1,3},
                T20{2,0},
                T21{2,1},
                T22{2,2},
                T23{2,3},
                T30{3,0},
                T31{3,1},
                T32{3,2},
                T33{3,3};

// Outputs
::std::ostream& operator<<(std::ostream& os, const Traceback& t)
{
    os << "{" << static_cast<unsigned short>(t.x_offset) << ','
        << static_cast<unsigned short>(t.y_offset) << '}';
    return os;
}

template<typename A, typename B>
::std::ostream& operator<<(std::ostream& os, const std::pair<A,B> p)
{
    os << "(" << p.first << ',' << p.second << '}';
    return os;
}

template<typename A>
::std::ostream& operator<<(std::ostream& os, const std::vector<A> v)
{
    os << "(";
    for(const auto &a : v)
        os << a << ',';
    os << ')';
    return os;
}

template<typename T>
::std::ostream& operator<<(std::ostream& os, const codonalign::Matrix<T>& m)
{
    for(size_t i = 0; i < m.nrows(); i++) {
        for(size_t j = 0; j < m.ncols(); j++) {
            os << m(i, j) << '\t';
        }
        os << std::endl;
    }
    return os;  // for multiple << operators.
}

/// returns vector, where v[i] = translation of codon v[i-2:i]
std::vector<seqan::AminoAcid> _translations(const seqan::IupacString s)
{
    const unsigned l = seqan::length(s);
    std::vector<seqan::AminoAcid> result(l, 'X');
    for(unsigned i = 2; i < l; ++i)
        result[i] = translate_codon(s[i-2], s[i-1], s[i]);
    return result;
}

template <typename TScore>
struct CodonAlignment
{
    seqan::Align<seqan::IupacString> dna_alignment;
    seqan::Align<seqan::Peptide> aa_alignment;
    TScore max_score;

    CodonAlignment(const seqan::IupacString& s1, const seqan::IupacString& s2) {
        using namespace seqan;
        resize(rows(dna_alignment), 2);
        assignSource(row(dna_alignment, 0), s1);
        assignSource(row(dna_alignment, 1), s2);
    }
};

template<typename TIter>
size_t max_index(TIter b, TIter e)
{
    return std::max_element(b, e) - b;
}

template<typename TScore>
CodonAlignment<TScore>
traceback(const seqan::IupacString& s1, const seqan::IupacString& s2,
          const codonalign::Matrix<TScore>& d,
          const codonalign::Matrix<Traceback>& trace)
{
    using namespace seqan;
    using std::pair;
    using std::vector;

    const vector<AminoAcid> t1 = _translations(s1),
                            t2 = _translations(s2);
    Peptide aa1, aa2;
    vector<size_t> aa1_gaps, aa2_gaps;
    CodonAlignment<TScore> result(s1, s2);

    // Find maximum score
    pair<size_t, size_t> max_pos = d.rev_index(max_index(d.begin(), d.end()));
    result.max_score = d(max_pos.first, max_pos.second);

    setEndPosition(row(result.dna_alignment, 0), max_pos.first);
    setEndPosition(row(result.dna_alignment, 1), max_pos.second);

    int i = max_pos.first, j = max_pos.second;
    assert(i > 0);
    assert(j > 0);
    assert(i <= length(s1));
    assert(j <= length(s2));

    while(i > 0 || j > 0) {
        assert(i >= 0);
        assert(j >= 0);
        const Traceback& tb = trace(i, j);
        if(tb.stop()) break;

        assert(tb.x_offset + tb.y_offset > 0);

        if(tb.x_offset == 3)
            append(aa1, t1[i-1]);
        if(tb.y_offset == 3)
            append(aa2, t2[j-1]);

        const short diff = tb.y_offset - tb.x_offset;

        if(tb.y_offset == 3 && tb.x_offset != 3) { // Y consumed amino acid
            aa1_gaps.push_back(length(aa1));
            insertGaps(row(result.dna_alignment, 0), i, diff);
        } else if (tb.x_offset == 3 && tb.y_offset != 3) { // X consumed amino acid
            aa2_gaps.push_back(length(aa2));
            insertGaps(row(result.dna_alignment, 1), j, -diff);
        } else if(diff) { // Frameshift in x and/or y - add a codon
            if(tb.x_offset > 0) {
                append(aa1, 'X');
            } else {
                aa1_gaps.push_back(length(aa1));
            }
            if(tb.y_offset > 0) {
                append(aa2, 'X');
            } else {
                aa2_gaps.push_back(length(aa2));
            }
            insertGaps(row(result.dna_alignment, 0), i, 3 - tb.x_offset);
            insertGaps(row(result.dna_alignment, 1), j, 3 - tb.y_offset);
        }

        i -= tb.x_offset;
        j -= tb.y_offset;
    }


    setBeginPosition(row(result.dna_alignment, 0), i);
    setBeginPosition(row(result.dna_alignment, 1), j);

    reverse(aa1);
    reverse(aa2);
    resize(rows(result.aa_alignment), 2);
    assignSource(row(result.aa_alignment, 0), aa1);
    assignSource(row(result.aa_alignment, 1), aa2);
    for(const size_t g : aa1_gaps)
        insertGap(row(result.aa_alignment, 0), length(aa1) - g);
    for(const size_t g : aa2_gaps)
        insertGap(row(result.aa_alignment, 1), length(aa2) - g);

    if(length(row(result.dna_alignment, 0)) != length(row(result.dna_alignment, 1))) {
        std::cerr << length(row(result.dna_alignment, 0)) << ' ' << length(row(result.dna_alignment, 1)) << std::endl;
        std::cerr << result.dna_alignment;
        assert(false);
    }
    if(length(row(result.aa_alignment, 0)) != length(row(result.aa_alignment, 1))) {
        std::cerr << length(row(result.aa_alignment, 0)) << ' ' << length(row(result.aa_alignment, 1)) << std::endl;
        std::cerr << result.aa_alignment;
        assert(false);
    }

    return result;
}

/** \brief Align s1 to s2
 *
 * Calculates the following recursion:
 * \f[
 * P_{i,j} = \max \begin{cases} d_{i-3,j} + stop\_s1 + g_o,\\
 *                              p_{i-3,j} + stop\_s1 + g_e,\\
 *                              \\
 *                              d_{i-2,j} + stop\_s1 + \delta + g_o,\\
 *                              p_{i-2,j} + stop\_s1 + \delta + g_e,\\
 *                              \\
 *                              d_{i-1,j} + stop\_s1 + \delta + g_o,\\
 *                              p_{i-1,j} + stop\_s1 + \delta + g_e
 *                 \end{cases}
 * \f]
 *
 * \f[
 * Q_{i,j} = \max \begin{cases} d_{i,j-3} + stop\_s2 + g_o, \\
 *                              q_{i,j-3} + stop\_s2 + g_e, \\
 *                              \\
 *                              d_{i,j-2} + stop\_s2 + \delta + g_o, \\
 *                              q_{i,j-2} + stop\_s2 + \delta + g_e, \\
 *                              \\
 *                              d_{i,j-1} + stop\_s2 + \delta + g_o, \\
 *                              q_{i,j-1} + stop\_s2 + \delta + g_e
 *                 \end{cases}
 * \f]
 * \f[
 * D_{i,j} = \max \begin{cases}
 *                   0 \\
 *                   d_{i-3,j-3} + \sigma(a_i, b_j),\\
 *                   \\
 *                   d_{i-3,j-2} + stop\_s1 + \delta,\\
 *                   d_{i-3,j-1} + stop\_s1 + \delta,\\
 *                   \\
 *                   d_{i-2,j-3} + stop\_s2 + \delta,
 *                   d_{i-1,j-3} + stop\_s2 + \delta,
 *                   \\
 *                   d_{i-1,j-1} + 2\delta,\\
 *                   d_{i-1,j-2} + 2\delta,\\
 *                   d_{i-2,j-1} + 2\delta,\\
 *                   d_{i-2,j-2} + 2\delta,\\
 *                   \\
 *                   p_{i,j},\\
 *                   q_{i,j}
 *                 \end{cases}
 * \f]
 * Where \f$\delta\f$ is the frame shift penalty, and \f$sigma(a, b)\f$ is the substitution cost between amino acids
 * \f$a\f$ and \f$b\f$.
 */
template<typename TScore, typename TMatrix>
CodonAlignment<TScore> codon_align_sw(const seqan::IupacString& s1, const seqan::IupacString& s2, const ScoringScheme<TScore,TMatrix>& score)
{
    using namespace seqan;
    using codonalign::Matrix;
    using std::vector;

    const unsigned l1 = seqan::length(s1), l2 = seqan::length(s2);
    const std::vector<AminoAcid> t1 = _translations(s1),
                                 t2 = _translations(s2);

    Matrix<TScore> d(l1 + 1, l2 + 1), // Match
                   p(l1 + 1, l2 + 1), // Deletion in query
                   q(l1 + 1, l2 + 1); // Deletion in reference
    Matrix<Traceback> trace(l1 + 1, l2 + 1);
    trace.fill(T00);

    for(long i = 1; i <= l1; i++) {
        for(long j = 1; j <= l2; j++) {
            AminoAcid aa1 = 'X', aa2 = 'X';

            // TODO: CHECK
            if(i > 0) aa1 = t1[i-1];
            if(j > 0) aa2 = t2[j-1];

            TScore stop_s1 = 0, stop_s2 = 0, subst_aa = -50000;

            // Handle stop codons
            if(aa1 == '*') stop_s1 = score.stop;
            if(aa2 == '*') stop_s2 = score.stop;
            if(aa1 == '*' || aa2 == '*') {
                subst_aa = stop_s1 + stop_s2;
            } else {
                // No stop codons
                if(i >= 3 && j >= 3)
                    subst_aa = seqan::score(score.substitution, aa1, aa2);
            }

            assert(stop_s1 <= 0);
            assert(stop_s2 <= 0);
            assert(score.frameshift < 0);
            assert(score.gapopen < 0);
            assert(score.gapextend < 0);

            // Select a value and traceback
            typedef std::pair<TScore,Traceback> TPair;
            typedef typename std::vector<TPair>::iterator TPairIt;
            auto first_comp = [](const TPair& p1, const TPair& p2) -> bool {
                return p1.first < p2.first;
            };

            // Calculate scores
            std::vector<TPair> scores;
            scores.reserve(15);

            // Update P matrix
            scores.clear();

            // Whole AA insertion
            if(i >= 3) {
                scores.emplace_back(d(i-3, j) + stop_s1 + score.gapopen, T30);
                scores.emplace_back(p(i-3, j) + stop_s1 + score.gapextend, T30);
            }
            // Frameshift
            if(i >= 2) {
                scores.emplace_back(d(i-2, j) + stop_s1 + score.frameshift + score.gapopen, T20);
                scores.emplace_back(p(i-2, j) + stop_s1 + score.frameshift + score.gapextend, T20);
            }
            if(i >= 1) {
                scores.emplace_back(d(i-1, j) + stop_s1 + score.frameshift + score.gapopen, T10);
                scores.emplace_back(p(i-1, j) + stop_s1 + score.frameshift + score.gapextend, T10);
            }

            TPair p_max = *std::max_element(std::begin(scores), std::end(scores), first_comp);
            p(i, j) = p_max.first;

            // Update Q matrix
            scores.clear();
            // Whole AA insertion
            if(j >= 3) {
                scores.emplace_back(d(i, j-3) + stop_s1 + score.gapopen, T03);
                scores.emplace_back(q(i, j-3) + stop_s1 + score.gapextend, T03);
            }
            if(j >= 2) {
                scores.emplace_back(d(i, j-2) + stop_s1 + score.frameshift + score.gapopen, T02);
                scores.emplace_back(q(i, j-2) + stop_s1 + score.frameshift + score.gapextend, T02);
            }
            // Frameshift
            if(j >= 1) {
                scores.emplace_back(d(i, j-1) + stop_s1 + score.frameshift + score.gapopen, T01);
                scores.emplace_back(q(i, j-1) + stop_s1 + score.frameshift + score.gapextend, T01);
            }
            TPair q_max = *std::max_element(std::begin(scores), std::end(scores), first_comp);
            q(i, j) = q_max.first;

            // Update D matrix
            scores.clear();
            scores.emplace_back(0, T00);  // Stop

            if(i >= 3 && j >= 3) {
                scores.emplace_back(d(i - 3, j - 3) + subst_aa, T33);
            }

            // Frameshifts on s2
            if(i >= 3) {
                if(j >= 2) scores.emplace_back(d(i-3, j-2) + stop_s1 + score.frameshift, T32);
                if(j >= 1) scores.emplace_back(d(i-3, j-1) + stop_s1 + score.frameshift, T31);
            }

            // Frameshifts on s1
            if(j >= 3) {
                if(i >= 2) scores.emplace_back(d(i-2, j-3) + score.frameshift + stop_s2, T23);
                if(i >= 1) scores.emplace_back(d(i-1, j-3) + score.frameshift + stop_s2, T13);
            }

            // frameshifts in both
            if(i >= 1) {
                if(j >= 1) scores.emplace_back(d(i-1, j-1) + 2*score.frameshift, T11);
                if(j >= 2) scores.emplace_back(d(i-1, j-2) + 2*score.frameshift, T12);
            }
            if(i >= 2) {
                if(j >= 1) scores.emplace_back(d(i-2, j-1) + 2*score.frameshift, T21);
                if(j >= 2) scores.emplace_back(d(i-2, j-2) + 2*score.frameshift, T22);
            }

            scores.push_back(p_max);
            scores.push_back(q_max);

            TPair d_max = *std::max_element(std::begin(scores), std::end(scores), first_comp);

            d(i, j) = d_max.first;
            trace(i, j) = d_max.second;
        }
    }

    return traceback(s1, s2, d, trace);
}

} // namespace codonalign

#endif // CODON_SMITH_WATERMAN_GOTOH_H
