Aligns $s_1$ and $s_2$ using the following recursion on three matrices:

* $P$ is filled with the best score ending in a deletion on $s_1$
* $Q$ is filled with best score ending in a deletion on $s_2$
* $D$ is filled the best score overall score from aligning residues from both sequences, $P$, $Q$,
  and 0 (to allow local alignments to terminate).

$$
p_{ij} = \max \begin{cases} d_{i-3,j} + \gamma_1 + g_o,\\
                             p_{i-3,j} + \gamma_1 + g_e,\\
                             \\
                             d_{i-2,j} + \gamma_1 + \delta + g_o,\\
                             p_{i-2,j} + \gamma_1 + \delta + g_e,\\
                             \\
                             d_{i-1,j} + \gamma_1 + \delta + g_o,\\
                             p_{i-1,j} + \gamma_1 + \delta + g_e
                \end{cases}
$$

$$
q_{ij} = \max \begin{cases} d_{i,j-3} + \gamma_2 + g_o, \\
                             q_{i,j-3} + \gamma_2 + g_e, \\
                             \\
                             d_{i,j-2} + \gamma_2 + \delta + g_o, \\
                             q_{i,j-2} + \gamma_2 + \delta + g_e, \\
                             \\
                             d_{i,j-1} + \gamma_2 + \delta + g_o, \\
                             q_{i,j-1} + \gamma_2 + \delta + g_e
                \end{cases}
$$

$$
d_{ij} = \max \begin{cases}
                  0 \\
                  \\
                  d_{i-3,j-3} + \sigma(a_i, b_j),\\
                  \\
                  d_{i-3,j-2} + \gamma_1 + \delta, \\
                  d_{i-3,j-1} + \gamma_1 + \delta,\\
                  \\
                  d_{i-2,j-3} + \gamma_2 + \delta,\\
                  d_{i-1,j-3} + \gamma_2 + \delta,\\
                  \\
                  d_{i-1,j-1} + 2\delta,\\
                  d_{i-1,j-2} + 2\delta,\\
                  d_{i-2,j-1} + 2\delta,\\
                  d_{i-2,j-2} + 2\delta,\\
                  \\
                  p_{i,j},\\
                  q_{i,j}
                \end{cases}
$$

Where:

* $g_o$ is the (amino acid) gap opening penalty
* $g_e$ is the (amino acid) gap extension penalty
* $\delta$ is the frame shift penalty, 
* $\sigma(a, b)$ is the substitution cost between amino acids $a$ and $b$
* $\gamma_i$ is equal to the stop-codon penalty if the current amino acid in
  sequence $i$ encodes a stop codon, 0 otherwise.
