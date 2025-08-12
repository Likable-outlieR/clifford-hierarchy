This repository contains the code for Appendix B of the paper "Permutation gates in the third level of the Clifford hierarchy" by Zhiyang He, Luke Robitaille, and Xinyu Tan.


The `gott_moch_conj_magma.txt` script checks our counterexample on seven qubits as shown in Section 5.2 and Figure 3. 

The `six_qubits_perm_search.cpp` script conducts the search for any non--semi-Clifford gates on all six-qubit permutations. Running the script should print
```console
starting search

search complete
```
which indicates that all six-qubit permutations are semi-Clifford. 
