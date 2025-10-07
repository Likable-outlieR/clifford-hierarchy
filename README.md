This repository contains code for the paper "Characterization of permutation gates in the third level of the Clifford hierarchy" by Zhiyang He, Luke Robitaille, and Xinyu Tan.


The `gott_moch_conj_magma.txt` script checks the equation FGF^{-1}=U_3 in Proposition 5.8.

The `six_qubits_perm_search.cpp` script conducts the search for any non--semi-Clifford gates on all six-qubit permutations, 
as discussed in Remark A.7. Running the script should print
```console
starting search

search complete
```
which indicates that all six-qubit permutations are semi-Clifford. 
