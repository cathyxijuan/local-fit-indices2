# local-fit-indices

This is a github repository for the project with Hao Wu

Abbreviations used in the code:
1) unstr: unstructured information matrix as the weight matrix. 
2) tri: Gamma is computed using the triple product formula (robust to misspecification).



Functions used in the code:
1) list.mean: find a matrix of means based on a list of matrices
2) list.sd: find a matrix of sds based on a list of matrices
3) srmr.adj: compute the new adjusted SRMR
4) srmr.adj.old: compute the old/default adjusted SRMR in step 2 of the two-stage procedure
5) c.adj.val: compute the adjust constant for the hypothesized model for computing RMSEA and CFI
6) rmsea.adj: compute the new adjusted RMSEA
7) cB.adj.val: compute adjust value for the baseline model for computing RMSEA and CFI
8) cfi.adj: compute the new adjusted CFI.
9) simu.fit: generate a list of matrices with fit indices for a given condition. Each simulated dataset has a matrix of fit indices with rows being different kinds of fit indices and columns being different path models. These matrices are combined as a list across repetitions. 
