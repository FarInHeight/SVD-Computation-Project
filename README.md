# Advanced-Numerical-Methods-Project
This project is based on an implementation of SVD decomposition through the $QR$ method. It is enriched by shifting, Hessemberg reduction, and Givens rotation matrices.

Basically, $$A = U \Sigma V^T$$
$V$ and $\Sigma$ are computed by the $QR$ method transforming $A^T A$ into hessemberg form and applying Givens rotation matrices at each step to obtain $Q$ as a multiplication of them and maintaining the hessemberg form throughout the process. 

Once $V$ is calculated as the product of $Q_i$ matrices, $U$ is computed as $U = A V \Sigma^{-1}$.

A complete overview of the project, with mathematics used, is available in [presentation.pdf](docs/presentation.pdf) (mathematical proofs are omitted; they can be found in the references).