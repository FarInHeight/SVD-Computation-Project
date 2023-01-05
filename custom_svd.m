function [U, S, S1, V] = custom_svd(A, toll)
%CUSTOM_SVD Function to compute the SVD factorization of a matrix
%   Given a matrix A as input and a tollerance for the stopping criterion
%   this function computes the SVD factorization of A.
    [U, S] = singular_vectors(A, toll, true);
    [V, S1] = singular_vectors(A, toll, false);
end

