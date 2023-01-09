function [U, S, V] = custom_svd(A, toll)
%CUSTOM_SVD Function to compute the SVD factorization of a matrix
%   Given a matrix A as input and a tollerance for the stopping criterion,
%   this function computes the SVD factorization of A.
    [V, S1] = singular_vectors(A, toll, false);

    [New_Diag, s_order] = sort(diag(S1), 'descend');
    S1 = diag( New_Diag );
    S1(length(S1), 1) = 0;
    
    V = V(s_order,:);

    S = S1;
    S1 = diag( 1./ diag( S1 ) );
    U = A*V*S1;
end