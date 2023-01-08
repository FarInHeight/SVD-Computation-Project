function [U, S, V] = custom_svd(A, toll)
%CUSTOM_SVD Function to compute the SVD factorization of a matrix
%   Given a matrix A as input and a tollerance for the stopping criterion,
%   this function computes the SVD factorization of A.
    [V, S1] = singular_vectors(A, toll, false);

    if min( diag( S1 ) ) < toll
        [U, S] = singular_vectors(A, toll, true);
        minimum = min([length(U), length(V)]);

        if minimum == length(U)
            S = S1(1:minimum, :);
            
            for k = 1:minimum
                B = A';
                constant = B(1, :) * U(:, k);
                if sign( constant ) ~= sign( V(1, k) )
                    V(:, k) = - V(:, k);
                end
            end
        else
            S = S(:, 1:minimum);
    
            for k = 1:minimum
                constant = A(1, :) * V(:, k);
                if sign( constant ) ~= sign( U(1, k) )
                    U(:, k) = - U(:, k);
                end
            end
        end
    else
        S = S1;
        S1 = diag( 1./ diag( S1 ) );
        U = A*V*S1;
    end
end