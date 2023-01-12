function [U, S, V] = custom_svd(A, toll)
%CUSTOM_SVD Function to compute the SVD factorization of a matrix
%   Given a matrix A as input and a tollerance for the stopping criterion,
%   this function computes the SVD factorization of A.
    [U, S] = singular_vectors(A, toll, true);
    [V, S1] = singular_vectors(A, toll, false);

    dims = size(A);
    minimum = min(dims);

    if minimum == dims(1)
        [U, V, S] = order_values(S, S1, U, V, 1);
        
        for k = 1:minimum
            B = A';
            constant = B(1, :) * U(:, k);
            if sign( constant ) ~= sign( V(1, k) )
                V(:, k) = - V(:, k);
            end
        end
    else
        [U, V, S] = order_values(S, S1, U, V, 2);

        for k = 1:minimum
            constant = A(1, :) * V(:, k);
            if sign( constant ) ~= sign( U(1, k) )
                U(:, k) = - U(:, k);
            end
        end
    end

    function [U1, V1, New_S] = order_values(S, S1, U, V, type)
        [New_Diag, s_order] = sort(diag(S), 'descend');
        
        S = diag(New_Diag);

        new_order = 1:dims(1);
        for h = 1:length(s_order)
            if h ~= s_order(h)
                new_order(h) = s_order(h);
            end
        end
    
        U1 = U(:, new_order);

        [New_Diag, s_order] = sort(diag(S1), 'descend');
    
        S1 = diag(New_Diag);
        
        new_order = 1:dims(2);
        for h = 1:length(s_order)
            if h ~= s_order(h)
                new_order(h) = s_order(h);
            end
        end
    
        V1 = V(:, new_order);

        if type == 1
            New_S = S1(1:minimum, :);
        else
            New_S = S(:, 1:minimum);
        end
    end
end