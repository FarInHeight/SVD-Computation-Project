function [U, S, V] = custom_svd(A, toll)
%CUSTOM_SVD Function to compute the SVD decomposition of a matrix
%   Given a matrix A as input and a tollerance for the stopping criterion,
%   this function computes the SVD decomposition of A.
    [m, n] = size(A);
    [V, S1] = singular_vectors(A, toll, false);

    % sort singular values and right singular vectors if they are not
    % already
    if ~ issorted( diag(S1) )
        [New_Diag, s_order] = sort(diag(S1), 'descend');
        
        S = diag(New_Diag);
        
        new_order = 1:n;
        for h = 1:length(s_order)
            if h ~= s_order(h)
                new_order(h) = s_order(h);
            end
        end
        
        V = V(:, new_order);
    end

    % take only non-zero diagonal entries
    arr = New_Diag;
    arr(arr > 0) = NaN;
    [zero, index] = max(arr, [], 'omitnan');
    if index ~= 1 
        arr = New_Diag(1 : index - 1);
    else
        arr = New_Diag;
    end

    % modify S1 = inv(S) to fit dimensions
    S1 = diag( 1./ arr );
    S1(n, 1) = 0;
    S1(1, n) = 0;

    % compute U
    U = A*V*S1;
end