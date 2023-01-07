function [M, H] = singular_vectors(A, toll, left)
%SINGULAR_VECTORS Function to compute the left or right singular vectors
%   Given a matrix A and a tollerance for the stopping criterion, this
%   function computes the left or right singular vectors of A and its
%   corrisponding singular values. If left is true then U is calculated, V
%   otherwise.
    if left
        [H, P] = hessemberg(A * A.');
    else
        [H, P] = hessemberg(A.' * A);
    end
    
    n = length(H);

    G = zeros(n - 1, 2);
    G_aux = zeros(2);

    M = P;
    err1 = toll + 1;
    err2 = toll + 1;
    
    while err1 > toll || err2 > toll
        H1 = H;
        M1 = M;
        for k = 1:n-1
            [G(k, 1), G(k, 2)] = givens(H, k, k+1);
            G_aux(1, 1) = G(k, 1);
            G_aux(1, 2) = - G(k, 2); 
            G_aux(2, 1) = G(k, 2);
            G_aux(2, 2) = G(k, 1);

            H(k:k+1, k:n) = G_aux * H(k:k+1, k:n);
        end

        for k = 1:n-1
            G_aux(1, 1) = G(k, 1);
            G_aux(1, 2) = G(k, 2); 
            G_aux(2, 1) = - G(k, 2);
            G_aux(2, 2) = G(k, 1);

            H(1:k+1, k:k+1) = H(1:k+1, k:k+1) * G_aux;
            M(1:n, k:k+1) = M(1:n, k:k+1) * G_aux;
        end
        
        err1 = norm(diag(H - H1), 1);
        err2 = norm(max(M - M1), 1);
    end

    H = sqrt( diag( diag( H ) ) );
    if all( imag( diag(H) ) < toll )
        H = real(H);
    end
    if all( imag( diag(M) ) < toll )
        M = real(M);
    end
end