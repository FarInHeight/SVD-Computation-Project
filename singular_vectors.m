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

    % G will hold c and s for each step
    G = zeros(n - 1, 2);
    % G_aux will be equal to [c, -s; s, c] for each step
    G_aux = zeros(2);

    % M will contain the eigenvectors
    M = P;
    err = toll + 1;
    
    while err > toll 
        H1 = H;
        T = H(n, n) * eye(n);

        % shifting
        H = H - T;

        for k = 1:n-1
            [G(k, 1), G(k, 2)] = givens(H, k, k+1);
            G_aux(1, 1) = G(k, 1);
            G_aux(1, 2) = - G(k, 2); 
            G_aux(2, 1) = G(k, 2);
            G_aux(2, 2) = G(k, 1);

            % nullify the elements below the diagonal
            H(k:k+1, k:n) = G_aux * H(k:k+1, k:n);
        end

        for k = 1:n-1
            G_aux(1, 1) = G(k, 1);
            G_aux(1, 2) = G(k, 2); 
            G_aux(2, 1) = - G(k, 2);
            G_aux(2, 2) = G(k, 1);

            % reconstruct the Hessemberg form
            H(1:k+1, k:k+1) = H(1:k+1, k:k+1) * G_aux;
            % update M as M * Q = M * G12' * G23' * ...
            M(1:n, k:k+1) = M(1:n, k:k+1) * G_aux;
        end
        
        % shifting
        H = H + T;

        err = norm(diag(H - H1), 1);
    end

    % calculate square root for singular values
    H = sqrt( diag( diag( H ) ) );

    % discard imaginary parts if too small
    discard_imag = 10^-5;

    if all( imag( diag(H) ) < discard_imag )
        H = real(H);
    end
    if all( imag( diag(M) ) < discard_imag )
        M = real(M);
    end
end