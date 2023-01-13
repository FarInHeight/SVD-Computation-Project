function [A, P] = hessemberg(A)
%HESSEMBERG Fuction to compute the Hessember reduction of a matrix
%   Given a matrix A, this function trasform A in a similar hessemberg
%   matrix. It will be produced as output the hessember matrix and the 
%   matrix P used for the transformation.
    n = size(A, 2);
    P = eye(n);
    for k=1:n-2
        sigma = sign(A(k+1, k)) * norm(A(k+1:n, k));
        v = [sigma + A(k+1, k); A(k+2:n, k)]; 
        beta = 1 / (sigma * (sigma + A(k+1, k))); 
        for j=k:n
            tau = beta * (v' * A(k+1:n, j)); 
            A(k+1:n, j) = A(k+1:n, j) - tau * v;
        end
        for j=1:n
            tau = beta * (A(j, k+1:n) * v);
            A(j, k+1:n) = A(j, k+1:n) - tau * v';
            tau = beta * (P(j, k+1:n) * v); 
            P(j, k+1:n) = P(j, k+1:n) - tau * v';
        end
    end
end

