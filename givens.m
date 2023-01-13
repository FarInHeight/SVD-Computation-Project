function [c, s] = givens(A, i, j)
%GIVENS Function to compute cos and sin of Gij
%   Given a matrix A, this function returns a vector containing the values
%   of givens rotation matrix Gij such that Gij * A is equal to A except
%   for the element (j, i), which will be zero.
     c = abs( A(i, i) ) / sqrt( A(i, i)^2 + A(j, i)^2 );
     s = - sign( A(j, i) / A(i, i) ) * abs( A(j, i) ) / sqrt( A(i, i)^2 + A(j, i)^2 );
end