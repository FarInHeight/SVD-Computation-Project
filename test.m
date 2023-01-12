clear;

A = readmatrix('datasets/DATA.csv');
A = A(2:end, 2:end);
[U, S, V] = custom_svd(A, 10^-5);
plot_singular_values( diag(S) );