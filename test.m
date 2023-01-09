clear;

load './datasets/dataset.mat';
[U, S, V] = custom_svd(A, 10^-5);
err = norm(A - U*S*V', "fro");
plot_singular_values( diag(S) );