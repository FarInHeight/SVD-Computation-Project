clear;

A = readmatrix('datasets/DATA.csv');
A = A(2:end, 2:end);
[U, S, V] = custom_svd(A, 10^-8);
plot_singular_values( diag(S) );

dim = min( size(A) );
err = zeros(dim, 1); 

for k = 1:dim
    S3 = S(1:k, 1:k);
    U3 = U(:, 1:k);
    V3 = V(:, 1:k);
    err(k) = norm(A - U3*S3*V3', "fro");
end

figure;
plot(err, '--r');
grid on;
xlabel('NUMBER OF RETAINED SINGULAR VALUES');
ylabel('ERROR - FROBENIUS NORM');