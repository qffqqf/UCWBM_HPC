function [mode_shapes, eigvals] = quadeigs_log(K, C, M, NBranch, sigma)

nDOF = size(K,1);
A_1 = cat(1, K, -C);
A_2 = cat(1, 0*C, -eye(nDOF));
A = cat(2, A_1, A_2);

B_1 = cat(1, 0*C, M);
B_2 = cat(1, eye(nDOF), 0*C);
B = cat(2, B_1, B_2);

A = logm(B\A);
[mode_shapes, eigvals] = eigs(A, eye(2*nDOF), NBranch, sigma);
mode_shapes = mode_shapes(1:nDOF, :);
mode_shapes_norm = ones(nDOF,1)* sqrt(sum(conj(mode_shapes).*mode_shapes, 1));
mode_shapes = mode_shapes./mode_shapes_norm;
eigvals = (diag(eigvals));
eigvals = exp(diag(eigvals));




