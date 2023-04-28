function [modes_xo, lambda_x] = postprocess_bloch_waves(modes_xo, lambda_x, trunc)

if trunc == 100
    nWaves = numel(lambda_x);
else
    nWaves = ceil(trunc*numel(lambda_x)/100);
end
keep_vals_1 = 1:nWaves;
lambda_x = lambda_x(keep_vals_1);
modes_xo = modes_xo(:,keep_vals_1);
