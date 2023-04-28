function value = reconstruct_value(resp_loc, lambdas, phases, ...
                                              bloch_modes, c_wbm)
%% Propagation constants
response_UC_x = resp_loc(1);
response_UC_y = resp_loc(2);
nWave = size(lambdas,1);
prop_const = zeros(nWave,1);
for iWave = 1:nWave
    lambda_x = lambdas(iWave,1);
    lambda_y = lambdas(iWave,2);
    phase_x = phases(iWave,1);
    phase_y = phases(iWave,2);
    index_x = phase_x+response_UC_x-1;
    index_y = phase_y+response_UC_y-1;
    prop_const(iWave) = (lambda_x^index_x)*(lambda_y^index_y);
end

bloch_modes_weighted = bloch_modes* diag(c_wbm);
value = bloch_modes_weighted * prop_const;

