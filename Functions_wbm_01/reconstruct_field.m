function field = reconstruct_field(N_UC_x, N_UC_y, lambdas, phases, ...
                                              bloch_modes, c_wbm)

%% Propagation constants
nWave = size(lambdas,1);
prop_const = zeros(N_UC_x*N_UC_y, nWave);
for iWave = 1:nWave
    lambda_x = lambdas(iWave,1);
    lambda_y = lambdas(iWave,2);
    phase_x = phases(iWave,1);
    phase_y = phases(iWave,2); 
    [yy,xx] = meshgrid(phase_y:N_UC_y-1+phase_y, phase_x:N_UC_x-1+phase_x);
    prop_const(:,iWave) = (lambda_x.^xx(:)).*(lambda_y.^yy(:));
end

% cubic = ktensor(c_wbm, {bloch_modes, prop_const});
% field = double(tensor(cubic));
% field = field(:);

bloch_modes_weighted = bloch_modes* diag(c_wbm);
field = bloch_modes_weighted * prop_const.';
field = field(:);

