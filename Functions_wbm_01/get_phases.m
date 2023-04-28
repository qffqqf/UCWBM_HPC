function phases = get_phases(N_UC_x, N_UC_y, lambdas)

nWave = size(lambdas,1);
phases = zeros(nWave, 2);
for iWave = 1:nWave
    lambda_x = lambdas(iWave,1);
    lambda_y = lambdas(iWave,2);
    if real(log(lambda_x)) > eps
        phase_x = -N_UC_x;
    else
        phase_x = 0;
    end
    if real(log(lambda_y)) > eps
        phase_y = -N_UC_y;
    else
        phase_y = 0;
    end
    phases(iWave,:) = [phase_x, phase_y];
end