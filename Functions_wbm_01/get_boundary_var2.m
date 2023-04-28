function modes_bd = get_boundary_var2(N_UC_x, N_UC_y, lambdas, phases, ...
                                                    bloch_modes, UC_dofs)
nWave = size(lambdas,1);
%% Get local modes
modes_local_B = bloch_modes([UC_dofs.BL, UC_dofs.B],:);
modes_local_T = bloch_modes([UC_dofs.B, UC_dofs.BR],:);
modes_local_L = bloch_modes([UC_dofs.L, UC_dofs.TL],:);
modes_local_R = bloch_modes([UC_dofs.BL, UC_dofs.L],:);
%% Get global modes
for iWave = 1:nWave
    lambda_x = lambdas(iWave,1);
    lambda_y = lambdas(iWave,2);
    phase_x = phases(iWave,1);
    phase_y = phases(iWave,2);
    index = [phase_x: N_UC_x-1+phase_x]';
    prop_B = (lambda_x.^index) * (lambda_y^phase_y);
    modes_global_B(:,iWave) = kron(prop_B, modes_local_B(:,iWave));
    index = [phase_x: N_UC_x-1+phase_x]';
    prop_T = (lambda_x.^index) * (lambda_y^(N_UC_y+phase_y));
    modes_global_T(:,iWave) = kron(prop_T, modes_local_T(:,iWave));
    index = [phase_y:N_UC_y-1+phase_y]';
    prop_L = (lambda_y.^index) * (lambda_x^phase_x);
    modes_global_L(:,iWave) = kron(prop_L, modes_local_L(:,iWave));
    index = [phase_y:N_UC_y-1+phase_y]';
    prop_R = (lambda_y.^index) * (lambda_x^(N_UC_x+phase_x));
    modes_global_R(:,iWave) = kron(prop_R, modes_local_R(:,iWave));
end

modes_bd = [modes_global_B; modes_global_T; modes_global_L; modes_global_R];

