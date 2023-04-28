clear all
% close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath([pwd,'/Functions_wbm_01'])
addpath([pwd,'/UCs/vbplates/tet'])
addpath([pwd,'/UCs/vbplates/hex'])

tic 
%% UC model setup
% Assuming renumbered nodes [1:max_number_nodes], corresponding to
% the matrices!
filename = 'uclvb_spcoarse';
import_option = 'mat'; %'mat' or 'op4' (note that op4 requires binary reader)
n_dofnodes = 3; % Amount of DOFs per node: shells = 6; solids = 3; solids & shells: = 6
z_id = 3;
savename = 'vbplate_10x10_wbm_lwn';

%% Load single UC model (identical UCs)
disp('Loading UC model ...')
load([filename,'.mat']);
% F_UC_FOM = comsol_matrices.forc;
[K_UC_FOM, M_UC_FOM, C_UC_FOM, UC_nodes, UC_coordinates, L_UC_x, L_UC_y] = import_FE_UC3(filename);
F_UC_FOM = zeros(size(K_UC_FOM,1),1);

%% UC model order reduction setup
reduction_I_tf = 1; 
reduction_A_tf = 1; 
method_ort = 'SVD'; 
tol_SVD = 1e-15; 
n_modes_I = 30; 
n_modes_A = 50; 
damping = 1e-2;  

%% WFEM 
trunc = 100;

%% Loewener MOR
nSamp_init = 4;
nSamp = 2; % Sampling frequency per iteration
nStep = 3;
tol_mf = 1e-2;

%% Finite structure setup
% Composition amount of unit cells
N_UC_x = 10;% Amount of UCs in X-direction
N_UC_y = 10;% Amount of UCs in Y-direction

%% Analysis setup
% Forced response analysis frequencies
freq = 0:2:1000;
omega_range = 2*pi*freq; % [rad/s]

% Single input point force (uz)
force_UC_x = 4; 
force_UC_y = 4; 
picker_F = @(r) norm(r-[L_UC_x,L_UC_y,0])<1e-10;
force_UC_node = pick_nodes(mesh_data, picker_F);
force_UC_node = force_UC_node(1);
% force_UC_node = 2083; % Which node in the specific UC is excited 2083
forc_loc = [force_UC_x, force_UC_y, force_UC_node];

% % Single/Multiple response node (uz)
response_UC_x = 6; % Which UC in the x-direction
response_UC_y = 6; % Which UC in the y-direction
picker_R = @(r) norm(r-[L_UC_x,L_UC_y,0])<1e-10;
response_UC_node = pick_nodes(mesh_data, picker_R);
% response_UC_node = 34; % Which node in the specific UC is the uz response assessed of
resp_loc = [response_UC_x, response_UC_y, response_UC_node];

% Post-processing
checkplots = 0; % 0 or 1, to verify matrix structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE UC MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!! Hard-coded: assign structural damping !!!!!
K_UC_FOM = K_UC_FOM*(1+damping*1i);

% Inspect original UC matrix structure
if checkplots
    figure
    subplot(121)
    spy(K_UC_FOM)
    title('K UC FOM')
    subplot(122)
    spy(M_UC_FOM)
    title('M UC FOM')
end

UC_dofs = calculate_dof_indices(UC_nodes,mesh_data);
n_nodes_UC = length(struct2array(UC_nodes));
n_dofs_UC_FOM = length(K_UC_FOM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE UC MOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply BMS
if reduction_I_tf
    disp('Interior modal reduction UC ...')
    [M_UC_redI, K_UC_redI, C_UC_redI, UC_dofs_redI, Phi_I, Psi_IA, V_UC_I] = ...
        interior_modal_reduction(M_UC_FOM,K_UC_FOM,C_UC_FOM,UC_dofs,n_modes_I);
end

%% Apply GBMS (only applicable if BMS is applied first!)
if reduction_I_tf&&reduction_A_tf % if BMS has already been performed
    % Note here: additional cleaning of the eigenvectors possible by
    % removing zero rows and re-inserting them
    disp('Boundary modal reduction UC ...')
    [M_UC_red,K_UC_red,C_UC_red,UC_dofs_red,V_UC,L_A] = ...
        boundary_modal_reduction(M_UC_FOM,K_UC_FOM,C_UC_FOM,UC_dofs,...
        M_UC_redI,K_UC_redI,UC_dofs_redI,n_modes_A,Phi_I,Psi_IA,method_ort,tol_SVD);
elseif reduction_I_tf&&~reduction_A_tf% If only BMS is to be performed
    M_UC_red = M_UC_redI;
    K_UC_red = K_UC_redI;
    C_UC_red = C_UC_redI;
    UC_dofs_red = UC_dofs_redI;
    V_UC = V_UC_I;
elseif ~reduction_I_tf&&~reduction_A_tf
    V_UC = speye(n_dofs_UC_FOM);
elseif ~reduction_I_tf&&reduction_A_tf
    error('Enable interior modal reduction to apply boundary modal reduction!')
end

%% Continue with FOM or ROM matrices
if reduction_I_tf||reduction_A_tf % UC ROM
    K_UC = K_UC_red;
    M_UC = M_UC_red;
    C_UC = C_UC_red;
    UC_dofs = UC_dofs_red;
    n_dofs_UC = length(K_UC);
else % UC FOM
    K_UC = K_UC_FOM;
    M_UC = M_UC_FOM;
    C_UC = C_UC_FOM;
    n_dofs_UC = n_dofs_UC_FOM;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FRF calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparation for WBM calculations 
UC_dofs.nDOF = size(K_UC,1);
master_slave = get_master_slave_mx(UC_dofs);
freq_check = linspace(freq(1),freq(end),3*numel(freq));
uz_1 = 0*freq_check;
uz_0 = rand(size(freq_check));
error_all = [];
error = tol_mf*1e6;
new_index = round(linspace(5,numel(freq_check)-5,nSamp_init));
omega_LR_index = [];
omega_LR = [];
uz_sample = [];
opts1.SYM = false;
opts2.SYM = true;

%% WBM frequency loop
while error > tol_mf
    for iWave = 1:numel(new_index)
        omega = 2*pi* freq_check(new_index(iWave));
        D_UC = -omega^2*M_UC + K_UC;
        %% WFEM for x & y boundaries
        bloch_modes = [];
        lambdas = [];
        for iUCy = 0:N_UC_y-1
            lambda_y = exp(1i* 2*iUCy*pi/N_UC_y);
            [modes, lambda_x] = get_Bloch_modes_x(UC_dofs, D_UC, lambda_y, master_slave, trunc);
            bloch_modes = [bloch_modes, modes];
            lambdas = [lambdas; [lambda_x, lambda_y*ones(size(lambda_x))]];
        end
        for iUCx = 0:N_UC_x-1
            lambda_x = exp(1i* 2*iUCx*pi/N_UC_x);
            [modes, lambda_y] = get_Bloch_modes_y(UC_dofs, D_UC, lambda_x, master_slave, trunc);
            bloch_modes = [bloch_modes, modes];
            lambdas = [lambdas; [lambda_x*ones(size(lambda_y)), lambda_y]];
        end
        phases = get_phases(N_UC_x, N_UC_y, lambdas);
        %% DFT: particular solution
        [modes_par, lambdas_par, phases_par] = ...
            compute_particular_solution(master_slave, D_UC, V_UC, forc_loc, ...
                                            N_UC_x, N_UC_y, mesh_data, F_UC_FOM, 1);
        %% Pack up wave functions 
        nBloch = size(bloch_modes, 2);
        nPar = size(modes_par, 2);
        modes_all = [bloch_modes, modes_par];
        lambdas_all = [lambdas; lambdas_par];
        phases_all = [phases; phases_par];
        %% Map wave basis
        modes_bd = get_boundary_var(N_UC_x, N_UC_y, lambdas_all, phases_all, ...
                                                     modes_all, UC_dofs);
        %% Solve
        K_wbm = modes_bd(:, 1:nBloch);
        F_wbm = - sum(modes_bd(:, nBloch+1:end), 2);
        if size(K_wbm,1) == size(K_wbm,2)
            [c_wbm,cond] = linsolve(K_wbm,F_wbm, opts1); 
        else
            F_wbm = K_wbm'*F_wbm;
            K_wbm = K_wbm'*K_wbm;
            [c_wbm,cond] = linsolve(K_wbm,F_wbm, opts2); 
        end

        %% Reconstruct
        response_UC_dof = nonzeros(mesh_data.nv_data(response_UC_node,:));
        modes_physical = V_UC(response_UC_dof(end),:)* modes_all;
        c_physical = [c_wbm; ones(nPar,1)];
        uz_sample = [uz_sample, reconstruct_value(resp_loc, lambdas_all, phases_all, ...
                                                    modes_physical, c_physical)];
    end
    omega_LR_index = [omega_LR_index, new_index];
    omega_LR = [omega_LR, 2*pi* freq_check(new_index)];
    uz_1 = frequency_interpolation(omega_LR, 2*pi* freq_check, uz_sample);
    [new_index, error] = adaptive_selection(omega_LR_index, uz_1, uz_0, nSamp, nStep);
    uz_0 = uz_1;
    error_all = [error_all, error];
    %% Plot error
    semilogy(error_all, '-','Linewidth', 2)
    hold on;
    semilogy([1, numel(error_all)+1], [tol_mf,tol_mf], '--', ...
              'Color', [0.6055, 0.6484, 0.6914],'Linewidth', 3)
    hold off;
    title('Matrix-free algorithm')
    xlabel('Iteration step [-]');ylabel('relative difference [-]')
    set(gca, 'FontSize', 20)
    xlim([1, numel(error_all)+1])
    ylim([1e-2*tol_mf, 1e2])
    yticks(10.^[-6:2:2])
    drawnow
end

uz = frequency_interpolation(omega_LR, 2*pi* freq, uz_sample);
nSamp = numel(uz_sample)

%% Post process FRF plot uz component
figure
semilogy(freq,abs(uz), 'Linewidth', 2)
title('Point to point FRF')
xlabel('Frequency [Hz]');ylabel('|u_z| [m]')
set(gca, 'FontSize', 15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timing = toc
uz_wbm = uz;
save([savename,'.mat'],'freq','uz_wbm','timing','nSamp');
