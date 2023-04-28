function [modes, lambda_y] = get_Bloch_modes_y_(UC_dofs, D_UC, lambda_x, master_slave, tol, digit)

%% Building the master-slave matrix
Ams = master_slave.Ams;
ms_dof = master_slave.master_dofs;
MS = Ams.c + lambda_x*Ams.x + lambda_x*Ams.xy + Ams.y;
% D_UC = D_UC/(norm(D_UC,'fro')/100);

%% Build eigen value problem
ID_c = [UC_dofs.I, UC_dofs.L, UC_dofs.R, UC_dofs.B, ...
        UC_dofs.BL, UC_dofs.BR];
ID_ys = [UC_dofs.T, UC_dofs.TR, UC_dofs.TL];
ID_i = [ms_dof.I, ms_dof.L];
ID_yo = [ms_dof.B, ms_dof.BL];
A_11 = MS(ID_c, ID_i);
A_12 = MS(ID_c, ID_yo);
A_22 = MS(ID_ys, ID_yo);
D_11 = D_UC(ID_c, ID_c);
D_12 = D_UC(ID_c, ID_ys);
D_21 = D_UC(ID_ys, ID_c);
D_22 = D_UC(ID_ys, ID_ys);

B_11 = A_11'*D_11*A_11;
B_12_c = A_11'*D_11*A_12;
B_12_y = A_11'*D_12*A_22;
B_21_c = A_12'*D_11*A_11;
B_21_my = A_22'*D_21*A_11;
B_22_c = A_12'*D_11*A_12 + A_22'*D_22*A_22;
B_22_y = A_12'*D_12*A_22;
B_22_my = A_22'*D_21*A_12;

C0 = B_22_c - B_21_c*inv(B_11)*B_12_c - B_21_my*inv(B_11)*B_12_y;
Cy = B_22_y - B_21_c*inv(B_11)*B_12_y;
Cmy = B_22_my - B_21_my*inv(B_11)*B_12_c;

scalar = norm(C0, 'fro');
% scalar = 1;

%% Solve 
% mp.Digits(200);
C0_ = mp(C0);
Cy_ = mp(Cy);
Cmy_ = mp(Cmy);
[modes_yo, lambda_y] = polyeig(Cmy_/scalar, C0_/scalar, Cy_/scalar);
modes_yo = double(modes_yo);
lambda_y = double(lambda_y);
[modes_yo, lambda_y] = postprocess_bloch_waves(modes_yo, lambda_y, tol, digit);

% nAllEigs = size(Cmy,1);
% [modes_yo_large, lambda_y_large] = quadeigs(Cmy*scalar, C0*scalar, Cy*scalar, nAllEigs, 'largestabs');
% [modes_yo_large, lambda_y_large] = postprocess_bloch_waves(modes_yo_large, lambda_y_large, tol, digit);
% IDs = find(abs(lambda_y_large)>max(abs(lambda_y)));
% lambda_y_large = lambda_y_large(IDs);
% modes_yo_large = modes_yo_large(:,IDs);
% 
% modes_yo = [modes_yo_large, modes_yo];
% lambda_y = [lambda_y_large; lambda_y];





% nAllEigs = size(Cmy,1);
% [modes_yo_xlarge, lambda_y_xlarge] = quadeigs(Cmy, C0, Cy, nAllEigs, 'largestabs');
% [modes_yo_xlarge, lambda_y_xlarge] = postprocess_bloch_waves(modes_yo_xlarge, lambda_y_xlarge, tol, digit);
% IDs = find(abs(lambda_y_xlarge)>max(abs(lambda_y)));
% lambda_y_xlarge = lambda_y_xlarge(IDs);
% modes_yo_xlarge = modes_yo_xlarge(:,IDs);
% modes_yo = [modes_yo_xlarge, modes_yo];
% lambda_y = [lambda_y_xlarge; lambda_y];

%% plot
mus = log(lambda_y(:))/1i;
% figure;
% s = 50*ones(size(mus));
% scatter(real(mus), imag(mus), s, 'filled')
% xlabel('Re(\mu)')
% ylabel('Im(\mu)')
% zlabel('Frequency')
% title('Dispersion curve')
% xlim([-pi,pi])
% ylim([-tol,tol])
% set(gca, 'FontSize', 20)

%% Reconstruct
MS_c = MS;
MS_c(UC_dofs.T, :) = 0;
MS_c(UC_dofs.TR, :) = 0;
MS_c(UC_dofs.TL, :) = 0;
modes_m = zeros(numel(ID_i), numel(lambda_y));
modes_m(ID_i,:) = -inv(B_11)*(B_12_c*modes_yo + B_12_y*modes_yo*diag(lambda_y));
modes_m(ID_yo,:) = modes_yo;
modes = MS_c*modes_m + (lambda_x*Ams.xy + Ams.y)*modes_m*diag(lambda_y);

modes_m_ref = zeros(numel(ID_i), numel(lambda_y));
modes_yo_ref = ones(size(modes_yo))/sqrt(numel(modes_yo(:,1)));
modes_m_ref(ID_i,:) = -inv(B_11)*(B_12_c*modes_yo_ref + B_12_y*modes_yo_ref*diag(lambda_y));
modes_m_ref(ID_yo,:) = modes_yo_ref;
modes_ref = MS_c*modes_m_ref + (lambda_x*Ams.xy + Ams.y)*modes_m_ref*diag(lambda_y);

%% accuracy check
resi_ = []; resi_abs = [];
nDOF = size(D_UC,1);
lambda_y_array = [];

for nMode = 1:numel(lambda_y)
    lambda_y_ = lambda_y(nMode);
    A_L = Ams.c + 1/lambda_x*Ams.x + 1/(lambda_x*lambda_y_)*Ams.xy + 1/lambda_y_*Ams.y;
    resi = A_L.'*D_UC* modes(:,nMode);
	ref_resi = D_UC* ones(nDOF,1)/nDOF;
	resi_ = [resi_, norm(resi)/norm(ref_resi)];
    resi_abs = [resi_abs, norm(resi)];
    lambda_y_array = [lambda_y_array, lambda_y_];
end
figure
loglog(abs(lambda_y_array), resi_, '*', 'Linewidth', 2)
xlabel('|\lambda_x|')
ylabel('Error')
title('Error of wave functions')
set(gca, 'FontSize', 20)
mus
mus = mus;
