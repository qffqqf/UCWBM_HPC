function [modes, lambda_x] = get_Bloch_modes_x_(UC_dofs, D_UC, lambda_y, master_slave, tol, digit)

%% Building the master-slave matrix
Ams = master_slave.Ams;
ms_dof = master_slave.master_dofs;
MS = Ams.c + lambda_y*Ams.y + lambda_y*Ams.xy + Ams.x;

%% Build eigen value problem
ID_c = [UC_dofs.I, UC_dofs.L, UC_dofs.B, UC_dofs.T, ...
        UC_dofs.BL, UC_dofs.TL];
ID_xs = [UC_dofs.R, UC_dofs.BR, UC_dofs.TR];
ID_i = [ms_dof.I, ms_dof.B];
ID_xo = [ms_dof.L, ms_dof.BL];
A_11 = MS(ID_c, ID_i);
A_12 = MS(ID_c, ID_xo);
A_22 = MS(ID_xs, ID_xo);
D_11 = D_UC(ID_c, ID_c);
D_12 = D_UC(ID_c, ID_xs);
D_21 = D_UC(ID_xs, ID_c);
D_22 = D_UC(ID_xs, ID_xs);

B_11 = A_11'*D_11*A_11;
B_12_c = A_11'*D_11*A_12;
B_12_x = A_11'*D_12*A_22;
B_21_c = A_12'*D_11*A_11;
B_21_mx = A_22'*D_21*A_11;
B_22_c = A_12'*D_11*A_12 + A_22'*D_22*A_22;
B_22_x = A_12'*D_12*A_22;
B_22_mx = A_22'*D_21*A_12;

C0 = B_22_c - B_21_c*inv(B_11)*B_12_c - B_21_mx*inv(B_11)*B_12_x;
Cx = B_22_x - B_21_c*inv(B_11)*B_12_x;
Cmx = B_22_mx - B_21_mx*inv(B_11)*B_12_c;

% scalar = 1;
scalar = norm(C0, 'fro');
%% Solve 
% mp.Digits(200);
C0_ = mp(C0);
Cx_ = mp(Cx);
Cmx_ = mp(Cmx);
[modes_xo, lambda_x] = polyeig(Cmx_/scalar, C0_/scalar, Cx_/scalar);
modes_xo = double(modes_xo);
lambda_x = double(lambda_x);
[modes_xo, lambda_x] = postprocess_bloch_waves(modes_xo, lambda_x, tol, digit);

% nAllEigs = round(size(Cmx,1));
% [modes_xo_large, lambda_x_large] = quadeigs(Cmx*scalar, C0*scalar, Cx*scalar, nAllEigs, 'largestabs');
% [modes_xo_large, lambda_x_large] = postprocess_bloch_waves(modes_xo_large, lambda_x_large, tol, digit);
% IDs = find(abs(lambda_x_large)>max(abs(lambda_x)));
% lambda_x_large = lambda_x_large(IDs);
% modes_xo_large = modes_xo_large(:,IDs);
% modes_xo = [modes_xo_large, modes_xo];
% lambda_x = [lambda_x_large; lambda_x];


% [modes_xo, lambda_x] = quadeigs(Cmx/scalar, C0/scalar, Cx/scalar, nAllEigs, 'largestabs');
% [modes_xo, lambda_x] = postprocess_bloch_waves(modes_xo, lambda_x, tol, digit);
% IDs = find(abs(lambda_x)>exp(1));
% lambda_x_large = lambda_x(IDs);
% modes_xo_large = modes_xo(:,IDs);
% 
% modes_xo = [modes_xo_large, modes_xo_small];
% lambda_x = [lambda_x_large; lambda_x_small];


% nAllEigs = size(Cmx,1);
% [modes_xo_xlarge, lambda_x_xlarge] = quadeigs(Cmx, C0, Cx, nAllEigs, 'largestabs');
% [modes_xo_xlarge, lambda_x_xlarge] = postprocess_bloch_waves(modes_xo_xlarge, lambda_x_xlarge, tol, digit);
% IDs = find(abs(lambda_x_xlarge)>max(abs(lambda_x)));
% lambda_x_xlarge = lambda_x_xlarge(IDs);
% modes_xo_xlarge = modes_xo_xlarge(:,IDs);
% modes_xo = [modes_xo_xlarge, modes_xo];
% lambda_x = [lambda_x_xlarge; lambda_x];

% modes_xo = [modes_xo_xlarge, modes_xo_large, modes_xo_small];
% lambda_x = [lambda_x_xlarge; lambda_x_large; lambda_x_small];

% [modes_xo, lambda_x] = quadeigs(Cmx, C0, Cx, 100, 'smallestabs');
% [modes_xo, lambda_x] = postprocess_bloch_waves(modes_xo, lambda_x, tol, digit);
% IDs = find(abs(lambda_x) < min(abs(lambda_x_small)));
% lambda_x_xsmall = lambda_x(IDs);
% modes_xo_xsmall = modes_xo(:,IDs);
% 
% modes_xo = [modes_xo_xlarge, modes_xo_large, modes_xo_small, modes_xo_xsmall];
% lambda_x = [lambda_x_xlarge; lambda_x_large; lambda_x_small; lambda_x_xsmall];

%% plot 
mus = log(lambda_x(:))/1i;
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
% hold on;

%% Reconstruct
MS_c = MS;
MS_c(UC_dofs.TR, :) = 0;
MS_c(UC_dofs.R, :) = 0;
MS_c(UC_dofs.BR, :) = 0;
modes_m = zeros(numel(ID_i), numel(lambda_x));
modes_m(ID_i,:) = -inv(B_11)*(B_12_c*modes_xo + B_12_x*modes_xo*diag(lambda_x));
modes_m(ID_xo,:) = modes_xo;
modes = MS_c*modes_m + (lambda_y*Ams.xy + Ams.x)*modes_m*diag(lambda_x);

modes_m_ref = zeros(numel(ID_i), numel(lambda_x));
modes_xo_ref = ones(size(modes_xo))/sqrt(numel(modes_xo(:,1)));
modes_m_ref(ID_i,:) = -inv(B_11)*(B_12_c*modes_xo_ref + B_12_x*modes_xo_ref*diag(lambda_x));
modes_m_ref(ID_xo,:) = modes_xo_ref;
modes_ref = MS_c*modes_m_ref + (lambda_y*Ams.xy + Ams.x)*modes_m_ref*diag(lambda_x);

%% accuracy check
resi_ = []; resi_abs = [];
lambda_x_array = [];
nDOF = size(D_UC,1);
for nMode = 1:numel(lambda_x)
    lambda_x_ = lambda_x(nMode);
    A_L = Ams.c + 1/lambda_y*Ams.y + 1/(lambda_y*lambda_x_)*Ams.xy + 1/lambda_x_*Ams.x;
    resi = A_L.'*D_UC* modes(:,nMode);
%     ref_resi = A_L.'*D_UC* (modes_ref(:,nMode));
    ref_resi = D_UC* ones(nDOF,1)/nDOF;
	resi_ = [resi_, norm(resi)/norm(ref_resi)];
    resi_abs = [resi_abs, norm(resi)];
    lambda_x_array = [lambda_x_array, lambda_x_];
end

figure
loglog(abs(lambda_x_array), resi_, '*', 'Linewidth', 2)
xlabel('|\lambda_x|')
ylabel('Error')
title('Error of wave functions')
set(gca, 'FontSize', 20)

mus
mus = mus;
