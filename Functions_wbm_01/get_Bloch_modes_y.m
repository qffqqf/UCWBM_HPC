function [modes, lambda_y] = get_Bloch_modes_y(UC_dofs, D_UC, lambda_x, master_slave, trunc)

%% Building the master-slave matrix
Ams = master_slave.Ams;
ms_dof = master_slave.master_dofs;
MS = Ams.c + lambda_x*Ams.x + lambda_x*Ams.xy + Ams.y;

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

B_11 = (A_11'*D_11*A_11);
B_12_c = (A_11'*D_11*A_12);
B_12_y = (A_11'*D_12*A_22);
B_21_c = (A_12'*D_11*A_11);
B_21_my = (A_22'*D_21*A_11);
B_22_c = (A_12'*D_11*A_12 + A_22'*D_22*A_22);
B_22_y = (A_12'*D_12*A_22);
B_22_my = (A_22'*D_21*A_12);

C0 = B_22_c - B_21_c*(B_11\B_12_c) - B_21_my*(B_11\B_12_y);
Cy = B_22_y - B_21_c*(B_11\B_12_y);
Cmy = B_22_my - B_21_my*(B_11\B_12_c);

%% Solve 
scalar = norm(C0);
[modes_yo, lambda_y] = quadeigs_cayley(Cmy/scalar, C0/scalar, Cy/scalar);
[modes_yo, lambda_y] = postprocess_bloch_waves(modes_yo, lambda_y, trunc);

%% Reconstruct
MS_c = MS;
MS_c(UC_dofs.T, :) = 0;
MS_c(UC_dofs.TR, :) = 0;
MS_c(UC_dofs.TL, :) = 0;
modes_m(ID_i,:) = -B_11\(B_12_c*modes_yo + B_12_y*modes_yo*diag(lambda_y));
modes_m(ID_yo,:) = modes_yo;
modes = MS_c*modes_m + (lambda_x*Ams.xy + Ams.y)*modes_m*diag(lambda_y);

