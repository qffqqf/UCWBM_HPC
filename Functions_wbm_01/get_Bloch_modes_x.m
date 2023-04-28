function [modes, lambda_x] = get_Bloch_modes_x(UC_dofs, D_UC, lambda_y, master_slave, trunc)

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

B_11 = (A_11'*D_11*A_11);
B_12_c = (A_11'*D_11*A_12);
B_12_x = (A_11'*D_12*A_22);
B_21_c = (A_12'*D_11*A_11);
B_21_mx = (A_22'*D_21*A_11);
B_22_c = (A_12'*D_11*A_12 + A_22'*D_22*A_22);
B_22_x = (A_12'*D_12*A_22);
B_22_mx = (A_22'*D_21*A_12);

C0 = B_22_c - B_21_c*(B_11\B_12_c) - B_21_mx*(B_11\B_12_x);
Cx = B_22_x - B_21_c*(B_11\B_12_x);
Cmx = B_22_mx - B_21_mx*(B_11\B_12_c);

%% Solve 
scalar = norm(C0);
[modes_xo, lambda_x] = quadeigs_cayley(Cmx/scalar, C0/scalar, Cx/scalar);
[modes_xo, lambda_x] = postprocess_bloch_waves(modes_xo, lambda_x, trunc);

%% Reconstruct     
MS_c = MS;
MS_c(UC_dofs.TR, :) = 0;
MS_c(UC_dofs.R, :) = 0;
MS_c(UC_dofs.BR, :) = 0;
modes_m(ID_i,:) = -B_11\(B_12_c*modes_xo + B_12_x*modes_xo*diag(lambda_x));
modes_m(ID_xo,:) = modes_xo;
modes = MS_c*modes_m + (lambda_y*Ams.xy + Ams.x)*modes_m*diag(lambda_x);

