function [u_UC, lambdas, phases] = compute_particular_solution(master_slave, D_UC, V_UC, forc_loc, ...
                                                N_UC_x, N_UC_y, mesh_data, F_UC_FOM, scalar)

%% construct master-slave
MS_c = master_slave.Ams.c;
MS_x = master_slave.Ams.x;
MS_y = master_slave.Ams.y;
MS_xy = master_slave.Ams.xy;

%% get the solution
force_UC_x = forc_loc(1);
force_UC_y = forc_loc(2);
force_UC_node = forc_loc(3);
if force_UC_node == -1
    f = F_UC_FOM;
else
    force_UC_dof = nonzeros(mesh_data.nv_data(force_UC_node,:));
    f = zeros(size(V_UC,1), 1);
    f(force_UC_dof(end)) = scalar;
end
f_UC = V_UC'*f;
u_UC = [];
lambdas = [];
phases = [];
for iN = 0:2*N_UC_x-1
    for jN = 0:N_UC_y-1
        lambda_x = exp(2*1i*pi*iN/2/N_UC_x);
        lambda_y = exp(2*1i*pi*jN/N_UC_y);
        force_ij = (1/lambda_x)^(force_UC_x-1) *(1/lambda_y)^(force_UC_y-1)* f_UC/(2*N_UC_x*N_UC_y)...
                 - (1/lambda_x)^(N_UC_x+force_UC_x-1) *(1/lambda_y)^(force_UC_y-1)* f_UC/(2*N_UC_x*N_UC_y);
%                  - (1/lambda_x)^(force_UC_x-1) *(1/lambda_y)^(N_UC_y+force_UC_y-1)* f_UC/(4*N_UC_x*N_UC_y)...
%                  + (1/lambda_x)^(N_UC_x+force_UC_x-1) *(1/lambda_y)^(N_UC_y+force_UC_y-1)* f_UC/(4*N_UC_x*N_UC_y);

        if norm(force_ij)>eps
            A_R = MS_c + lambda_x*MS_x + lambda_y*MS_y + lambda_x*lambda_y*MS_xy;
            A_L = MS_c + 1/lambda_x*MS_x + 1/lambda_y*MS_y + 1/(lambda_x*lambda_y)*MS_xy;
            u_UC = [u_UC, (A_R* ((A_L.'*D_UC*A_R)\(A_L.'*force_ij)))];
            lambdas = [lambdas; [lambda_x, lambda_y]];
            phases = [phases; [0,0]];
        end
    end
end
