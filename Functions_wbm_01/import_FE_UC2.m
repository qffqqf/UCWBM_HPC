function [K_UC,M_UC,C_UC,UC_nodes,UC_coordinates,L_UC_x,L_UC_y] = import_FE_UC2(filename, import_option)

load([filename,'.mat']);
M_UC = comsol_matrices.mass;
K_UC = comsol_matrices.stfK;
C_UC = 0* comsol_matrices.mass;

K_UC = (K_UC + K_UC')/2;
M_UC = (M_UC + M_UC')/2;

L_UC_x = max(mesh_data.nd_data(:,1));
L_UC_y = max(mesh_data.nd_data(:,2));

[uc_nodes, ~] = pick_uc_boundary(mesh_data, [L_UC_x, L_UC_y]);
UC_nodes.I = uc_nodes.I';
UC_nodes.L = uc_nodes.L';
UC_nodes.R = uc_nodes.R';
UC_nodes.B = uc_nodes.B';
UC_nodes.T = uc_nodes.T';
UC_nodes.BL = uc_nodes.BL';
UC_nodes.BR = uc_nodes.BR';
UC_nodes.TL = uc_nodes.TL';
UC_nodes.TR = uc_nodes.TR';
UC_coordinates = [[1:size(mesh_data.nd_data,1)]', mesh_data.nd_data];
