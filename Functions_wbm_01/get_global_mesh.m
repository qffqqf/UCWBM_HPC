function mesh_data_global = get_global_mesh(N_UC_x, N_UC_y, mesh_data, coords_disc)

%% nodal coords
mesh_data_global.nd_data = coords_disc.all(:,2:4);
mesh_data_global.element_type = mesh_data.element_type;
nNode = size(mesh_data_global.nd_data,1);
mesh_data_global.nv_data = reshape([1:3*nNode], [3, nNode])';
mesh_data_global.nDOF = 3*nNode;

%% elements
ele_disc = cell(1,N_UC_x*N_UC_y);
n_nd_UC = size(mesh_data.nd_data,1);
for i = 1:N_UC_x*N_UC_y
    ele_disc{i} = mesh_data.en_data + (i-1)*n_nd_UC;
end
mesh_data_global.en_data = cat(1, ele_disc{:});

%% surf elements
% ele_disc = cell(1,N_UC_x*N_UC_y);
% for i = 1:N_UC_x*N_UC_y
%     ele_disc{i} = mesh_data.en_data_surf + (i-1)*n_nd_UC;
% end
% mesh_data_global.en_data_surf = cat(1, ele_disc{:});

%% set type
% mesh_data_global.element_type = mesh_data.element_type;