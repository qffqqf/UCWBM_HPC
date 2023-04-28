function mesh = convert_mesh(mesh_data, U, n_dofnodes)

mesh.nodes = mesh_data.nd_data + [U(1:n_dofnodes:end),U(2:n_dofnodes:end),U(3:n_dofnodes:end)+1e-10*rand(size(mesh_data.nd_data,1),1)];
mesh.field1 = U(3:n_dofnodes:end);
mesh.field2 = sqrt(U(1:n_dofnodes:end).^2 + U(2:n_dofnodes:end).^2 + U(3:n_dofnodes:end).^2);
mesh.element_type = mesh_data.element_type;
switch mesh_data.element_type
    case "C3D10"
        reorder = [1,3,6,10,2,5,4,7,8,9];
        mesh.connectivity = mesh_data.en_data(:, reorder);
    case "C3D4"
        mesh.connectivity = mesh_data.en_data;
    case "C3D27"
        reorder = [1,7,9,3,19,25,27,21,...
                   4,8,6,2,22,26,24,20,...
                   10,16,18,12];
        mesh.connectivity = mesh_data.en_data(:, reorder);
    case "S9" 
        reorder = [1,3,9,7,2,6,8,4];
        mesh.connectivity = mesh_data.en_data(:, reorder);
end

