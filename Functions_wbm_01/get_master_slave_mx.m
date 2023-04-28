function master_slave = get_master_slave_mx(UC_dofs)

%% construct master-slave matrices
del_coln = [UC_dofs.R, UC_dofs.T, UC_dofs.BR, UC_dofs.TR, UC_dofs.TL];
master_slave.Ams.c = zeros(UC_dofs.nDOF, UC_dofs.nDOF);
master_slave.Ams.c(UC_dofs.I, UC_dofs.I) = eye(numel(UC_dofs.I));
master_slave.Ams.c(UC_dofs.L, UC_dofs.L) = eye(numel(UC_dofs.L));
master_slave.Ams.c(UC_dofs.B, UC_dofs.B) = eye(numel(UC_dofs.B));
master_slave.Ams.c(UC_dofs.BL, UC_dofs.BL) = eye(numel(UC_dofs.BL));
master_slave.Ams.c(:, del_coln) = [];
master_slave.Ams.x = zeros(UC_dofs.nDOF, UC_dofs.nDOF);
master_slave.Ams.x(UC_dofs.R, UC_dofs.L) = eye(numel(UC_dofs.L));
master_slave.Ams.x(UC_dofs.BR, UC_dofs.BL) = eye(numel(UC_dofs.BL));
master_slave.Ams.x(:, del_coln) = [];
master_slave.Ams.y = zeros(UC_dofs.nDOF, UC_dofs.nDOF);
master_slave.Ams.y(UC_dofs.T, UC_dofs.B) = eye(numel(UC_dofs.B));
master_slave.Ams.y(UC_dofs.TL, UC_dofs.BL) = eye(numel(UC_dofs.BL));
master_slave.Ams.y(:, del_coln) = [];
master_slave.Ams.xy = zeros(UC_dofs.nDOF, UC_dofs.nDOF);
master_slave.Ams.xy(UC_dofs.TR, UC_dofs.BL) = eye(numel(UC_dofs.BL));
master_slave.Ams.xy(:, del_coln) = [];

%% get master ids
master_ids = [1:UC_dofs.nDOF];
master_ids(UC_dofs.I) = 1;
master_ids(UC_dofs.L) = 2;
master_ids(UC_dofs.B) = 3;
master_ids(UC_dofs.BL) = 4;
master_ids(del_coln) = []; 
master_slave.master_dofs.I = find(master_ids==1);
master_slave.master_dofs.L = find(master_ids==2);
master_slave.master_dofs.B = find(master_ids==3);
master_slave.master_dofs.BL = find(master_ids==4);
