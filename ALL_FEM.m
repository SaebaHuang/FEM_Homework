%% Problem Information
f = @(x, y) 2 .* sin(x) .* sin(y);
G = [0, pi; 0, pi]; % Rectangle area
n = 8;

%% Make Nodes 1
[nodes_1, node_idxs_1, s_1] = MakeNodes(G, 2*n);

%% Make Units 1
units_1 = MakeUnits_1(2*n);

%% FEM1
u_FEM1 = FEM1(nodes_1, node_idxs_1, units_1, 2*n, s_1, f);

%% Make Nodes 2
[nodes_2, node_idxs_2, s_2] = MakeNodes(G, 2 * n);

%% Make Units 2
units_2 = MakeUnits_2(n);

%% FEM2
u_FEM2 = FEM2(nodes_2, node_idxs_2, units_2, n, 4 * s_2, f);
