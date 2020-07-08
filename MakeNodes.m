%% Make Mesh for space G
function [nodes, node_idxs, s] = MakeNodes(G, n)
nodes = cell(n + 1, n + 1);
node_idxs = cell(n + 1, n + 1);
count = 1;
x_st = G(1, 1);
x_end = G(1, 2);
y_st  = G(2, 1);
y_end = G(2, 2);
s = (x_end - x_st) * (y_end - y_st) / (2 * n * n); 
for i = 1:n + 1 
    for j = 1:n + 1
        nodes{j, i} = [x_st + (x_end - x_st) * (i - 1) / n, y_st + (y_end - y_st) * (j - 1) / n];
        if ~(i == 1 || i == n + 1 || j == 1 || j == n + 1)
            node_idxs{j, i} = count;
            count = count + 1;
        else
            node_idxs{j, i} = -1;
        end
    end
end
end