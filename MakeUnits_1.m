%% Make Units
function units = MakeUnits_1(n)
units = cell(2 * n * n, 1);
tri_id = 1;
for i = 1:n
    for j = 1:n
        p1 = [i, j];
        p2 = [i, j + 1];
        p3 = [i + 1, j];
        p4 = [i + 1, j + 1];
        units{tri_id} = [p1; p2; p3];
        units{tri_id + 1} = [p3; p2; p4];
        tri_id = tri_id + 2;
    end
end
end