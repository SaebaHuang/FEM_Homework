%% Make Units
function units = MakeUnits_2(n)
units = cell(2 * n * n, 1);
tri_id = 1;
for i = 1:2:2*n-1
    for j = 1:2:2*n-1
        p1 = [i, j];
        p2 = [i, j+1];
        p3 = [i, j+2];
        p4 = [i+1, j];
        p5 = [i+1, j+1];
        p6 = [i+1, j+2];
        p7 = [i+2, j];
        p8 = [i+2, j+1];
        p9 = [i+2, j+2];
        units{tri_id} = [p1; p3; p7; p5; p4; p2];
        units{tri_id + 1} = [p7; p3; p9; p6; p8; p5];
        tri_id = tri_id + 2;
    end
end
end