%% Triangle Version FEM-1 with degree aligned first bound condition
function u_FEM1 = FEM1(nodes, node_idxs, units, n, s, f)
A = zeros((n - 1)^2);
b = zeros((n - 1)^2, 1);
% Cal A and b
for i = 1:numel(units)
    triangle = units{i, 1};
    p1 = triangle(1, :);
    p2 = triangle(2, :);
    p3 = triangle(3, :);
    idx_1 = node_idxs{p1(1), p1(2)};
    idx_2 = node_idxs{p2(1), p2(2)};
    idx_3 = node_idxs{p3(1), p3(2)};
    cord_1 = nodes{p1(1), p1(2)};
    cord_2 = nodes{p2(1), p2(2)};
    cord_3 = nodes{p3(1), p3(2)};
    a1 = det([cord_2; cord_3]); b1 = -det([1, cord_2(2); 1, cord_3(2)]); c1 = det([1, cord_2(1); 1, cord_3(1)]);
    a2 = -det([cord_1; cord_3]); b2 = det([1, cord_1(2); 1, cord_3(2)]); c2 = -det([1, cord_1(1); 1, cord_3(1)]);
    a3 = det([cord_1; cord_2]); b3 = -det([1, cord_1(2); 1, cord_2(2)]); c3 = det([1, cord_1(1); 1, cord_2(1)]);
    a_int_1_1 = (b1^2 + c1^2) / (4 * s);
    a_int_1_2 = (b1*b2 + c1*c2) / (4 * s);
    a_int_1_3 = (b1*b3 + c1*c3) / (4 * s);
    a_int_2_2 = (b2^2 + c2^2) / (4 * s);
    a_int_2_3 = (b2*b3 + c2*c3) / (4 * s);
    a_int_3_3 = (b3*b3 + c3*c3) / (4 * s);
    if ~(idx_1 == -1 || idx_2 == -1)
        A(idx_1, idx_2) = A(idx_1, idx_2) + a_int_1_2;
        A(idx_2, idx_1) = A(idx_2, idx_1) + a_int_1_2;
    end
    if ~(idx_1 == -1 || idx_3 == -1)
        A(idx_1, idx_3) = A(idx_1, idx_3) + a_int_1_3;
        A(idx_3, idx_1) = A(idx_3, idx_1) + a_int_1_3;
    end
    if ~(idx_3 == -1 || idx_2 == -1)
        A(idx_2, idx_3) = A(idx_2, idx_3) + a_int_2_3;
        A(idx_3, idx_2) = A(idx_3, idx_2) + a_int_2_3;
    end
    if ~(idx_1 == -1)
        A(idx_1, idx_1) = A(idx_1, idx_1) + a_int_1_1;
    end
    if ~(idx_2 == -1)
        A(idx_2, idx_2) = A(idx_2, idx_2) + a_int_2_2;
    end
    if ~(idx_3 == -1)
        A(idx_3, idx_3) = A(idx_3, idx_3) + a_int_3_3;
    end
    int_1 = @(x, y) f(x, y) .* (a1 + b1 .* x + c1 .* y);
    int_2 = @(x, y) f(x, y) .* (a2 + b2 .* x + c2 .* y);
    int_3 = @(x, y) f(x, y) .* (a3 + b3 .* x + c3 .* y);
    if cord_1(2) == cord_2(2)
        y_min = cord_1(2);
        k_max = (cord_3(2) - cord_2(2)) / (cord_3(1) - cord_2(1));
        y_max = @(x) k_max * (x - cord_1(1)) + cord_3(2);
    end
    if cord_1(2) == cord_3(2)
        k_min = (cord_2(2) - cord_1(2)) / (cord_2(1) - cord_1(1));
        y_min = @(x) k_min * (x - cord_1(1)) + cord_1(2);
        y_max = cord_1(2);
    end
    if ~(idx_1 == -1)
        b(idx_1) = b(idx_1) + integral2(int_1, cord_1(1), cord_2(1), y_min, y_max) / (2 * s);
    end
    if ~(idx_2 == -1)
        b(idx_2) = b(idx_2) + integral2(int_2, cord_1(1), cord_2(1), y_min, y_max) / (2 * s);
    end
    if ~(idx_3 == -1)
        b(idx_3) = b(idx_3) + integral2(int_3, cord_1(1), cord_2(1), y_min, y_max) / (2 * s);
    end
end
u = A \ b;
u_FEM1 = zeros(n+1);
for i = 1:n+1
    for j = 1:n+1
        idx = node_idxs{i, j};
        if ~(idx == -1)
            u_FEM1(i, j) = u(idx);
        end
    end
end
end

