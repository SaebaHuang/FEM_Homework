%% Triangle Version FEM-2 with degree aligned first bound condition
function u_FEM2 = FEM2(nodes, node_idxs, units, n, s, f)
A = zeros((2*n - 1)^2);
b = zeros((2*n - 1)^2, 1);
% Cal A and b
for i = 1:numel(units)
    triangle = units{i, 1};
    
    p1 = triangle(1, :);
    p2 = triangle(2, :);
    p3 = triangle(3, :);
    p4 = triangle(4, :); 
    p5 = triangle(5, :);
    p6 = triangle(6, :);
    
    idx_1 = node_idxs{p1(1), p1(2)};
    idx_2 = node_idxs{p2(1), p2(2)};
    idx_3 = node_idxs{p3(1), p3(2)};
    idx_4 = node_idxs{p4(1), p4(2)};
    idx_5 = node_idxs{p5(1), p5(2)};
    idx_6 = node_idxs{p6(1), p6(2)};
    idxs = [idx_1, idx_2, idx_3, idx_4, idx_5, idx_6];
    
    cord_1 = nodes{p1(1), p1(2)};
    cord_2 = nodes{p2(1), p2(2)};
    cord_3 = nodes{p3(1), p3(2)};
    
    a1 = det([cord_2; cord_3]); b1 = -det([1, cord_2(2); 1, cord_3(2)]); c1 = det([1, cord_2(1); 1, cord_3(1)]);
    a2 = -det([cord_1; cord_3]); b2 = det([1, cord_1(2); 1, cord_3(2)]); c2 = -det([1, cord_1(1); 1, cord_3(1)]);
    a3 = det([cord_1; cord_2]); b3 = -det([1, cord_1(2); 1, cord_2(2)]); c3 = det([1, cord_1(1); 1, cord_2(1)]);
    
    syms x y
    L1 = (a1 + b1*x + c1*y) / (2 * s);
    L2 = (a2 + b2*x + c2*y) / (2 * s);
    L3 = (a3 + b3*x + c3*y) / (2 * s);
    phi_1 = L1*(2*L1-1);
    phi_2 = L2*(2*L2-1);
    phi_3 = L3*(2*L3-1);
    phi_4 = 4*L2*L3;
    phi_5 = 4*L1*L3;
    phi_6 = 4*L1*L2;
    
    phi_list = {phi_1, phi_2, phi_3, phi_4, phi_5, phi_6};
    
     % Integral Area
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
    
    for m_tmp = 1:6
        for n_tmp = 1:6
            if ~(idxs(m_tmp) == -1 || idxs(n_tmp) == -1)
                phi_m = phi_list{m_tmp};
                phi_n = phi_list{n_tmp};
                phi_xx = diff(phi_m, 'x') * diff(phi_n, 'x');
                phi_yy = diff(phi_m, 'y') * diff(phi_n, 'y');
                int_func_sym = phi_xx + phi_yy + x*y;
                int_func_tmp = matlabFunction(int_func_sym);
                int_func = @(x,y) int_func_tmp(x, y) - x .* y;
                A(idxs(m_tmp), idxs(n_tmp)) = A(idxs(m_tmp), idxs(n_tmp)) + integral2(int_func, cord_1(1), cord_2(1), y_min, y_max);
            end
        end
    end
    
    phi_i_tmp = phi_1 + x * y;
    phi_i = matlabFunction(phi_i_tmp);
    int_1 = @(x, y) f(x, y) .* (phi_i(x, y) - x .* y);
    
    phi_i_tmp = phi_2 + x * y;
    phi_i = matlabFunction(phi_i_tmp);
    int_2 = @(x, y) f(x, y) .* (phi_i(x, y) - x .* y);
    
    phi_i_tmp = phi_3 + x * y;
    phi_i = matlabFunction(phi_i_tmp);
    int_3 = @(x, y) f(x, y) .* (phi_i(x, y) - x .* y);
    
    phi_i_tmp = phi_4 + x * y;
    phi_i = matlabFunction(phi_i_tmp);
    int_4 = @(x, y) f(x, y) .* (phi_i(x, y) - x .* y);
    
    phi_i_tmp = phi_5 + x * y;
    phi_i = matlabFunction(phi_i_tmp);
    int_5 = @(x, y) f(x, y) .* (phi_i(x, y) - x .* y);
    
    phi_i_tmp = phi_6 + x * y;
    phi_i = matlabFunction(phi_i_tmp);
    int_6 = @(x, y) f(x, y) .* (phi_i(x, y) - x .* y);
    
    if ~(idx_1 == -1)
        b(idx_1) = b(idx_1) + integral2(int_1, cord_1(1), cord_2(1), y_min, y_max);
    end
    if ~(idx_2 == -1)
        b(idx_2) = b(idx_2) + integral2(int_2, cord_1(1), cord_2(1), y_min, y_max);
    end
    if ~(idx_3 == -1)
        b(idx_3) = b(idx_3) + integral2(int_3, cord_1(1), cord_2(1), y_min, y_max);
    end
    if ~(idx_4 == -1)
        b(idx_4) = b(idx_4) + integral2(int_4, cord_1(1), cord_2(1), y_min, y_max);
    end
    if ~(idx_5 == -1)
        b(idx_5) = b(idx_5) + integral2(int_5, cord_1(1), cord_2(1), y_min, y_max);
    end
    if ~(idx_6 == -1)
        b(idx_6) = b(idx_6) + integral2(int_6, cord_1(1), cord_2(1), y_min, y_max);
    end
end

u = A \ b;
u_FEM2 = zeros(2*n+1);
for i = 1:2*n+1
    for j = 1:2*n+1
        idx = node_idxs{i, j};
        if ~(idx == -1)
            u_FEM2(i, j) = u(idx);
        end
    end
end
end