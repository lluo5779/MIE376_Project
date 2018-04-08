function [x_optimal, fval] = Solver_mat(xi, initial_funding, G)

    [no_splits, no_assets] = size(xi);

    f = -1*[zeros(1, no_assets*(no_splits+1)) repmat([1 -100], 1, no_splits^2)]/(no_splits^2);
    total_no_var = size(f, 2);
    
    ones_1asset = ones(1, no_assets);
    one_set_of_yw = kron(eye(no_splits), [-1 1]);
    Aeq = [ones_1asset zeros(1,total_no_var - no_assets);
        xi -1*kron(eye(no_splits), ones_1asset) zeros(no_splits, 2*no_splits^2);
        zeros(no_splits^2, no_assets) kron(eye(no_splits), xi) kron(eye(no_splits), one_set_of_yw)];
    beq = [initial_funding; zeros(no_splits, 1); G*ones(no_splits^2, 1)];
    
    [x_optimal, fval] = linprog(f, [], [], Aeq, beq, zeros(1, total_no_var), [])  ;
    
end
