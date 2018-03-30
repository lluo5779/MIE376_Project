[weigths, optimal_value] = Solver[f1, f2];

f1 = [22 13 9 17 8 23];

% shortage penalty and surplus penalty
f2 = [16 19 14 28 21 33]./3;

% capacity constraint 
A = [% initial deterministic variables
    ones(1,3) zeros(1, 24-3);
    zeros(1,3) ones(1,3) zeros(1, 24-6)]

b = [%initial constraint 65 80];

% three scenaros 
Aeq = [repmat(eye(3),3,2) kron(eye(n),[1, -1])] % n is the number of equality constraint

beq = [%xi vector corresponding to every variables]

[x, fval] = linprog([f1 f2 f2 f2], A, b, Aeq, beq, zeros(24,1), [])