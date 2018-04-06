[x, fval] = Solver(xi, b, G);

%number of scenarios
[r, c] = size(xi);

%objective function with all y and w at the end with probability 1/r 
obj = [zeros(1, 20+r*20), repmat([5, -10], 1, r)]./r;

%first constraint
FC = [ones(1,20), zeros(1,r*23)];

%second constraint 
SC = zeros(3, 20+20*r+3*r);

for i = 1:r
    SC(i,1:20) = xi(i,:);
    SC(i,20*i+1:20+20*i) = ones(1,20);
end

%third constraint
TC = zeros(3, 20+20*r+3*r);

for i = 1:r
    TC(i,20*i+1:20+20*i) = xi(i,:);
    TC(i,19+20*r+2*i:19+20*r+2*i+1) = [-1,1];
end

[x, fval] = linprog(-obj, [], [], [FC;SC;TC], [b; zeros(r,1); repmat(5,[r,1])], zeros(20+23*r,1),[]);
