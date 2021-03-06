function [x, fval] = Solver(xi, b, G)

%number of scenarios
[r, c] = size(xi);

%objective function with all y and w at the end with probability 1/r 
obj = [zeros(1, 20), repmat([2, -5], 1, r)]./r;

%first constraint
FC = [ones(1,20), zeros(1,r*2)];

% %second constraint 
% SC = zeros(r, 20+22*r);
% 
% for i = 1:r
%     SC(i,1:20) = xi(i,:);
%     SC(i,20*i+1:20+20*i) = -1*ones(1,20);
% end

%third constraint
TC = zeros(r, 20+2*r);

for i = 1:r
    TC(i,(1:20)) = xi(i,:);
    TC(i,19+2*i:19+2*i+1) = [-1,1];
end

[x, fval] = linprog(-obj, [], [], [FC;TC], [b; repmat(G,[r,1])], zeros(20+2*r,1),[]);

% MULTI-STAGE STOCHASTIC PROGRAM
% 
% %number of scenarios
% [r, c] = size(xi);
% 
% %objective function with all y and w at the end with probability 1/r 
% obj = [zeros(1, 20+r*20), repmat([2, -5], 1, r)]./r;
% 
% %first constraint
% FC = [ones(1,20), zeros(1,r*22)];
% 
% %second constraint 
% SC = zeros(r, 20+22*r);
% 
% for i = 1:r
%     SC(i,1:20) = xi(i,:);
%     SC(i,20*i+1:20+20*i) = -1*ones(1,20);
% end
% 
% %third constraint
% TC = zeros(r, 20+22*r);
% 
% for i = 1:r
%     TC(i,(20*i+1):(20+20*i)) = xi(i,:);
%     TC(i,19+20*r+2*i:19+20*r+2*i+1) = [-1,1];
% end
% 
% [x, fval] = linprog(-obj, [], [], [FC;SC;TC], [b; zeros(r,1); repmat(G,[r,1])], zeros(20+22*r,1),[]);
