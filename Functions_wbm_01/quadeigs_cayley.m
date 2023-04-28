function [modes, eigvals] = quadeigs_cayley(varargin)

Cm1 = varargin{1};
C0 = varargin{2};
C1 = varargin{3};
%% Linearization
A = [-Cm1, -C0; 0*C0, eye(size(C0,1))];
C = [0*C0, C1; eye(size(C0,1)), 0*C0];
%% Transformed EVP
[modes_ln, cosi] = eig(A+C,A-C,'vector');
cosi_inf = find(isinf(cosi));
ind = abs(cosi-1)>1e-20;
eigvals = (cosi+1)./(cosi-1).*ind + 1e20.*(~ind);
eigvals(cosi_inf) = 1;
[~,sort_id] = sort(abs(log(abs(eigvals))));

%% For each eigenvalue, extract the eigenvector from whichever portion
% of the big eigenvector matrix X gives the smallest normalized residual.
n = size(C0,1);
p = 2;
np = n*p;
V = zeros(n,p);
for j = 1:np
   V(:) = modes_ln(:,j);
   R = varargin{p+1};
   if ~isinf(eigvals(j))
       for k = p:-1:1
           R = varargin{k} + eigvals(j)*R;
       end
   end
   R = R*V;
   res = sum(abs(R))./ sum(abs(V));  % Normalized residuals.
   [~,ind] = min(res);
   modes_ln(1:n,j) = V(:,ind)/norm(V(:,ind),2);  % Eigenvector with unit 2-norm.
end
modes = modes_ln(1:n,:);

%% sort
eigvals = eigvals(sort_id);
modes = modes(:,sort_id);