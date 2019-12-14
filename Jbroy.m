function [J] = Jbroy(x)

n = length(x); 
J = zeros(n,n);

J(1,1) = 3 - 4*x(1);
J(1,2) = -2;

for k = 2:n-1
    J(k,k-1) = -1;
    J(k,k) = 3 - 4*x(k);
    J(k,k+1) = -2;
end

J(n,n-1) = -1;
J(n,n) = 3 - 4*x(n);

