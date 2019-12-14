function [J] = Jpvc(y, L)
h = 1/(L+1);
n = L;
J = zeros(n,n);

J(1,1) = -h^2*(sin(y(1))+ h + y(1)*cos(y(1))) - 2;
J(1,2) = 1;

for k = 2:n-1
    J(k,k-1) = 1;
    J(k,k) = -h^2*(sin(y(k))+ k*h + y(k)*cos(y(k))) - 2;
    J(k,k+1) = 1;
end

J(n,n-1) = 1;
J(n,n) = -h^2*(sin(y(n))+ n*h + y(n)*cos(y(n))) - 2;