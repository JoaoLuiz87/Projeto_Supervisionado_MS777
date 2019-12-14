function [F] = broy(x)
n = length(x);
F = zeros(n,1);

F(1) = (3-2*x(1))*x(1) - 2*x(2) + 1;

for i = 2:n-1
    F(i) = (3-2*x(i))*x(i) - x(i-1) - 2*x(i+1) + 1;
end

F(n) = (3-2*x(n))*x(n) - x(n-1)  + 1 ;
