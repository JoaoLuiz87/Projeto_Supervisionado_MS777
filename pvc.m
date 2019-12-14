function [F] = pvc(y, L)
h = 1/(L+1);
n = L;
%LEMBRAR QUE Y =(Y1,Y2,...,YN) SÃO AS VARÁVEIS QUE QUEREMOS
F = zeros(n,1);
F(1) = 1 - y(1)*(h^2 * sin(y(1)) + h + 2) + y(2);

for i = 2:n-1
    F(i) = y(i-1) - y(i)*(h^2 * sin(y(i)) + i*h + 2) + y(i+1);
end

F(n) = y(n-1) - y(n)*(h^2 * sin(y(n)) + (n-1)*h + 2) + 5;