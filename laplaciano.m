function f = laplaciano(x, L)
h = 1/(L+1);
n = L*L;
ndl = n - L;
k = 1;

%Primeira equa��o i = 1, j = 1
f(1) = -(-4*x(1) + x(2) + x(1+L));

%Equa��es 2 a L-1, ou seja, 1 < i < L-1  e  j = 1 [ x(k-L) = 0 ]:
for k = 2:L-1
   f(k) = -(-4*x(k) + x(k+1) + x(k-1) + x(k+L));
end

% L-�sima equa��o: i = L  e j = 1  [x(L+1) = 0]:
f(L) = -(-4*x(L) + x(L-1) + x(L+L));

k = L;
j = 1;

while k < ndl
    k = k+1;
%   j = j+1;
    
    % Equa��o com i = 1  e  1 < j < l  [ x(k-1) = 0 ]:

    f(k) = -(-4*x(k) + x(k+1) + x(k+L) + x(k-L));

    % Equa��es com 1 < i,j < L  [Nenhum termo no bordo]
    for i = 2:L-1
        k = k+1;
        f(k) = -(-4*x(k) + x(k+1) + x(k-1) + x(k+L) + x(k-L));
    end 
    k = k+1;
    %Equa��o com i = L e 1 < j < L  [Borda direita]
    f(k) = -(-4*x(k) + x(k-1) + x(k+L) + x(k-L));
end
k = k+1; %Indo para a ultima fileira da discretiza��o

%Equa��o i = 1  e  j = L [x(k-1) = x(k+L) = 0]
f(k) = -(-4*x(k) + x(k+1) + x(k-L));

%Equa��es 1 < i < L  e  j = L [x(k+L) = 0]
for i = 2 : L-1
   k = k + 1;
   f(k)=-(-4*x(k)+x(k+1)+x(k-1)+x(k-L));
end

%�ltima equa��o da discretiza��o i = L  e  j = L:
k = k+1;
f(k) = -(-4*x(k) + x(k-1) + x(k-L));

f = (1/(h*h))*f;