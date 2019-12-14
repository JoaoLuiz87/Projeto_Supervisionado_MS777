%Esta fun��o calcula o vetor F para o problema da Convec��o-Difus�o dado
%por -laplaciano(u) + lambda*u*( u_s + u_t) = f(s,t).
%Primeiramente iremos calcular o termo das derivadas parciais de primeira
%ordem e, em seguida, iremos somar com o laplaciano calculado por outra
%fun��o.

function F = F_ConvDif(x, L, indep, lambda)

h = 1/(L+1);
n = L*L;
hh = h*h;
k = 1;
ndl = n - L;

%Primeira equa��o: i = 1, j = 1

F(1) = x(1)*( x(2) + x(1+L) );

%Equa��es 2 a L-1, ou seja, 1 < i < L-1  e  j = 1 [ x(k-L) = 0 ]:
for k = 2:L-1
   F(k) = x(k)*( x(k+1) - x(k-1) + x(k+L));
end

% L-�sima equa��o: i = L  e j = 1  [x(L+1) = 0]:
F(L) = x(L)*( -x(L-1) + x(L+L) );

k = L;
j = 1;

while k < ndl
    k = k+1;
%   j = j+1;
    
    % Equa��o com i = 1  e  1 < j < l  [ x(k-1) = 0 ]:

    F(k) = x(k)*( x(k+1) + x(k+L) - x(k-L));

    % Equa��es com 1 < i,j < L  [Nenhum termo no bordo]
    for i = 2:L-1
        k = k+1;
        F(k) = x(k)*( x(k+1) - x(k-1) + x(k+L) - x(k-L));
    end 
    k = k+1;
    %Equa��o com i = L e 1 < j < L  [Borda direita]
    F(k) = x(k)*( -x(k-1) + x(k+L) - x(k-L) );
end
k = k+1; %Indo para a ultima fileira da discretiza��o

%Equa��o i = 1  e  j = L [x(k-1) = x(k+L) = 0]
F(k) = x(k)*( x(k+1) - x(k-L));

%Equa��es 1 < i < L  e  j = L [x(k+L) = 0]
for i = 2 : L-1
   k = k + 1;
   F(k)= x(k)*( x(k+1) - x(k-1) - x(k-L) );
end

%�ltima equa��o da discretiza��o i = L  e  j = L:
k = k+1;
F(k) = x(k)*( - x(k-1) - x(k-L));

F = lambda*(1/(2*h))*F;

%Avaliando o laplaciano
lap = feval('laplaciano', x, L);

%Definindo a equa��o do PVC
F = F' + lap';
F = F - indep;
