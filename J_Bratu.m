%Calcula o Jacobiano da função F do problema de Bratu

function J = J_Bratu(x, L, lambda)

h = 1/(L+1);
hh = h*h;
n = L*L;
c1 = 4./hh;
c2 = -1./hh;
nL = n - L;

%Primeiramente vamos determinar o Jacobiano do Laplaciano
A = zeros(n,n);

%Derivadas parciais da equação 1
A(1,1) = c1;
A(1,2) = c2;
A(1,1+L) = c2;

%Derivadas parcias das equações 2 a L-1:
for k = 2:L-1
    A(k,k-1) = c2;
    A(k,k) = c1;
    A(k,k+1) = c2;
    A(k,k+L) = c2;
end

%Derivada parcial para a equação L:
A(L,L-1) = c2;
A(L,L) = c1;
A(L,L+L) = c2;

%	Derivadas parciais das equações (L+1) até ( n-L )
k = L;

while( k < nL)
    k = k+1;
    
    %Derivadas parciais das equações do tipo i = 1  e  1 < j < L:
    A(k,k-L) = c2;
    A(k,k) = c1;
    A(k,k+1) = c2;
    A(k,k+L) = c2;

    %Derivadas parciais das equações do tipo 1 < i < L  e  1 < j < L
	for i = 2:L-1
        k = k+1;
        A(k,k-L) = c2; 
        A(k,k-1) = c2;
        A(k,k) = c1;
        A(k,k+1) = c2;
        A(k,k+L) = c2;
    end
    
    %Derivadas parciais da equação do tipo i = L e 1 < j < L
    k = k+1;
    A(k,k-L) = c2;
 	A(k,k-1)  = c2;
    A(k,k) = c1;
	A(k,k+L) = c2; 
end

%Derivadas parciais da equação do tipo i = 1 e j = L
    k = k+1;
	A(k,k-L) = c2;
	A(k,k) = c1;
	A(k,k+1) = c2;
    
%Derivadas parciais das equações do tipo 1 < i < L e j = L
%
	for i = 2: L-1
	   k = k+1;
	   A(k,k-L) = c2;
	   A(k,k-1) = c2;
	   A(k,k) = c1;
	   A(k,k+1) = c2;
    end
    
%Derivadas parciais da equacao tipo i = L e j = L
    k = k + 1;
    A(k,k-L) = c2;
    A(k,k-1) = c2;
    A(k,k) = c1;

%Gerando a matriz do termo não linear
B = zeros(n,n);

for i=1:n
    B(i,i) = -lambda*(exp(x(i)));
end

%Matriz jacobiana de dimensão n x n do sistema L x L.
J = A + B;

    