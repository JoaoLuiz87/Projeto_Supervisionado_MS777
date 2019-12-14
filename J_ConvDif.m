%Calcula a matriz jacobiana para o termo não linear da equação da
%Difusão-Convecção.

function J = J_ConvDif(x, L, lambda)


h = 1/(L+1);
hh = 2*h; 
n = L*L;
nL = n - L;
A = zeros(n,n);

%Derivadas da primeir%a equação

A(1,1) = (x(2) + x(L+1));
A(1,2) = x(1);
A(1,1+L) = x(1);

%Derivadas parcias das equações 2 a L-1:
for k = 2:L-1
    A(k,k-1) = -x(k);
    A(k,k) = x(k+1) - x(k-1) + x(k+L) ;
    A(k,k+1) = x(k);
    A(k,k+L) = x(k);
end

%Derivafas parciais da equação L:

A(L,L-1) = -x(L);
A(L,L) = -x(L-1) + x(L+L);
A(L,L+L) = x(L);

%	Derivadas parciais das equações (L+1) até ( n-L )
k = L;

while( k < nL)
    k = k+1;
    
    %Derivadas parciais das equações do tipo i = 1  e  1 < j < L:
    A(k,k-L) = -x(k);
    A(k,k) = ( x(k+1) + x(k+L) - x(k-L) );
    A(k,k+1) = x(k);
    A(k,k+L) = x(k);

    %Derivadas parciais das equações do tipo 1 < i < L  e  1 < j < L
	for i = 2:L-1
        k = k+1;
        A(k,k-L) = -x(k); 
        A(k,k-1) = -x(k);
        A(k,k) = x(k+1) - x(k-1) + x(k+L) - x(k-L) ;
        A(k,k+1) = x(k);
        A(k,k+L) = x(k);
    end
    
    %Derivadas parciais da equação do tipo i = L e 1 < j < L
    k = k+1;
    A(k,k-L) = -x(k);
 	A(k,k-1)  = -x(k);
    A(k,k) = -x(k-1) - x(k-L) + x(k+L);
	A(k,k+L) = x(k); 
end

%Derivadas parciais da equação do tipo i = 1 e j = L
    k = k+1;
	A(k,k-L) = -x(k);
	A(k,k) = -x(k-L) + x(k+1);
	A(k,k+1) = x(k);
    
%Derivadas parciais das equações do tipo 1 < i < L e j = L
%
    for i = 2: L-1
	   k = k+1;
	   A(k,k-L) = -x(k);
	   A(k,k-1) = -x(k);
	   A(k,k) = x(k+1) - x(k-1) - x(k-L);
	   A(k,k+1) = x(k);
    end
    
    %Derivadas parciais da equacao tipo i = L e j = L
    k = k + 1;
    A(k,k-L) = -x(k);
    A(k,k-1) = -x(k) ;
    A(k,k) = -x(k-1) - x(k-L);
    
    J = A;
    J = (lambda/(hh))*J;
    
    J = J + J_laplaciano(L);