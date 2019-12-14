% Programa para testar e gerar as soluções do problema de Convecção-Difusão dado por
% -laplaciano(u) + lambda*u*[u_s + u_t] - f(s,t) = 0
% onde u = u(s,t)  e u = 0 no bordo do quadrado unitário.

%Aqui utilizaremos o chute inicial x0 = (0,0,...,0) e usaremos as
%diferentes escolhas para o termo forçante do Newton-GMRES.

%Saída:
%t = tempo de execução
%iter_in = nº de iterações internas
%iter_ex = nº de iterações externas
%==================================================================================

close all
tic
clear all
L = 63;
h = 1/(L+1);
hh = h*h;
n = L*L;
c1 = 4./hh;
c2 = -1./hh;
xsol = zeros(n,1);

%Parâmetro lambda-----
lambda = 50;
%---------------------

%Chute inicial------
x = zeros(n,1);
%-------------------

u = zeros(L,L);
k = 1;
i = 1;
[X,Y] = meshgrid(h:h:(1-h)); 

%Primeiramente vamos gerar a solução exata para determinar f(s,t)

for s = h : h: (1-h)
    for t = h : h : (1-h)
        %xsol(k) = 10*t*s*(1-t)*(1-s)*exp(s^4.5);
        %xsol(k)=1000*t*s*(1-t)*(1-s)*(s-0.5)*(t-0.5)*exp(s^4.5);
        xsol(k) = (s*2 - s^3)*(sin(3*pi*t));
        k = k + 1;
    end
end

indep = zeros(n,1);
indep = F_ConvDif(xsol , L, indep, lambda);

%Agora aplicamos o método de Newton Inexato a partir do chute inicial
f = @F_ConvDif;
j = @J_ConvDif;
F = f(x, L, indep, lambda);
J = j(x, L, lambda);

F_norm = norm(F)
normvec = [F_norm];
%Inicializando os contadores de iterações externas e internas
k=0;
iter_in = 0;
alpha = 10^-4;

%Definindo os valores para os termos forçantes de acordo com as escolhas
eta_k = 0.01;

while k < 300 

%Plotando a superfície gerada pela solução exata
if k>0 && k<10
    subplot(3,3,k)
    a = 1;
    for i = 1:L
        for r= 1:L
            u(i,r) = x(a);
            a = a+1;
        end
    end
    surf(X, Y, u)
    title(['Iteração ',num2str(k)]);
end    
    %Verificando a tolerância do método
    if F_norm < 10^-4
        break
    end
   
    %Verificando se o termo forçante está próximo da resolução do problema
    if eta_k < 2*(10^-4)
        eta_k = 0.8*(10^-4)*norm(f(x, L, indep, lambda));
    end
    
    eta_k = 1./(2^(k+1));
    
    %Encontrando a direção de descida e calculando o número de iterações
    [s, ~, ~, iter] = gmres(J, -F, 30, eta_k);  
    iter_in = iter_in + iter(1)*iter(2);   
    
    %Verificando a condição de decréscimo
    t=1;
    while norm(f( x + (t*s), L, indep, lambda)) >  norm(f(x, L, indep, lambda))*(1 - (t*alpha)) 
        t = 0.5*t;
    end
    
    %Iteração do método
    xaux = x;
    x = x + t*s;  %Iteração do método utilizando tamanho do passo t
     
   %eta_k = (norm(f(x, L, indep, lambda))/norm(f(xaux, L, indep, lambda)))^1.61803;
    
   
    F = f(x, L, indep, lambda);
    J = j(x, L, lambda);
    F_norm = norm(F)
    normvec = [normvec, F_norm];
    k = k+1;
end
toc
normvec = normvec';
fprintf('O método de Newton Inexato convergiu com %d iterações externas e \n%d iterações internas, para uma solução tal que a norma de F é %d\n', k, iter_in,F_norm);
