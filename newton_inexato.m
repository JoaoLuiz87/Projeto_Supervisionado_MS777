%Fun��o Newton-GMRES.

function [x, normvec, iter_in, k] = newton_inexato(x, func, jfunc, L)
tic
%Recendo os valores iniciais da fun��o e da Jacobiana 
F = func(x, L);
J = jfunc(x, L);

%C�lculo da norma-2 inicial de F 
F_norm = norm(F);
normvec = [F_norm];
%Inicializando os contadores de itera��es externas e internas
k=0;
iter_in = 0;
alpha = 10^-4;

%Definindo os valores para os termos for�antes de acordo com as escolhas
eta_k = 0.1;
while k < 300 
    
    if F_norm < 10^-4
        break
    end
    
    if eta_k < 2*(10^-4)
        eta_k = 0.8*(10^-4)*norm(func(x, L));
    end
    
    [s, ~, ~, iter] = gmres(J, -F, 30, eta_k);  %Encontrando a dire��o de descida
    iter_in = iter_in + iter(1)*iter(2);   % N�mero de itera��es internas
    
    t=1;
    while norm(func( x + (t*s), L)) >  norm(func(x, L))*(1 - (t*alpha))
        t = 0.5*t;
    end
    
    xaux = x;
    x = x + t*s;  %Itera��o do m�todo utilizando tamanho do passo t
    
   %eta_k = 1./(k+2);
   %eta_k = (norm(func(x, L))/norm(func(xaux, L)))^1.61803398875
    
    F = func(x, L);
    J = jfunc(x, L);
    F_norm = norm(F);
    normvec = [normvec, F_norm];
    k = k+1;
end
toc
normvec = normvec';
fprintf('O m�todo de Newton Inexato convergiu com %d itera��es externas e \n%d itera��es internas, para uma solu��o tal que a norma de F � %d\n', k, iter_in,F_norm);

    
    