%Arquivo que avalia a fun��o do problema dado por
% -laplaciano(u) + lambda*exp(u) - f(s,t) = 0, onde f(s,t) = indep

function F = F_Bratu(x, L, indep, lambda)

%Abaliando o Laplaciano
lap = feval('laplaciano', x, L);

%Avaliando a fun��o h
h = -lambda*(exp(x));

%C�lculo da fun��o F
F = lap' + h;
F = F - indep;

end