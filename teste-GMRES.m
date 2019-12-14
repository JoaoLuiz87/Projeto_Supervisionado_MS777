% -----------------------------------------------
% |      TESTES COMPUTACIONAIS COM GMRES        |
% -----------------------------------------------

dim = [20,50,80,100,300,500]';     %Vetor com a dimens�o das matrizes de teste
n_teste = length(dim);    %Tamanho do vetor dimens�o

%Vetores auxiliares na montagem da tabela
flag_vec = zeros(n_teste,1);
iter_vec = zeros(n_teste,1);
error_vec = zeros(n_teste,1);
cond_vec = zeros(n_teste,1);

%Iterando a cada dimens�o
for k=1:n_teste
    
    %Definindo a matriz do sistema linear e o vetor independente
    
    e =  round(10*rand(dim(k))-5);
    A = spdiags(e);
    
    %A = hilb(dim(k)); %Valores no intervalo [-50,50]
    
    x_sol = round(10*rand(dim(k),1)-5); %Valores no intervalo[-5,5]
    b = A*x_sol;
    
    %Inicializando o m�todo GMRES
    [x, flag, relres, iter, resvec] = gmres(A,b,[],[],dim(k));
    
    %C�lculo do erro entre a solu��o original e a calculada pelo GMRES
    erro = norm(x - x_sol);
    
    %Guardando os valores do m�todo para cria��o da tabela
    flag_vec(k,1) = flag;
    iter_vec(k,1) = iter(2);
    error_vec(k,1) = erro;
    cond_vec(k,1) = cond(A);
    
    %Vetor utilizado para plotar o gr�fico das itera��es (come�ando do 0)
    aux = 0:iter(2);
    resvec = resvec/norm(b);
     length(resvec)
     length(aux)
    
    %Plot do gr�fico da norma-2 res�duo a cada itera��o
    plot(aux, resvec,'-o', 'LineWidth','2');
    hold on
        
 
end
xlabel('N�mero de Itera��es')
ylabel('Norma Relativa do Res�duo')
legend('n=20','n=50','n=80','n=100', 'n=300', 'Location', 'Northeast')
hold off

%Gera��o da tabela
%T = table(dim,flag_vec,cond_vec, error_vec, iter_vec)

