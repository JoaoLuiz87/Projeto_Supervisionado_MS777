% -----------------------------------------------
% |      TESTES COMPUTACIONAIS COM GMRES        |
% -----------------------------------------------

dim = [20,50,80,100,300,500]';     %Vetor com a dimensão das matrizes de teste
n_teste = length(dim);    %Tamanho do vetor dimensão

%Vetores auxiliares na montagem da tabela
flag_vec = zeros(n_teste,1);
iter_vec = zeros(n_teste,1);
error_vec = zeros(n_teste,1);
cond_vec = zeros(n_teste,1);

%Iterando a cada dimensão
for k=1:n_teste
    
    %Definindo a matriz do sistema linear e o vetor independente
    
    e =  round(10*rand(dim(k))-5);
    A = spdiags(e);
    
    %A = hilb(dim(k)); %Valores no intervalo [-50,50]
    
    x_sol = round(10*rand(dim(k),1)-5); %Valores no intervalo[-5,5]
    b = A*x_sol;
    
    %Inicializando o método GMRES
    [x, flag, relres, iter, resvec] = gmres(A,b,[],[],dim(k));
    
    %Cálculo do erro entre a solução original e a calculada pelo GMRES
    erro = norm(x - x_sol);
    
    %Guardando os valores do método para criação da tabela
    flag_vec(k,1) = flag;
    iter_vec(k,1) = iter(2);
    error_vec(k,1) = erro;
    cond_vec(k,1) = cond(A);
    
    %Vetor utilizado para plotar o gráfico das iterações (começando do 0)
    aux = 0:iter(2);
    resvec = resvec/norm(b);
     length(resvec)
     length(aux)
    
    %Plot do gráfico da norma-2 resíduo a cada iteração
    plot(aux, resvec,'-o', 'LineWidth','2');
    hold on
        
 
end
xlabel('Número de Iterações')
ylabel('Norma Relativa do Resíduo')
legend('n=20','n=50','n=80','n=100', 'n=300', 'Location', 'Northeast')
hold off

%Geração da tabela
%T = table(dim,flag_vec,cond_vec, error_vec, iter_vec)

