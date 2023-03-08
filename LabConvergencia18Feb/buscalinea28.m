function[xmin, iter, X, F] = buscalinea28(f, df, d2f, x0, method)
%parámetros
EPS = 1e-8; 
MAXITER =1e4; 
n= size(x0, 1);
X = zeros(n, MAXITER); % 
F = zeros(1, MAXITER); %
DF = zeros(n, MAXITER); %guarda la historia de la derivada
%Inicializacion 
X(:,1) = x0; 
k = 1; 
lp = 1; 
DF(:, 1)= ones(size(x0)); 
df1=norm(DF(:,1)); 
B= eye(n); 
while norm(df1) > EPS && k <= MAXITER
    xk = X(:, k);
    DF(:, k) = df(xk); 
    switch method
        case 'GD'
            Bk = B; 
        case 'Newton'
            Bk = d2f(xk);
       case "SR1"
            if k == 1
                Bk = B;
            else
                sk = X(:,k) - X(:,k-1);
                yk = DF(:,k) - DF(:,k-1);
                vk = yk - Bk*sk;
                Bk = Bk + (vk*vk')/(vk'*sk);
            end
        case "BFGS"
            if k == 1
                Bk = B;
            else
                sk = X(:,k) - X(:,k-1);
                yk = DF(:,k) - DF(:,k-1);
                vk =  Bk*sk;
                Bk = Bk - (vk*vk')/(vk'*sk) + (yk*yk')/(yk'*sk);
            end
    end
    %DF(:, k) = df(xk); 
    df1 = DF(:, k);
    p= -linsolve(Bk, DF(:,k)); 
    alpha= bisect(f, xk, p); 
    k=k+1;
    X(:,k) = X(:, k-1) + alpha*p;
    F(1, k) = f(xk);
   
end
xmin = X(:, k); 
iter = k; 
X = X(:, 1:k); 
F = F(1, 1:k);

