function[xmin, X, F] = regconf(f, df, d2f, x0, method,delta, etha)
%parámetros
EPS = 1e-8; 
deltak=delta - EPS; 
MAXITER =1e4; 
n= size(x0, 2);
X = zeros(n, MAXITER); % 
F = zeros(1, MAXITER); %
DF = zeros(n, MAXITER); %guarda la historia de la derivada
%Inicializacion 
X(:,1) = x0; 
k = 1; 
%lp = 1; 
DF(:, 1)= ones(size(x0)); 
df1=norm(DF(:,1)); 
B= eye(n); 
while norm(df1) > EPS && k <= MAXITER
    xk = X(:, k);
    DF(:, k) = df(xk); 
    sk = X(:,k) - X(:,k-1);
    yk = DF(:,k) - DF(:,k-1);
    vk = yk - Bk*sk;
    Bk = Bk + (vk*vk')/(vk'*sk);
    dgk = norm(DF(:,k));
    dk= DF(:,k);
    tao = 1; 
    switch method
        case 'cauchy'
            if dk'*Bk*dk <= 0
                tao = min(dgk^3/deltak*dk'*Bk*dk , 1); 
            end
            ps = -deltak*dk/dgk; 
            pc = tao*ps; 
            pk= pc;
        case 'dogleg'
            pb = -(Bk^-1)*dk; 
            pu = -(dk'*dk/dk'*Bk*dk)*dk; 
            b= 2*pu'*(pb-pu); 
            a= (pb-pu)'(pb-pu); 
            c= pu'*pu - deltak^2;
            tao= ((-b + (b^2 - 4*a*c)^(1/2))/2*a)-1;
            if 0<=tao<=1
                pd = tao*pu;
            else
                pd = pu + (tao-1)*(pb-pu); 
            end
            pk=pd;           
    end
    %DF(:, k) = df(xk); 
    df1 = DF(:, k);
    p= -linsolve(Bk, DF(:,k)); 
    alpha= bisect(f, xk, p); 
    k=k+1;
    X(:,k) = X(:, k-1) + alpha*p;
    F(1, k) = norm(f(xk));
   
end

