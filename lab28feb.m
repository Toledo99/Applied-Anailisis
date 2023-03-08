n = 10; 
u = gallery('orthog', n);
b= 10*rand(n,1); 

for k = 1:6
    lambda= diag((1:10).^k);
    Q=(u*lambda*u');
    f= @(x) x'*Q*x + b'*x; 
    df= @(x) 2*Q*x + b; 
    d2f= @(x) 2*Q;
    x0= zeros(n,1);
    [xmin, iter, Xk, F] = buscalinea28(f,df,d2f,x0,'GD');
    [xmin1, iter1, Xk1, F1] = buscalinea28(f,df,d2f,x0,'Newton');
    [xmin2, iter2, Xk2, F2] = buscalinea28(f,df,d2f,x0,'SR1');
    [xmin3, iter3, Xk3, F3] = buscalinea28(f,df,d2f,x0,'BFGS');
    
    subplot(2,3,k)
    hold on
    %plot(F,'r-') 
    plot(F1,'b--')
    plot(F2,'g:')
    plot(F3,'k:')
    
end
 
 
 