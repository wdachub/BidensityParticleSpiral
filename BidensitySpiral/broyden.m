fun=@(x)[x(1)^3-2*x(1)+2;x(2)^2-1];

x0=[1;2];
x1=[2;2];
f0=fun(x0);
f1=fun(x1);
Niter=0;
invJ0=eye(2);
while  (max(abs(f1))>1e-10)

deltax=x1-x0;
deltaf=f1-f0;

invJ1=invJ0+(deltax-invJ0*deltaf)/(deltax'*invJ0*deltaf)*deltax'*invJ0;
x2=x1-invJ1*fun(x1);
invJ0=invJ1;
x0=x1;
f0=f1;
x1=x2;
f1=fun(x1);   
Niter=Niter+1;    
end

x1
f1
Niter