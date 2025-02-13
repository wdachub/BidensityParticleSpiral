function [z,sol,f1,x1,Niter] =bidensitySolver1(Re,Ri,rhop,alpha,R,r,Kv,Kc,phim, hr,phitotal,absTol)

%alternative way to implement bidensitySolver. 
%the performance of this method needs furthur investigation. 

%% defulat value
if ~exist('absTol','var')
    absTol=1e-3;
end

nIter=0;
%%
temp=2*r*R/(9*alpha*Kc)+1/(rhop(1)-1);
phic1=min(phim, 0.5*(sqrt(temp^2+(8*r*R)/(9*alpha*Kc))-temp));
temp=2*r*R/(9*alpha*Kc)+1/(rhop(2)-1);
phic2=min(phim, 0.5*(sqrt(temp^2+(8*r*R)/(9*alpha*Kc))-temp));


%%
ODEoptions = odeset('Events', @(t,y)EventFunction(t,y,phim),'RelTol', 1e-5, 'AbsTol', 1e-7);
zspan = [0 hr];
sigma0= -alpha*Re*Ri/(r*R)*(hr+(rhop(1)-1)*phitotal(1)+(rhop(2)-1)*phitotal(2));

ratio=log(phitotal(1)/(phitotal(1)+phitotal(2)));
phitotal=sum(phitotal);

% if phitotal<max(phic2,phic1)
%     x0=[phitotal/hr;ratio];
%     x1=[phitotal/hr*0.95;ratio*1.1];
% else
%     x0=[phitotal/hr;ratio];
%     x1=[min(phitotal/hr*1.05,phim);ratio*1.1];
%     
% end

%%
x0=[phitotal/hr;ratio];
x1=[phitotal/hr;ratio*1.1];
    
initialValue = [x0; sigma0;0];
[z,sol] =shotting(zspan,initialValue,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim,hr,ODEoptions);
f0=sol(end,3);%sigma

initialValue = [x1; sigma0;0];
[z,sol] =shotting(zspan,initialValue,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim,hr,ODEoptions);
f1=sol(end,3);%sigma

%%
while(abs(f1)>1e-2)
xnew=(x0(2)*f1-x1(2)*f0)/(f1-f0);
x0(2)=x1(2);
f0=f1;
x1(2)=xnew;
initialValue = [x1;sigma0;0];
[z,sol] =shotting(zspan,initialValue,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim,hr,ODEoptions);
f1=sol(end,3);%sigma
nIter=nIter+1;
end
x2=x1;
%%
x0=[phitotal/hr*0.98;ratio];
x1=[phitotal/hr*0.98;ratio*1.1];
    
initialValue = [x0; sigma0;0];
[z,sol] =shotting(zspan,initialValue,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim,hr,ODEoptions);
f0=sol(end,3);%sigma

initialValue = [x1; sigma0;0];
[z,sol] =shotting(zspan,initialValue,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim,hr,ODEoptions);
f1=sol(end,3);%sigma

%%
while(abs(f1)>1e-2)
xnew=(x0(2)*f1-x1(2)*f0)/(f1-f0);
x0(2)=x1(2);
f0=f1;
x1(2)=xnew;
initialValue = [x1;sigma0;0];
[z,sol] =shotting(zspan,initialValue,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim,hr,ODEoptions);
f1=sol(end,3);%sigma
nIter=nIter+1;
end
x0=x2;

%%

initialValue = [x0; sigma0;0];
[z,sol] =shotting(zspan,initialValue,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim,hr,ODEoptions);
f0=[sol(end,3);sol(end,4)-phitotal];


initialValue = [x1; sigma0;0];
[z,sol] =shotting(zspan,initialValue,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim,hr,ODEoptions);
f1=[sol(end,3);sol(end,4)-phitotal];


Niter=0;
invJ0=eye(2);
while  (max(abs(f1))>1e-3)
    deltax=x1-x0;
    deltaf=f1-f0;
    
    invJ1=invJ0+(deltax-invJ0*deltaf)/(deltax'*invJ0*deltaf)*deltax'*invJ0;
    x2=x1-invJ1*f1;
    
    %
    % if  x2(1)>1
    %     x2(1)=1;
    % elseif x2(1)<0
    %     x2(1)=0;
    % end
    % if  x2(2)>0
    %  x2(2)=0;
    % end
    
%prevent the solution goes outside the domain
    while (x2(1)>1||x2(1)<0||x2(2)>0||norm(deltax,2)>0.2)   
        deltax=deltax/2;
        x1=x0+deltax;
        initialValue = [x1; sigma0;0];
        [z,sol] =shotting(zspan,initialValue,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim,hr,ODEoptions);
        f1=[sol(end,3);sol(end,4)-phitotal];
            deltaf=f1-f0;
        invJ1=invJ0+(deltax-invJ0*deltaf)/(deltax'*invJ0*deltaf)*deltax'*invJ0;
        x2=x1-invJ1*f1;

    end
    
    
    
    x0=x1;
    f0=f1;
    x1=x2;
    invJ0=invJ1;
    
    initialValue = [x1; sigma0;0];
    [z,sol] =shotting(zspan,initialValue,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim,hr,ODEoptions);
    f1=[sol(end,3);sol(end,4)-phitotal];
    
    
    
    Niter=Niter+1;
end


%%
end

function  [z,sol] =shotting(zspan,initialValue,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim,hr,ODEoptions)


[z,sol] = ode15s(@(z,phiSigma)particleEquation(z,phiSigma,Re,Ri,rhop,alpha,R,r,Kv,Kc,phim), zspan, initialValue,ODEoptions);
if z(end)<hr%solver terminate before reaching hr
    sol=[sol;
        sol(end,1),sol(end,2),...
        sol(end,3)+alpha*Re*Ri/(r*R)*(1+(rhop(1)-1)*sol(end,1)+rhop(2)*sol(end,2))*(hr-z(end)),...
        sol(end,4)];
    z=[z;1];
end


end