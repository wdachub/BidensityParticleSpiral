function dh = heightEquation(r,h,Re,Ri,rhop,alpha,R,Kv,Kc,phim)


temp=2*r*R/(9*alpha*Kc)+1/(rhop-1);
phi=min(phim, 0.5*(sqrt(temp^2+(8*r*R)/(9*alpha*Kc))-temp));

if phi>=phim
    dphi=0;
else
    temp1=r*R/(alpha*Kc);
    dphi=0.5*(-4/9*temp1+ (16/9*temp1+8/9*temp1*temp)/(2*sqrt(8/9*r*temp1+temp^2)));
    
end

rho=1+(rhop-1)*phi;
drho=(rhop-1)*dphi;

lambda=alpha/R;
Gamma=1+lambda^2/r^2;

if phi>0.99*phim
    %close to the max volume fraction
    %viscosity is infinity
dh=-h*(3/2*lambda^2/Gamma/r.^3+3/8*drho/rho);
else
 mu=(1-phi/phim)^(-2);
dh=-h*(3/2*lambda^2/Gamma/r.^3+3/8*drho/rho)+...
    6/35*lambda^2*Re^2*Ri*h^4*rho^2/(Gamma^3*r^3*mu^2);
end
    

% dh=-h*(3/2*lambda^2/Gamma/r.^3)+...
%     6/35*lambda^2*Re^2*Ri*h^4/(Gamma^3*r^3);

end

