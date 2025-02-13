function psi = streamfunction(z,r,h,Re,Ri,rhop,alpha,R,Kv,Kc,phim)

if abs(rhop-1)<eps
    phi=0;
    dphi=0;
else
    temp=2*r*R/(9*alpha*Kc)+1/(rhop-1);
    phi=min(phim, 0.5*(sqrt(temp.^2+(8*r*R)/(9*alpha*Kc))-temp));
    
    
    
    if phi>=phim
        dphi=0;
    else
        temp1=r*R/(alpha*Kc);
        dphi=0.5*(-4/9*temp1+ (16/9*temp1+8/9*temp1.*temp)./(2*sqrt(8/9*r.*temp1+temp.^2)));
        
    end
end
rho=1+(rhop-1)*phi;
drho=(rhop-1)*dphi;
mu=(1-phi./phim).^(-2);

lambda=alpha/R;
Gamma=1+lambda^2./r.^2;
psi=Re*Ri*z.^2.*(2*z-3*h).*(z-h).*(Gamma.*r.^3.*drho+4*lambda^2*rho)./(48*Gamma.^3.*r.^2.*mu)...
    +alpha^2*Re^3*Ri^2*z.^2.*rho.^3.*(z-2*h).^2.*(z-h).*(z.^2-2*h.*(2*h+z))...
    ./(840*R^2*Gamma.^5.*r.^2.*mu.^3);


%
% psi=h.^4.*Re*Ri*zn.^2.*(2*zn-3).*(zn-1).*(Gamma.*r^3*drho+4*lambda^2*rho)/(48*Gamma.^3*r^2*mu)...
%     +h.^5.*alpha^2*Re^3*Ri^2*zn^2.*rho^3.*(zn-2).^2.*(zn-1).*(zn.^2-2*(2+zn))...
%     /(840*R^2*Gamma.^5.*r^2*mu^3);
end

