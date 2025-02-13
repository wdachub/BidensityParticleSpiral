function phic = critical_volume(r,rhop,alpha,R,Kc,phim)

temp=2.*r.*R./(9.*alpha.*Kc)+1./(rhop-1);
phic=min(phim, 0.5*(sqrt(temp.^2+(8.*r.*R)./(9.*alpha.*Kc))-temp));
end
