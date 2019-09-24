mu2=-5:1:5;
sigma2=-5:1:5;
es=sort(-mu2);
stdd=sort(sqrt(1+sigma2.^2));
es=unique(es);
stdd=unique(stdd);
ii=length(es);
jj=length(stdd);
z=zeros(ii,jj);
for i=1:ii
    for j=1:jj
        z(i,j)=I_ij([0 -es(i)],[1 0;0 sqrt(stdd(j)^2-1)],1);
    end
end
contour(z.')