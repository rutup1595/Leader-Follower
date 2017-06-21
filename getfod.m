function [K,L,T]=getfod(G,method)
 K=dcgain(G);
 if nargin==1
 [Kc,Pm,wc,wcp]=margin(G); ikey=0; L=1.6*pi/(3*wc); T=0.5*Kc*K*L;
 if isfinite(Kc)
     x0=[L;T];
 while ikey==0
 u=wc*x0(1); v=wc*x0(2);
 FF=[K*Kc*(cos(u)-v*sin(u))+1+v^2; sin(u)+v*cos(u)];
 J=[-K*Kc*wc*sin(u)-K*Kc*wc*v*cos(u), -K*Kc*wc*sin(u)+2*wc*v;
 wc*cos(u)-wc*v*sin(u), wc*cos(u)];
 x1=x0-inv(J)*FF;
 if norm(x1-x0)<1e-8
     ikey=1; 
 else
     x0=x1;
 end
end
 L=x0(1); T=x0(2);
 end
 elseif nargin==2 && method==1
 [n1,d1]=tfderv(G.num{1},G.den{1}); [n2,d2]=tfderv(n1,d1);
 K1=dcgain(n1,d1); K2=dcgain(n2,d2);
 Tar=-K1/K; T=sqrt(K2/K-Tar^2); L=Tar-T;
 end
 function [e,f]=tfderv(b,a)
 f=conv(a,a); na=length(a); nb=length(b);
 e1=conv((nb-1:-1:1).*b(1:end-1),a);
 e2=conv((na-1:-1:1).* a(1:end-1),b); maxL=max(length(e1),length(e2));
 e=[zeros(1,maxL-length(e1)) e1]-[zeros(1,maxL-length(e2)) e2];