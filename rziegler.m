function [Gc,Kp,Ti,Td,bta,H]=rziegler(vars)
 K=vars(1); L=vars(2); T=vars(3); N=vars(4); a=K*L/T; Kp=1.2/a;
 Ti=2*L; Td=L/2; Kc=vars(5); Tc=vars(6); kappa=Kc*K; tau=L/T; H=[];
 if (kappa > 2.25 && kappa<15) || (tau>0.16 && tau<0.57)
 bta=(15-kappa)/(15+kappa);
 elseif (kappa<2.25 && kappa>1.5) || (tau<0.96 && tau>0.57)
 mu=4*kappa/9; bta=8*(mu-1)/17; Ti=0.5*mu*Tc;
 elseif (kappa>1.2 && kappa<1.5)
 Kp=5*(12+kappa)/(6*(15+14*kappa)); Ti=0.2*(4*kappa/15+1); bta=1;
 end
 Gc=tf(Kp*[bta*Ti,1],[Ti,0]); 
 nH=[Ti*Td*bta*(N+2-bta)/N,Ti+Td/N,1];
 dH=conv([Ti*bta,1],[Td/N,1]); H=tf(nH,dH); 