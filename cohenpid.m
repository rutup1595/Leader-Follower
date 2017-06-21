function [Gc,Kp,Ti,Td,H]=cohenpid(key,vars)
 K=vars(1); L=vars(2); T=vars(3); N=vars(4);
 a=K*L/T; tau=L/(L+T); Ti=[]; Td=[];
 switch key
 case 1,Kp=(1+0.35*tau/(1-tau))/a;
 case 2,
 Kp=0.9*(1+0.92*tau/(1-tau))/a; Ti=(3.3-3*tau)*L/(1+1.2*tau);
 case {3,4}, Kp=1.35*(1+0.18*tau/(1-tau))/a;
 Ti=(2.5-2*tau)*L/(1-0.39*tau); Td=0.37*(1-tau)*L/(1-0.81*tau);
 case 5
 Kp=1.24*(1+0.13*tau/(1-tau))/a; Td=(0.27-0.36*tau)*L/(1-0.87*tau);
 end
 [Gc,H]=writepid(Kp,Ti,Td,N,key);