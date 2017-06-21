function [Gc,Kp,Ti,Td,H]=chrpid(key,tt,vars)
 K=vars(1); L=vars(2); T=vars(3); N=vars(4); a=K*L/T; Ti=[]; Td=[];
 ovshoot=vars(5); 
 if (tt==1)
     TT=T; 
 else
     TT=L; tt=2;
 end
 if ovshoot==0,
 KK=[0.3,0.35,1.2,0.6,1,0.5; 0.3,0.6,4,0.95,2.4,0.42];
 else
 KK=[0.7,0.6,1,0.95,1.4,0.47; 0.7,0.7,2.3,1.2,2,0.42];
 end
 switch key
 case 1, Kp=KK(tt,1)/a;
 case 2, Kp=KK(tt,2)/a; Ti=KK(tt,3)*TT;
 case {3,4}, Kp=KK(tt,4)/a; Ti=KK(tt,5)*TT; Td=KK(tt,6)*L;
 end
[Gc,H]=writepid(Kp,Ti,Td,N,key);