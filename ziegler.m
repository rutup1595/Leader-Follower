function [Gc,Kp,Ti,Td,H]=ziegler(key,vars)
 Ti=[]; Td=[]; H=1;
 if length(vars)==4,
 K=vars(1); L=vars(2); T=vars(3); N=vars(4); a=K*L/T;
 if key==1, Kp=1/a;
 elseif key==2, Kp=0.9/a; Ti=3.33*L;
 elseif key==3 || key==4, Kp=1.2/a; Ti=2*L; Td=L/2;
 end
 elseif length(vars)==3,
 K=vars(1); Tc=vars(2); N=vars(3);
 if key==1, Kp=0.5*K;
 elseif key==2, Kp=0.4*K; Ti=0.8*Tc;
 elseif key==3 || key==4, Kp=0.6*K; Ti=0.5*Tc; Td=0.12*Tc; 
 end
 elseif length(vars)==5,
 K=vars(1); Tc=vars(2); rb=vars(3); N=vars(5);
 pb=pi*vars(4)/180; Kp=K*rb*cos(pb);
 if key==2, Ti=-Tc/(2*pi*tan(pb));
 elseif key==3||key==4, Ti=Tc*(1+sin(pb))/(pi*cos(pb)); Td=Ti/4; 
 end
 end
[Gc,H]=writepid(Kp,Ti,Td,N,key);