function [Gc,H]=writepid(Kp,Ti,Td,N,key)
 switch key
 case 1, Gc=Kp;
 case 2, Gc=tf(Kp*[Ti,1],[Ti,0]); H=1;
 case 3, nn=[Kp*Ti*Td*(N+1)/N,Kp*(Ti+Td/N),Kp];
 dd=Ti*[Td/N,1,0]; Gc=tf(nn,dd); H=1;
 case 4, d0=sqrt(Ti*(Ti-4*Td)); Ti0=Ti; Kp=0.5*(Ti+d0)*Kp/Ti;
 Ti=0.5*(Ti+d0); Td=Ti0-Ti; Gc=tf(Kp*[Ti,1],[Ti,0]);
 nH=[(1+Kp/N)*Ti*Td, Kp*(Ti+Td/N), Kp];
 H=tf(nH,Kp*conv([Ti,1],[Td/N,1]));
 case 5, Gc=tf(Kp*[Td*(N+1)/N,1],[Td/N,1]); H=1;
 end