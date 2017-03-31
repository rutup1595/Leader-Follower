clc; clear all; close all;
%
vl=1.3;
tfinal=40;
%
i=1;
vf(i)=0;
t(i)=0;
%f=0.55;
wl(i)=1;%-sin(1.1*t(i));
%wl(i)=6;
xl(i)=5;
yl(i)=2;
thtl(i)=pi/2;
xf(i)=3;
yf(i)=1.5;
thtf(i)=pi/6;
phi=-pi/3;
d=1.2;
%
h=0.001;
%
wfmax=pi;
K=0.7692;
vlmax=2;
vfmax=4;
%
%

l(i)=0;
etan(i)=0;
alphan(i)=0;
btan(i)=0;
wf(i)=0;    
%
%

while t(i)<tfinal
    %Leader
        %disp('1')
        %Follower
        l(i)=sqrt(((xl(i)-(xf(i)+0.2*cos(thtf(i))))^2)+((yl(i)-(yf(i)+0.2*sin(thtf(i))))^2));
        etan(i)=atan2(yl(i)-(yf(i)+0.2*sin(thtf(i))),xl(i)-(xf(i)+0.2*cos(thtf(i))));
        alphan(i)=etan(i)-thtl(i);
        btan(i)=etan(i)-thtf(i);
        vf(i)=-(d-l(i))*cos(btan(i))+(phi-alphan(i))*l(i)*sin(btan(i))-l(i)*wl(i)+cos(alphan(i)-btan(i));%cos(btan(i))-l(i)*sin(btan(i))*wl(i)-vl*cos(alphan(i)+btan(i));
       % wf(i)=(-(d-l(i))*sin(btan(i))-(phi-alphan(i))*l(i)*cos(btan(i))+l(i)*cos(btan(i))*wl(i)+vl*sin(-alphan(i)+btan(i)))/0.2;
        wf(i)=((d-l(i))*0.5*(-sin(btan(i)))/0.2)+(sin(btan(i)-alphan(i)))/d+sin(btan(i))*wl(i)/l(i)-l(i)*cos(btan(i))*(phi-alphan(i))/0.2;
        xl(i+1) = xl(i) + h*vl*cos(thtl(i));
        yl(i+1) = yl(i) + h*vl*sin(thtl(i));
        thtl(i+1) = thtl(i) + h*wl(i);
        xf(i+1) = xf(i) + h*vf(i)*cos(thtf(i));
        yf(i+1) = yf(i) + h*vf(i)*sin(thtf(i));
        thtf(i+1)= thtf(i) + h*wf(i);    
        i=i+1;
        t(i)=t(i-1)+h;
        wl(i)=1%-sin(1.1*t(i));
       
end
%
plot(xl,yl,'r',xf,yf); axis equal; grid on;