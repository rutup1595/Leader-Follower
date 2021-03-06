clc; clear all; close all;
%
vl=1.3;
tfinal=100;
%f=1;
%
i=1;
vf(i)=0;
t(i)=0;
f=0.55;
wl(i)=-sin(f*t(i));
%wl(i)=6;
xl(i)=5;
yl(i)=2;
thtl(i)=pi/2;
xf(i)=3.5;
yf(i)=2.2;
thtf(i)=pi/6;
%a(i)=sin(f*t(i));

%
%phi(i)=atan((yl(i)-yf(i))/((xl(i)-xf(i))))-thtf(i);
%d(i)=sqrt((yl(i)-yf(i))^2+((xl(i)-xf(i)))^2);
phi=-pi/3;
d=0.7;
%
h=0.01;
%
wfmax=pi;
K=0.7692;
vlmax=1.5;
vfmax=4;

while t(i)<tfinal
    %Leader
      
    xl(i+1) = xl(i) + h*vl*cos(thtl(i));
    yl(i+1) = yl(i) + h*vl*sin(thtl(i));
    thtl(i+1) = thtl(i) + h*wl(i);
    %Follower
    
    %
   
    %
    %phi(i)=atan((yl(i)-yf(i))/((xl(i)-xf(i))))-thtf(i);
    %d(i)=sqrt((yl(i)-yf(i))^2+((xl(i)-xf(i)))^2);
    %vf(i)=(vl*cos(thtl(i)-thtf(i)-phi)/cos(phi));
    %
    %a(i)=sin(f*t(i));
    %a(i)=0.4;
   
    beta(i) = thtl(i)-thtf(i);
    alpha=asin(K*d*cos(phi));
    if beta(i)>0
        if beta(i)>alpha
            vf(i)=vfmax;
            wf(i)=wfmax;
        else 
            vf(i)=(vl*cos(beta(i)-phi))/(cos(phi));
            if vf(i)>vfmax
                vf(i)=vfmax;
            end
            wf(i)=(vl*sin(beta(i)))/(d*cos(phi));
        end
    else
        if beta(i) < alpha
            vf(i)=vfmax;
            wf(i)=-wfmax;
        else 
            vf(i)=(vl*cos(beta(i)-phi))/(cos(phi));
            if vf(i)>vfmax
                vf(i)=vfmax;
            end
            wf(i)=(vl*sin(beta(i)))/(d*cos(phi));
        end
    end
     %
    xf(i+1) = xf(i) + h*vf(i)*cos(thtf(i));
    yf(i+1) = yf(i) + h*vf(i)*sin(thtf(i));
    thtf(i+1)= thtf(i) + h*wf(i);  
    %
    i=i+1;
    t(i)=t(i-1)+h;
    wl(i)=-sin(f*t(i));
    %wl(i)=6;
end
%
plot(xl,yl,'r',xf,yf); axis equal; grid on;
%figure; plot(t,a); grid on;
figure; plot(t,thtf); grid on;
figure; plot(t,thtl); grid on;
%figure; plot(t, a); grid on;
%figure; plot(t, tht); grid on;
%
%tht_d(1)=0;
%for j=2: length(tht)
%   tht_d(j)=(tht(j)-tht(j-1))/h;
%end
%figure(2); hold on;  plot(t, tht_d, 'r');