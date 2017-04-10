clc;clear all;
syms t x 
%a=((((x(2)-x(5))*cos(x(4)-pi/3)-d+(x(3)-x(6))*sin(x(4)-pi/3)))+(1.3*cos(x(1)-x(4)+pi/3)))/(cos(pi/3));
 %b=((-(x(2)-x(5))*sin(x(4))-1.2*sin(-pi/3)+(x(3)-x(6))*cos(x(4)))+(1.3*sin(x(1)-x(4))))/(1.2*cos(-pi/3));
%g = @(t,x)[-sin(1.1*t);1.3*cos(x(1));1.3*sin(x(1));((-(x(2)-x(5))*sin(x(4))-1.2*sin(-pi/3)+(x(3)-x(6))*cos(x(4)))+(1.3*sin(x(1)-x(4))))/(1.2*cos(-pi/3));
%(((((x(2)-x(5))*cos(x(4)-pi/3)-1.2+(x(3)-x(6))*sin(x(4)-pi/3)))+(1.3*cos(x(1)-x(4)+pi/3)))/(cos(pi/3)))*sin(x(4))];
%[t,xa] = ode45(@(t,x) g(t,x),[0 10],[pi/2 5 2 pi/6 0 1]);
%h = @(t,xf)[];

 

%for n = 1:length(t)
    g = @(t,x)[0.6;1.3*cos(x(1));1.3*sin(x(1));((-(x(2)-x(5))*sin(x(4))-1.2*sin(-pi/3)+(x(3)-x(6))*cos(x(4)))+(1.3*sin(x(1)-x(4))))/(1.2*cos(-pi/3));
(((((x(2)-x(5))*cos(x(4)-pi/3)-1.2+(x(3)-x(6))*sin(x(4)-pi/3)))+(1.3*cos(x(1)-x(4)+pi/3)))/(cos(pi/3)))*cos(x(4));
(((((x(2)-x(5))*cos(x(4)-pi/3)-1.2+(x(3)-x(6))*sin(x(4)-pi/3)))+(1.3*cos(x(1)-x(4)+pi/3)))/(cos(pi/3)))*sin(x(4))];
[t,xa] = ode45(@(t,x) g(t,x),[0 20],[pi/2 5 2 pi/6 3.5 2.2]);
x1=xa(:,2);
x2=xa(:,3);x3=xa(:,5);x4=xa(:,6);
figure(1)  
filename = 'leader_odecircle.gif';

     for n=1:length(t)
     % multicomet(xa(:,[2 5]),xa(:,[3 6]))
    % plot(xa(:,2),xa(:,3),'g',xa(:,5),xa(:,6),'k');axis equal;grid on;
    plot(x1(1:n),x2(1:n),'r',x3(1:n),x4(1:n),'g');axis equal;grid on;title('Sinusoidal input');xlabel('x-coordinate');ylabel('y-coordinate');legend('leader','follower','Location','NorthWest');
      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif','DelayTime',0, 'Loopcount',inf);
      
          else
          imwrite(imind,cm,filename,'gif','DelayTime',0,'WriteMode','append');
      end
end


%plot(xa(:,2),xa(:,3),'g',xa(:,5),xa(:,6),'k');axis equal;grid on;
%figure;plot(t,xa(:,1));axis equal;grid on;