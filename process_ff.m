function[Ga_d]=process_ff(tfinal,fi,fd,a,h,c)
i=1;
k=0;
t_f(i)=0.5;
x=0.01;
K_process=1.2;
sigma=0.7;
fo(i)=K_process*h(i);
noise = fi(i) + sigma*randn(size(sin(t_f(i))));
fi(i)=fi(i);%+noise;
u(i)=(fi(i)+fd(i)-fo(i))/a;

% %%__________________fi model________________________________________________
% while t_f(i)<tfinal
%     
%          h(i+1)= h(i) + x*u(i);                                  %%h=height ;x=delta ;u=control function
%          i=i+1; 
%          t_f(i)=t_f(i-1)+x;
%          %noise = fi(1)+ 8*sigma*randn(size(sin(2*t_f(i))));
%          fi(i)=fi(1);%+noise;
%          fo(i)=K_process*h(i); 
%          u(i)=(fi(i)-fo(i))*1/a;
% 
% end

while t_f(i)<tfinal
    
         h(i+1)= h(i) + x*u(i);                                  %%h=height ;x=delta ;u=control function
         i=i+1; 
         t_f(i)=t_f(i-1)+x;
         noise = fi(1)+ 8*sigma*randn(size(sin(2*t_f(i))));
         fi(i)=fi(1);%+noise;
         fd(i)=fd(1);
         fo(i)=K_process*h(i); 
         u(i)=(fi(i)+fd(i)-fo(i))*1/a;

end
figure(6);plot(t_f,h,'r');title('added')
%%___________________________________________________________________
%%__________filter___________________________________________
% bta=0.2;                                                               %%filter: y=bta(current)+(1-bta)*previous
% for u=2:i
% h(u)=bta*h(u)+(1-bta)*h(u-1);
% end

%___________________________________________________________________
%%______________ReactionCurve__________________________________________
[modelf,controllerf]=ReactionCurve(t_f,h);
 fprintf('Process gain: %g, Time constant: %g, Time delay: %g\n',modelf.gain, modelf.time_constant, modelf.time_delay)
gf=modelf.gain;
tcf=modelf.time_constant;
tdf=modelf.time_delay;
% We can compare how good the approximation is.
Ga_f = tf(modelf.gain,[modelf.time_constant 1]);
set(Ga_f,'InputDelay',tdf);


%____________fd process reaction___________________________________________
j=1;
t_d(j)=1;
x=0.001;
K_process=1.2;
h_d(j)=h(1);
fo_d(j)=K_process*h_d(j);
u_d(j)=(fd(j)-fo_d(j))/a;
while t_d(j)<tfinal

          
         h_d(j+1)= h_d(j) + x*u_d(j);                                  %%h=height ;x=delta ;u=control function
         j=j+1; 
         t_d(j)=t_d(j-1)+x;
         fd(j)=fd(1);%*heaviside(t_d(j));
         fo_d(j)=K_process*h_d(j); 
         u_d(j)=(fd(j)-fo_d(j))*1/a;

end

%%____________________________________________________________________
%%______________ReactionCurve__________________________________________
[model,controller]=ReactionCurve(t_d,h_d);
 fprintf('Process gain: %g, Time constant: %g, Time delay: %g\n',model.gain, model.time_constant, model.time_delay)
gd=model.gain;
tcd=model.time_constant;
tdd=model.time_delay;
% We can compare how good the approximation is.
Ga_d = tf(model.gain,[model.time_constant 1]);
set(Ga_d,'InputDelay',tdd);
%%_________________________________________________________________________
%%___________________controller____________________________________________
switch (c)
    
    case 1,
        Kp=controllerf.PID;
        Ti=1;
        Td=1;
        G=Kp;
        
    case 2,
        [Ka,La,Ta]=getfod(Ga_f);
        [G,Kp,Ti,Td]=ziegler(3,[Ka,La,Ta,10]);
        
    case 3,
        
        [Ka,La,Ta]=getfod(Ga_f);
        [G,Kp,Ti,Td]=chrpid(3,1,[Ka,La,Ta,10,0]);
        
    case 4,
        [Ka,La,Ta]=getfod(Ga_f);
        [G,Kp,Ti,Td]=cohenpid(3,[Ka,La,Ta,10]);
        
    case 5,
        [Ka,La,Ta]=getfod(Ga_f);
        [G,Kp,Ti,Td]=wjcpid([Ka,La,Ta,10]);
        
    case 6,
        for ic=1:3
            [Ka,La,Ta]=getfod(Ga_f);
            [G,Kp,Ti,Td]=optPID(3,1,[Ka,La,Ta,10,ic]);
            
        end
end
%%_________________________________________________________________________
%%__________________actual reaction curve__________________________________

h_sp=3.5;
Gol=G*Ga_f;
s=tf('s');
%%control variable (output of sys)
Y=feedback(Ga_d,Gol)*fd(1)/s;
[A B C D]=ssdata(Y);
[num den]=ss2tf(A,B,C,D);
Y=tf(num,den);
[y1 t1]=step(Y);

%%________________________________________________________________________
%%___________feedforward___________________________________________________
    
ph=tdd-tdf;
s=tf('s');
Lff=ph;
if ph==0
    
    
    Gff=-gd*((tcf+tdf)*s+1)/(gf*(tcd*s+1))
    [A B C D]=ssdata(Gff);
    [num den]=ss2tf(A,B,C,D);
    Gff=tf(num,den);
elseif ph<0
    %bta=ph/(2*tc(1)*(1-exp(sqrt(-ph/tc(1)))));
    Tff=tcf-(ph+Lff)/1.7;
    Gff=gd*((1/1+s*tcd)-(1*exp(-s*(Lff+ph))/1+s*Tff))*exp(-s*tdd);
    Kff=(gd/gf)-gd*Kp*((ph+Lff*(1-1/1.7)))/Ti;
    [Kc, pm ,wc]=margin(Gff);
    L=1.6*pi/(3*wc); T=0.5*Kc*Kff*L;
    Gff=tf(Kff,[T 1],'InputDelay',L)
else
    Gff=-gd*((tcf*s+1)*exp(-ph)/(gf*(tcd*s+1)))
end
%%_________________________________________________________________________
%%error
s=tf('s');
err=h_sp-h(end);
%%controller ouput
fi_controller=err*G;
%%control variable
Gv=25/0.0833*s+1;
co=fi_controller*0.75-Gff*fd(1);
cv=(Ga_f*co/s)+Ga_d*fd(1)/s;
[y2 t2]=step(cv);
figure(3);step(Y,'g');figure(4);step(cv,'r');hold off;
 figure(5);step(co,'b');hold on; plot(t_f,fd(1),'y');
% 

%%___________things to include_____________________________________________
%              if h(i)==h(i-1)&& k<100                         %%plot tillcomplete t
%                  k=k+1;
%          
%              elseif h(i)==h(i-1)&& k==100
%                  t_f(i)=tfinal;
%              end
