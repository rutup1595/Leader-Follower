function[Ga]=dynamicmodel(tfinal,fi,fd,a,h,c)
i=1;
k=1;
g=[];
tc=[];
td=[];
x=0.01;
K_process=1.2;
sigma=0.56;
colorstring='bgr';
fo(i)=K_process*h(i);
%noise = fi(i) + sigma*randn(size(sin(t_f(i))));
%%____________________fi model_____________________________________________
while k<3
%while k<2 
if k==1||k~=1 && fd(i)==0
        i=1;
        fi(i)=fi(1);
        fd(i)=0;
        t(k,i)=0.5;   
elseif k~=1 && fd(i)~=0
        i=1;
        fd(i)=fd(1);
        fi(i)=0;
        t(k,i)=0.8;
end
u(k,i)=(fi(i)+fd(i)-fo(i))/a;


h(k,1)=h(1);
    while t(k,i)<tfinal
        
        h(k,i+1)= h(k,i) + x*u(i);                                  %%h=height ;x=delta ;u=control function
        i=i+1;
        t(k,i)=t(k,i-1)+x;
        %noise = fi(1)+ sigma*randn(size(sin(2*t(i))));
        if k==1|| k>1 && fd(1)==0
            
            fd(i)=0;
            fi(i)=fi(1);
            q=i;
            
        elseif k>1 && fd(1)~=0
            fd(i)=fd(1);%+noise;
            fi(i)=0;
           
        end
        fo(i)=K_process*h(k,i);
        u(k,i)=(fi(i)+fd(i)-fo(i))*1/a;
        
    end
    figure(3);plot(t(1,1:q),h(1,1:q),'b');legend('Height vs Time (original)');hold on;
   
    [model(k),controller(k)]=Reaction_Curve(t(k,1:q),h(k,1:q));
    fprintf('Process gain: %g, Time constant: %g, Time delay: %g\n',model(k).gain, model(k).time_constant, model(k).time_delay)
    g(k)=model(k).gain;
    tc(k)=model(k).time_constant;
    td(k)=model(k).time_delay;
    % We can compare how good the approximation is.
    
    
    if k==1
        q=i;
        Ga= tf(model(k).gain,[model(k).time_constant 1]);
        set(Ga,'InputDelay',td(k));
        
        switch (c)
            
            case 1,
                Kp=controllerf.PID;
                Ti=1; 
                Td=1;
                G=Kp;
                
            case 2,
                [Ka,La,Ta]=getfod(Ga(k));
                [G,Kp,Ti,Td]=ziegler(3,[Ka,La,Ta,10]);
                
            case 3,
                
                [Ka,La,Ta]=getfod(Ga(k));
                [G,Kp,Ti,Td]=chrpid(3,1,[Ka,La,Ta,10,0]);
                
            case 4,
                [Ka,La,Ta]=getfod(Ga(k));
                [G,Kp,Ti,Td]=cohenpid(3,[Ka,La,Ta,10]);
                
            case 5,
                [Ka,La,Ta]=getfod(Ga(k));
                [G,Kp,Ti,Td]=wjcpid([Ka,La,Ta,10]);
                
            case 6,
                for ic=1:3
                    [Ka,La,Ta]=getfod(Ga(k));
                    [G,Kp,Ti,Td]=optPID(3,1,[Ka,La,Ta,10,ic]);
                    
                end
        end
        e=[];
        h_sp=6;
        j=1;
        u_feedback=[];
        while j<i
            e(j)=h_sp-h(k,j);
            if j~=1
                add(j)=e(j)+e(j-1);
            else
                add(j)=e(j);
            end
            u_feedback(j)=Kp*(e(j)+add(j)/Ti)+1;
            j=j+1;
        end
    else
        Gd= tf(model(k).gain,[model(k).time_constant 1]);
        set(Gd,'InputDelay',td(k));
        
    end
    k=k+1;
    
end
k=k-1;
%%___________________________________________________________________
%%__________filter___________________________________________
% bta=0.2;                                                               %%filter: y=bta(current)+(1-bta)*previous
% for u=2:i
% h(u)=bta*h(u)+(1-bta)*h(u-1);
% end
%%________________________________________________________________________
%%___________feedforward___________________________________________________

ph=td(k)-td(k-1);
s=tf('s');
Lff=ph;
if ph==0  
    Gff=-g(k)*((tc(k-1)+td(k-1))*s+1)/(g(k-1)*(tc(k)*s+1));
    [A B C D]=ssdata(Gff);
    [num den]=ss2tf(A,B,C,D);
    Gff=tf(num,den);
elseif ph<0
    %bta=ph/(2*tc(1)*(1-exp(sqrt(-ph/tc(1)))));
    Tff=tc(k-1)-(ph+Lff)/1.7;
    Gff=g(k)*((1/1+s*tc(k))-(1*exp(-s*(Lff+ph))/1+s*Tff))*exp(-s*td(k));
    Kff=(g(k)/g(k-1))-g(k)*Kp*((ph+Lff*(1-1/1.7)))/Ti;
    [Kc, pm ,wc]=margin(Gff);
    L=1.6*pi/(3*wc); T=0.5*Kc*Kff*L;
    Gff=tf(Kff,[T 1],'InputDelay',L);
else
    Gff=-g(k)*((tc(k-1)*s+1)*exp(-ph)/(g(k-1)*(tc(k)*s+1)));
    Gff=-g(k)/g(k-1);
end
%%_________________________________________________________________________
j=1;
h_new=[];
u_feedforward=[];
u_plant=[];
u_output=[];
u_new=[];
h_new(1)=h(1,1);
 fi(j)=fi(1);
 fo(j)=K_process*h_new(j);
while j<q
    fd(1,j)=fd(1);
   %%_________________feedorward signal____________________________
   if fd(1)==0
   u_feedforward(1,j)=0;
   else
   u_feedforward(1,j)=Gff*fd(1,j);
   end
    %%_______________plant model signal_____________________________
    u_plant(1,j)=u_feedback(1,j)+u_feedforward(1,j);
    %%_____________output of plant model______________________________
    u_new(j)=(u_plant(1,j)-fo(j))*1/a;
    h_new(j+1)= h_new(j) + x*u_new(j);
   
    %%____________total ouput (plant +disturbance)____________________
    if fd(1)==0
    u_output(1,j)=u_plant(1,j);
    else
    u_output(1,j)=u_plant(1,j)+u(k,j);
    end
    %%________________________________________________________________
    j=j+1;
    fi(j)=fi(1);
    fo(j)=K_process*h_new(j);
end
figure(4);plot(t(k-1,1:q-1),u_feedback,'b');legend('controller input');figure(5);plot(t(k-1,1:q-1),u_plant,'g');legend('plant input');title('control signal');
figure(6);plot(t(k-1,1:q),h_new,'r',t(k-1,1:q-1),h_sp,'b');legend('new height','set point');title('set point height');legend('controlled variable')
figure(7);plot(t(k-1,1:q-1),u_output,'g');title('output');