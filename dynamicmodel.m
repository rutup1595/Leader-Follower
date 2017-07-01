function[Ga]=dynamicmodel(tfinal,fi,fd,a,h,c)
i=1;
k=1;

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
    else
        i=1;
        fd(i)=fd(1);
        fi(i)=0;
        t(k,i)=0.8;
    end
    u(k,i)=(fi(i)+fd(i)-fo(i))/a;
    
    %%-------------CALCULATE THE HEIGHT FROM THE PLANT MODEL------------------
    h(k,1)=h(1);
    while t(k,i)<tfinal
        
        h(k,i+1)= h(k,i) + x*u(i);                                  %%h=height ;x=delta ;u=control function
        i=i+1;
        t(k,i)=t(k,i-1)+x;
        if k==1|| k>1 && fd(1)==0
            
            fd(i)=0;
            fi(i)=fi(1);
            q=i;
            
        elseif k>1 && fd(1)~=0
            fd(i)=fd(1);
            fi(i)=0;
            
        end
        fo(i)=K_process*h(k,i);
        u(k,i)=(fi(i)+fd(i)-fo(i))*1/a;
        
    end
    figure(3);plot(t(1,1:q),h(1,1:q),'b');legend('Height vs Time (original)');hold on;
    %%-------CALCULATION OF PLANT TRANSFER FUNCTION----------------------------

    [model(1,k),controller(1,k)]=Reaction_Curve(t(k,:),h(k,:));
    fprintf('Process gain: %g, Time constant: %g, Time delay: %g\n',model(1,k).gain, model(1,k).time_constant, model(1,k).time_delay)
    g(1,k)=model(1,k).gain;
    tc(1,k)=model(1,k).time_constant;
    td(1,k)=model(1,k).time_delay;

    %--------SELECTION OF CONTROLLER------------------------------------------
    if k==1
     
        Ga= tf(model(1,k).gain,[model(1,k).time_constant 1]);
        set(Ga,'InputDelay',td(1,k));
        
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
        
%---------------CALCULATION OF ERROR ----------------------------------------
        
        h_sp=6;
        j=1;
        while j<i
            e(k,j)=h_sp-h(k,j);                       %% ERROR= SET POINT-FEEBACK SIGNAL
            if j==1
                u_feedback(k,j)=(Kp(k,1)*(e(k,j)+(e(k,j)/Ti(k,1))))+1;
                disp('hi')
            else
                u_feedback1(k,j)=(Kp(k,1)*(e(k,j)));
                e(k,j)=e(k,j)+e(k,j-1);
                u_feedback(k,j)=u_feedback1(k,j)+(Kp(k,1)*(e(k,j)/Ti(k,1)))+1;
            end
            j=j+1;
        end
    else
        Gd= tf(model(1,k).gain,[model(1,k).time_constant 1]);       %%DISTURBANCE TRANSFER FUNCTION
        set(Gd,'InputDelay',td(1,k));
        disp('paro disturb')
    end
    k=k+1;
    
end
k=k-1;

%%------CALCULATION OF FEEDFORWARD CONTROLLER------------------------------
ph=td(1,k)-td(1,k-1);
s=tf('s');
Lff=ph;
if ph==0
    Gff=-g(1,k)*((tc(1,k-1)+td(1,k-1))*s+1)/(g(1,k-1)*(tc(1,k)*s+1));
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
    % Gff=-g(k)*((tc(k-1)*s+1)*exp(-ph)/(g(k-1)*(tc(k)*s+1)));
    Gff=-g(1,k)/g(1,k-1);
    disp('paro')
end
%%_________________________________________________________________________
j=1;
% h_new=[];
% u_feedforward=[];
% u_plant=[];
% u_output=[];
% u_new=[];

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
   
    
    %%____________total ouput (plant +disturbance)____________________
%     if fd(1)==0
%         u_output(1,j)=u_plant(1,j);
%     else
%         u_output(1,j)=u_plant(1,j)+u(k,j);
%     end
    %%________________________________________________________________
  
end
j=1;
h_new(1,j)=h(1,1);
fo(1,j)=K_process*h_new(1,j);
u_new(1,j)=(u_plant(1,j)-fo(j))*1/a;
while j<q
    
    h_new(1,j+1)= h_new(1,j) + x*u_new(1,j);
    j=j+1;
    fo(j)=K_process*h_new(j);
    u_new(1,j)=(u_plant(1,j)-fo(j))*1/a;
  
end
figure(6);plot(t(k-1,1:q),h_new,'r',t(k-1,1:q-1),h_sp,'b');legend('new height','set point');title('set point height');
