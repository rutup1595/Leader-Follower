function[]=gain_schedule(tfinal,fi,fd,h,c)
i=1;
k=1;
g=[];
tc=[];
td=[];
Kp=[];
Ti=[];
Td=[];
x=0.01;
K_process=0.5;
r=3;
colorselect='gbry';
%alpha=H^3/4*pi*R^3;

while k<4
    %while k<2
    if k==1||k>1 && fd(i)==0
        i=1;
        fi(i)=fi(1)+(k-1)*0.5;
        fd(i)=0;
        t(k,i)=0.5;
    elseif k~=1 && fd(i)~=0
        i=1;
        fd(i)=fd(1);
        fi(i)=0;
        t(k,i)=0.8;
    end
    h(k,1)=h(1);
    a(i)=pi*(2*r*h(i)-h(i)^2);
    fo(i)=K_process*sqrt(2*9.8*h(i));
    u(k,i)=(fi(i)+fd(i)-fo(i))/a(i);
    
    
    
    while t(k,i)<tfinal
        
        h(k,i+1)= h(k,i) + x*u(i);                                  %%h=height ;x=delta ;u=control function
        i=i+1;
        t(k,i)=t(k,i-1)+x;
        %noise = fi(1)+ sigma*randn(size(sin(2*t(i))));
        if k==1|| k>1 && fd(1)==0
            
            fd(i)=0;
            fi(i)=fi(1)+(k-1)*0.5;
            q=i;
            
        elseif k==4 && fd(1)~=0
            fd(i)=fd(1);%+noise;
            fi(i)=0;
            
        end
        a(i)=pi*(2*r*h(i)-h(i)^2);
        fo(i)=K_process*sqrt(2*9.8*h(k,i));
        
        %fo(i)=K_process*alpha*(h(i)^-1.5);
        u(k,i)=(fi(i)+fd(i)-fo(i))/a(i);
        
        
        
               % figure(2);plot(fo,h(1,:),'g')
    end
     figure(1);plot(t(k,:),h(k,1:i),'color',colorselect(k));hold on;

        %%______regions____________________________________________________________
        [model(k),controller(k)]=Reaction_Curve(t(k,:),h(k,:));
        fprintf('Process gain: %g, Time constant: %g, Time delay: %g\n',model(k).gain, model(k).time_constant, model(k).time_delay)
        g(k)=model(k).gain;
        tc(k)=model(k).time_constant;
        td(k)=model(k).time_delay;
        % We can compare how good the approximation is.
        
        
        if k<4
            q=i;
            
            switch (c)
                
                case 1,
                    Kp(k)=controller(k).PID;
                    Ti(k)=1;
                    Td(k)=1;
                    G(k)=Kp;
                    
                case 2,
                    Ga= tf(model(k).gain,[model(k).time_constant 1]);
                    set(Ga,'InputDelay',td(k));
            
                    [Ka,La,Ta]=getfod(Ga);
                    [G,Kp(k),Ti(k),Td(k)]=ziegler(3,[Ka,La,Ta,10]);
                    
                case 3,
                    Ga= tf(model(k).gain,[model(k).time_constant 1]);
                    set(Ga,'InputDelay',td(k));
            
                    [Ka,La,Ta]=getfod(Ga);
                    [G,Kp,Ti,Td]=chrpid(3,1,[Ka,La,Ta,10,0]);
                    
                case 4,
                    Ga= tf(model(k).gain,[model(k).time_constant 1]);
                    set(Ga,'InputDelay',td(k));
            
                    [Ka,La,Ta]=getfod(Ga);
                    [G,Kp,Ti,Td]=cohenpid(3,[Ka,La,Ta,10]);
                    
                case 5,
                    Ga= tf(model(k).gain,[model(k).time_constant 1]);
                    set(Ga,'InputDelay',td(k));
            
                    [Ka,La,Ta]=getfod(Ga);
                    [G,Kp,Ti,Td]=wjcpid([Ka,La,Ta,10]);
                    
                case 6,
                    for ic=1:3
                        Ga= tf(model(k).gain,[model(k).time_constant 1]);
                        set(Ga,'InputDelay',td(k));
            
                        [Ka,La,Ta]=getfod(Ga);
                        [G,Kp,Ti,Td]=optPID(3,1,[Ka,La,Ta,10,ic]);
                        
                    end
            end
            e=[];
            h_sp=1.25*h(k,end);
            j=1;
            u_feedback=[];
            while j<i
                e(j)=h_sp-h(k,j);
                if j~=1
                    add(j)=e(j)+e(j-1);
                else
                    add(j)=e(j);
                end
                u_feedback(j)=Kp(k)*(e(j)+add(j)/Ti(k))+1;
                j=j+1;
            end
            %fi(1)=fi(1)+0.5;
        else
            Gd= tf(model(k).gain,[model(k).time_constant 1]);
            set(Gd,'InputDelay',td(k));
            
        end
        k=k+1;
    end
% %     %_________________________________________________________________________
% % j=1;
% % h_new=[];
% % u_feedforward=[];
% % u_plant=[];
% % u_output=[];
% % u_new=[];
% % h_new(1)=h(1,1);
% %  fi(j)=fi(1);
% %  fo(j)=K_process*h_new(j);
% % while j<q
% %     fd(1,j)=fd(1);
% %    %_________________feedorward signal____________________________
% %    if fd(1)==0
% %    u_feedforward(1,j)=0;
% %    else
% %    u_feedforward(1,j)=Gff*fd(1,j);
% %    end
% %     %_______________plant model signal_____________________________
% %     u_plant(1,j)=u_feedback(1,j)+u_feedforward(1,j);
% %     %_____________output of plant model______________________________
% %     u_new(j)=(u_plant(1,j)-fo(j))*1/a;
% %     h_new(j+1)= h_new(j) + x*u_new(j);
% %    
% %     %____________total ouput (plant +disturbance)____________________
% %     if fd(1)==0
% %     u_output(1,j)=u_plant(1,j);
% %     else
% %     u_output(1,j)=u_plant(1,j)+u(k,j);
% %     end
% %     %________________________________________________________________
% %     j=j+1;
% %     fi(j)=fi(1);
% %     fo(j)=K_process*h_new(j);
% % end
