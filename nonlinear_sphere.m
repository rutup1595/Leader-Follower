function []=nonlinear_sphere(tfinal,fi,fd,h,c)
i=1;
k=1;
g=[];
tc=[];
td=[];
Kp=[];
T=[];
Td=[];
x=0.01;
K_process=0.5;
r=3;
%alpha=H^3/4*pi*R^3;

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
    a(i)=pi*(2*r*h(i)-h(i)^2);
    fo(i)=K_process*sqrt(2*9.8*h(i));
    u(k,i)=(fi(i)+fd(i)-fo(i))/a(i);
    
    
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
        a(i)=pi*(2*r*h(i)-h(i)^2);
        fo(i)=K_process*sqrt(2*9.8*h(k,i));
        
        %fo(i)=K_process*alpha*(h(i)^-1.5);
        u(k,i)=(fi(i)+fd(i)-fo(i))/a(i);
        
    end
    
    k=k+1;
    
end
figure(1);plot(t(1,:),h(1,1:i),'r');
% figure(2);plot(fo,h(1,:),'g')

%%______regions____________________________________________________________
j=1;
a=[];
region=[];
while j<5
    if j==1
        a=find(h(1,:)<=2.2);
    elseif j==2
        a=find(h(1,:)<=2.4);
        
    elseif j==3
        a=find(h(1,:)<=2.6);
        
    elseif j==4
        a=find(h(1,:)<=2.8);
        
    end
    region(j)=a(end);
    if j==1
        b=1;
        [model(j) controller(j)]=Reaction_Curve(t(1,1:region(1)),h(1,1:region(1)));
    else
        b=region(j-1)+1;
        [model(j) controller(j)]=Reaction_Curve(t(1,region(j-1)+1:region(j)),h(1,region(j-1)+1:region(j)));
    end
    fprintf('Process gain: %g, Time constant: %g, Time delay: %g\n',model(j).gain, model(j).time_constant, model(j).time_delay)
    g(j)=model(j).gain;
    tc(j)=model(j).time_constant;
    td(j)=model(j).time_delay;
    
    Ga=tf(model(j).gain,[model(j).time_constant 1]);
    set(Ga,'InputDelay',td(j));
 %%____________________controller__________________________________________
    
    switch (c)
        
        case 1,
            Kp(j)=controllerf.PID;
            Ti(j)=1;
            Td(j)=1;
            G=Kp(j);
            
        case 2,
            [Ka,La,Ta]=getfod(Ga);
            [G,Kp(j),Ti(j),Td(j)]=ziegler(3,[Ka,La,Ta,10]);
            
        case 3,
            
            [Ka,La,Ta]=getfod(Ga);
            [G,Kp(j),Ti(j),Td(j)]=chrpid(3,1,[Ka,La,Ta,10,0]);
            
        case 4,
            [Ka,La,Ta]=getfod(Ga);
            [G,Kp(j),Ti(j),Td(j)]=cohenpid(3,[Ka,La,Ta,10]);
            
        case 5,
            [Ka,La,Ta]=getfod(Ga);
            [G,Kp(j),Ti(j),Td(j)]=wjcpid([Ka,La,Ta,10]);
            
        case 6,
            for ic=1:3
                [Ka,La,Ta]=getfod(Ga);
                [G,Kp(j),Ti(j),Td(j)]=optPID(3,1,[Ka,La,Ta,10,ic]);
                
            end
            
    end
%% set point ______________________________________________________________
    e=[];
    aa=[];
    b=[];
    h_sp=4;
    u_feedback=[];
    z=1;
    if j==1
        b(j)=1;
    else
        b(j)=region(j-1)+1;
        
    end
    while b(j)<=(region(j)) && z<5
        e(z)=h_sp-h(1,b(j));
        if b(j)==1 || (j~=1&&b(j)==region(j-1)+1)
            add(z)=e(z);
            
        else
            add(z)=e(z)+e(z-1);
        end
        u_feedback(j,b(j))=Kp(j)*(e(z)+(add(z)/Ti(j)))+1;
        z=z+1;
        b(j)=b(j)+1;
    end
    j=j+1;
end

aa=cat(2,u_feedback(1,1:b(1)),u_feedback(2,b(1)+1:b(2)),u_feedback(3,b(2)+1:b(3)),u_feedback(4,b(3)+1:b(4)));
disp(size(aa))
figure(2);plot(t(1,1:size(aa,2)),aa(1:end),'g');
%____________________________________________________________________-
i=1;
a(i)=pi*(2*r*h(i)-h(i)^2);
fo(i)=K_process*sqrt(2*9.8*h(i));
fi(i)=aa(1,i);
u(i)=(fi(i)-fo(i))/a(i);
h(k,1)=h(1);
while t(i)<tfinal
    
    h_new(i+1)= h(i) + x*u(i);                                  %%h=height ;x=delta ;u=control function
    i=i+1;
    t(i)=t(i-1)+x;
    %noise = fi(1)+ sigma*randn(size(sin(2*t(i))));
    a(i)=pi*(2*r*h(i)-h(i)^2);
    fo(i)=K_process*sqrt(2*9.8*h(k,i));
    fi(i)=aa(1,i);
    %fo(i)=K_process*alpha*(h(i)^-1.5);
    u(i)=(fi(i)-fo(i))/a(i);

end
figure(3);plot(t,h_new,'r')