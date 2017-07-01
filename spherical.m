function[]=spherical(tfinal,fi,fd,h,c)
i=1;
k=1;
g=[];
tc=[];
td=[];
x=0.01;
f=fi;
K_process=0.5;
r=3;
colorselect='gbryk';
%alpha=H^3/4*pi*R^3;

while k<4
    %while k<2
    if k==1||k>1 && fd(1)==0
        i=1;
        fi(k,i)=f+(k-1)*0.5;
        fd(i)=0;
        t(k,i)=10;
    elseif k~=1 && fd(1)~=0
        i=1;
        fd(i)=fd(1);
        fi(k,i)=0;
        t(k,i)=0.8;
    end
    h(k,1)=h(1,1);
    a(k,i)=pi*(2*r*h(k,i)-h(k,i)^2);
    fo(k,i)=K_process*sqrt(2*9.8*(h(k,i)-h(k,1)));
    u(k,i)=(fi(i)+fd(i)-fo(i))/a(k,i);
    
    
    
    while t(k,i)<tfinal
        
        h(k,i+1)= h(k,i) + x*u(k,i);                                  %%h=height ;x=delta ;u=control function
        i=i+1;
        t(k,i)=t(k,i-1)+x;
        %noise = fi(1)+ sigma*randn(size(sin(2*t(i))));
        if k==1|| k>1 && fd(1)==0
            
            fd(i)=0;
            fi(k,i)=f+(k-1)*0.5;
            q=i;
            
        elseif k==4 && fd(1)~=0
            fd(i)=fd(1);%+noise;
            fi(i)=0;
            
        end
        a(k,i)=pi*(2*r*h(k,i)-h(k,i)^2);
        fo(k,i)=K_process*sqrt(2*9.8*(h(k,i)-h(k,1)));
        %fo(i)=K_process*alpha*(h(i)^-1.5);
        u(k,i)=(fi(k,i)+fd(i)-fo(k,i))/a(k,i);
        
    end
        k=k+1;
end
figure(1);plot(fi(:,1),h(:,end),'color',colorselect(k));
figure(2);
coeffs = polyfit(fi(:,1), h(:,end), 5);
xfit = linspace(fi(1,1), fi(end,1), 80);
yfit = polyval(coeffs, xfit);
hold on;
figure(2); plot(xfit, yfit, 'ro-', 'LineWidth', 2);
grid on;

q=1;
n=10;
while q<6
    if n==10
    h_fit=yfit(1:n+1);
    x_fit=xfit(1:n+1);
    t_fit=linspace(10,20,11);
    else    
    h_fit=yfit(n-10:n);     
    x_fit=xfit(10*q:10*(q+1));
    t_fit=linspace(10*q,10*(q+1),11);
    end
    [model(q,1),controller(q,1)]=Reaction_Curve(t_fit,h_fit);
    g(q,1)=model(q,1).gain;
    tc(q,1)=model(q,1).time_constant;
    td(q,1)=10;
    j=q;
     switch (c)
        
        case 1,
            Kp(k,1)=controllerf.PID;
            Ti(k,1)=1;
            Td(k,1)=1;
            G=Kp(k,1);
            
        case 2,
            Ga=tf(model(j,1).gain,[model(j,1).time_constant 1]);
            set(Ga,'InputDelay',td(j,1));
            [Ka,La,Ta]=getfod(Ga);
            [G,Kp(j,1),Ti(j,1),Td(j,1)]=ziegler(3,[Ka,La,Ta,10]);
            
        case 3,
            Ga=tf(model(j,1).gain,[model(j,1).time_constant 1]);
            set(Ga,'InputDelay',td(j,1));
            [Ka,La,Ta]=getfod(Ga);
            [G,Kp(j,1),Ti(j,1),Td(j,1)]=chrpid(3,1,[Ka,La,Ta,10,0]);
            
        case 4,
            Ga=tf(model(j).gain,[model(j).time_constant 1]);
            set(Ga,'InputDelay',td(j));
            [Ka,La,Ta]=getfod(Ga);
            [G,Kp(j,1),Ti(j,1),Td(j,1)]=cohenpid(3,[Ka,La,Ta,10]);
            
        case 5,
            Ga=tf(model(j,1).gain,[model(j,1).time_constant 1]);
            set(Ga,'InputDelay',td(j,1));
            [Ka,La,Ta]=getfod(Ga);
            [G,Kp(j,1),Ti(j,1),Td(j,1)]=wjcpid([Ka,La,Ta,10]);
            
        case 6,
            Ga=tf(model(j,1).gain,[model(j,1).time_constant 1]);
            set(Ga,'InputDelay',td(j,1));
            for ic=1:3
                [Ka,La,Ta]=getfod(Ga);
                [G,Kp(j,1),Ti(j,1),Td(j,1)]=optPID(3,1,[Ka,La,Ta,10,ic]);
                
            end
            
     end
     uu=1;
     while uu<=length(h_fit)
     e(q,uu)=h_fit(end)-h_fit(uu);
     if n==10
           h_sp(q,1)=h_fit(end-1);
     else
           h_sp(q,1)=h_fit(end);
     end
     if uu==1
         u_feedback(q,uu)=(Kp(q,1)*(e(q,uu)+(e(q,uu)/Ti(q,1))))+1;
     else
         u_feedback1(q,uu)=(Kp(q,1)*(e(q,uu)));
         e(q,uu)=e(q,uu)+e(q,uu-1);
         u_feedback(q,uu)=u_feedback1(q,uu)+(Kp(q,1)*(e(q,uu)/Ti(q,1)))+1;
     end
     % u_feedback(k,b)=Kp(k,1)*(e(k,b));
     uu=uu+1;
     end 
   
    
         i=1;
        h_new(q,i)=h_fit(1);
        t_new(q,i)=t_fit(1);
        a(q,i)=pi*(2*r*h_new(q,i)-h_new(q,i)^2);
        fo(q,i)=K_process*sqrt(2*9.8*(h_new(q,i)-h_fit(1)));
        fi(q,i)=u_feedback(q,1);
        u(q,i)=(fi(q,i)-fo(q,i))/a(q,i);
        x=(t_fit(end)-t_fit(1))/10;
         
        while i<length(t_fit)
            
            h_new(q,i+1)= h_new(q,i) + x*u(q,i);                                  %%h=height ;x=delta ;u=control function
            i=i+1;
            t_new(q,i)=t_new(q,i-1)+x;
            %noise = fi(1)+ sigma*randn(size(sin(2*t(i))));
            a(q,i)=pi*(2*r*h_new(q,i)-h_new(q,i)^2);
            fo(q,i)=K_process*sqrt(2*9.8*(h_new(q,i)-h_fit(1)));
           
            fi(q,i)=u_feedback(q,i);
            u(q,i)=(fi(q,i)-fo(q,i))/a(q,i);
        end
 figure(3);plot(t_new(q,:),h_sp(q,:),'k');hold on;plot(xq,vq2,'o',t_new(q,:),h_new(q,:),'color',colorselect(q));
 n=n+10;    
 q=q+1;        
end

%%disp(h_sp(k-2))
    
    %figure(2);plot(t(q-2,:),h_sp(k-2),'k');hold on;plot(t(k-2,:),u_feedback(k-2,:),'color',colorselect(k-2));
    %figure(3);plot(t(k-1,:),h_sp(k-1),'k');hold on;plot(t(k-1,:),u_feedback(k-1,:),'color',colorselect(k-1));
    
    %figure(5);plot(t(k-2,:),u(k-2,:),'k')
%     figure(3);plot(t(q-2,:),h_sp(q-2,1),'k');hold on;plot(t_new(q-2,:),h_new(q-2,:),'color',colorselect(q-2));
%     figure(4);plot(t(q-1,:),h_sp(q-1,1),'k');hold on;plot(t_new(q-1,:),h_new(q-1,:),'color',colorselect(q-1));
    
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