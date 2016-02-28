clear all
clear;clc

M1=10;
M2=20;
m1=0.947;
cp1=4.18;
m2=1.25;
cp2=1.9;
% k=1/(1/0.85+1/1.7);
k=0.6;
s=11.511;

dt=1;
N=200;

%====================================================================
A=[1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)
      k*s*dt/(M2*cp2)     1-(m2*cp2+k*s)*dt/(M2*cp2)];
  
B=[m1*cp1*dt/(M1*cp1)   0
      0    m2*cp2*dt/(M2*cp2)];
  
x(:,1)=[305
      330];
u=[290
   350];

H=eye(4);

I=[1 0
   0 1];  
 


%==========================================================================
%get the state matrix
%[Uw_k      1      0                         0                     0                            0  0  0  0
% Ub_k   =  0      1                         0                     0                            0  0  0  0
% Tow_k     m1*cp1 0                   1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)              0  0  0  0
% Tob_k     0      m2*cp2*dt/(M2*cp2)  k*s*dt/(M2*cp2)             1-(m2*cp2+k*s)*dt/(M2*cp2)   0  0  0  0
% beta_wk   0      0                         0                     0                            1  0  0  0
% beta_bk   0      0                         0                     0                            0  1  0  0
% beta_o1k  0      0                         0                     0                            0  0  1  0
% beta_02k  0      0                         0                     0                            0  0  0  1
%]

A_star=[1                    0                         0                     0                            0  0  0  0;
        0                    1                         0                     0                            0  0  0  0;
        m1*cp1*dt/(M1*cp1)   0                   1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)              0  0  0  0;
        0                    m2*cp2*dt/(M2*cp2)  k*s*dt/(M2*cp2)             1-(m2*cp2+k*s)*dt/(M2*cp2)   0  0  0  0;
        0                    0                         0                     0                            1  0  0  0;
        0                    0                         0                     0                            0  1  0  0;
        0                    0                         0                     0                            0  0  1  0;
        0                    0                         0                     0                            0  0  0  1
        ];

    
    
Active_GE_pos=zeros(1,4);  

H_star=[eye(4),diag([0 0 0 0])];%this means no gross error.

%---------------------------------------------------------------|
Q_state=0.002*eye(4);%where [290,350,305,330] is the true value |
%---------------------------------------------------------------|

%-----------------------------------|
% w_state=sqrt(Q_state)*randn(4,N);%  |
w_state=normrnd(0,sqrt(Q_state(1)),4,N);
%-----------------------------------|

RY=diag([0.07,0.07,0.07,0.07]);%5*eye(4);the covariance of the measurement

% V_t=sqrt(RY)*randn(4,N);
V_t=normrnd(0,sqrt(RY(1)),4,N);
%===========get the random walk data=======================
%w_betat~N(0,Q_beta)
%suppose that N=50 and the covariance is supposed 3sigma~10sigma gross
%error

%-----------------
Q_beta=0.07*eye(4);
% Q_beta=5*eye(4);
% w_beta=sqrt(Q_beta)*randn(4,N);
w_beta=normrnd(0,sqrt(Q_beta(1)),4,N);
%-----------------

W=[w_state;w_beta];
%-----------------

X_state=zeros(8,N+1);
Y_measure=zeros(4,N+1);

%the initial value
Tin_1=290;
Tin_2=350;
To_1=305;
To_2=330;
beta_in1=0;
beta_in2=0;
beta_ot1=0;
beta_ot2=0;

X_state(:,1)=[Tin_1; Tin_2; To_1; To_2; beta_in1; beta_in2; beta_ot1; beta_ot2];
%----------------------------------------
Y_measure(:,1)=[Tin_1; Tin_2;To_1; To_2];
%----------------------------------------
X_true=zeros(8,N);
X_true(:,1)=X_state(:,1);
%get the measure of every sample time.
for i=1:N-1
    X_true(:,i+1)=A_star*X_state(:,i);
    X_state(:,i+1)=A_star*X_state(:,i)+W(i);
    Y_measure(:,i)=H_star*X_state(:,i)+V_t(i);
end
Y_measure(:,N)=H_star*X_state(:,N)+V_t(N);
%|-------------------------------------------------|
%|Y_measure also should be -betas(the static gross)|
%|Y_mea_mins_staG                                  |
%|-------------------------------------------------|
temp=Y_measure;
for i=10:N
    Y_measure(1,i)=temp(1,i)+10;
end

% 
% for i=50:N
%     Y_measure(3,i)=temp(3,i)+5;
% end
% for i=20:N
%     Y_measure(1,i)=temp(1,i)+X_state(5,i);
% end






%==========================================================================
%the process of Kalman Filter.

M_a=diag(Active_GE_pos);
%============================
Kerr=100;                     %=actually I dont know how to choose the Kerr
%============================
RY_a=(Kerr-1)*M_a*RY+RY;

X_star=zeros(8,N);



I_p=eye(8);

Q_a=diag([diag(Q_state)',diag(Q_beta)']);

X_star_0=X_state(:,1);
% P_t_min_1=cov(X_state');
P_t_min_1=0.05*eye(8);
%==========================================================================
delta_Rt=zeros(4,N);

Sigma_rt=zeros(4,4);

Z=zeros(4,N);

%---------------------|
beta_s=[0 0 0 0]' ;           %| it also should be a matrix(4,1)
%---------------------|
%---------------------------------------------
%filtering residual error


Z_tip_beta=zeros(4,N);

%obtaining univariate statics
sqr_sigma_0 =zeros(4,N);%every row is the diag value
sqr_sigma_1 =zeros(4,N);%every row is the diag value

%cal HMM
%--------------------------|
pi 	= [0.9 0.1;0.9 0.1;0.9 0.1;0.9 0.1]';%4×2
%--------------------------|
bool_value=zeros(4,N);



%========================
%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------|
%     sima_beta_0
%MSE(sima_beta_0) should be manipulated |
MSE=1*ones(4,1);           %|
%--------------------------|
%%%%%%%%%%%%%%%%%%%%%%%%%
%========================
%cal sigma_infinate
sigma_infinate=cal_sigma_infinate(H_star,A_star,RY_a,Q_a,P_t_min_1,N);

beta_c=3*ones(4,1);

[U_z_c S_z_c V_z_c]=svd(sigma_infinate);

Z_c=U_z_c*sqrt(S_z_c)*V_z_c\beta_c;






%==========================================================================



%==========================================================================
%cal Q_z
%--------------------------|
%Q_z should be manipulated |
Q_z0=10*eye(4);             %|
%--------------------------|

% P_z_t=zeros(4,4);
P_z_t=0*eye(4,4);
% P_z_t=0.5*eye(4); %%%%%   the initial of EWMAF


[Q_z,pz_infinate]=tuning_Q_Z(Q_z0,P_z_t,N)

%==========================================================================
count=0;
%the main part of kalman filter!





for t=1:N
    %start from the 2nd
    
    

    S_t=H_star*(A_star*P_t_min_1*A_star'+Q_a)*H_star'+RY_a;
   
    K_a=(A_star*P_t_min_1*A_star'+Q_a)*H_star'/S_t;
    
    X_star(:,t)=A_star*X_star_0+K_a*(Y_measure(:,t)-H_star*A_star*X_star_0);
    %----------------------------------------------------------------
    %cal the residual errror
    delta_Rt(:,t)=Y_measure(:,t)-beta_s-H_star*X_star(:,t);
    
    %get the covriance matrix of delta_Rt
    Sigma_rt=resiual_error_cov(H_star,A_star,RY_a,K_a,Q_a,P_t_min_1,Kerr);
    
    
    %standardizing the residual error
    [U S V]=svd(Sigma_rt);
    Z(:,t)=(U*sqrt(s)*V)\delta_Rt(:,t);
    
    %filtering the residual error
   
    if (t==1)
        temp_Z_tip_beta=Z(:,t);
    end
   
    
    [Z_tip_beta(:,t),temp_P_z_t]=EWMAF(Q_z,P_z_t,Z(:,t),temp_Z_tip_beta);
    
    [Q_z_t,Q_beta_t]=cal_Qz_Qbeta(Q_z,Z_tip_beta(:,t),sigma_infinate,Z_c);
    
    
    
    
    
    %get univariable statistics
    
    sqr_sigma_0(:,t)=diag(temp_P_z_t);
    
    %get sigma_one
    
    sqr_sigma_1(:,t)=Sigma_one(MSE,Sigma_rt,sqr_sigma_0(:,t));
    
    %cal the  HMM

    [bool_value(:,t),P_ht]=cal_HMM(sqr_sigma_0(:,t),sqr_sigma_1(:,t),pi,Z_tip_beta(:,t));
    
    pos=find(bool_value(:,t)==1);
    
    
    if (~isempty(pos))
        for j=1:length(pos)
            Active_GE_pos(pos(j))=1; 
        end  
    else
        Active_GE_pos=zeros(1,4);
    end
    %----------------------------------------------------------------
    %----------------------|
    %tuning Q_z and Q_beta-|___________________
    Q_z=Q_z_t;
    Q_beta=Q_beta_t;
    Q_a=diag([diag(Q_state)',diag(Q_beta)']);
    %----------------------|
    M_a=diag(Active_GE_pos);
    RY_a=(Kerr-1)*M_a*RY+RY;
    
    %----------------------|
    H_star=[H,diag(Active_GE_pos)];
    if(H_star==[eye(4),eye(4)])
        count=count+1;
    end
    
    X_star_0=X_star(:,t);
    
    pi=P_ht; %update the probability
    
    
    temp_Z_tip_beta=Z_tip_beta(:,t);
    
    
    P_z_t=temp_P_z_t;   %update the p_z_t in filtering the residual errors
    
    P_t=(I_p-K_a*H_star)*(A_star*P_t_min_1*A_star'+Q_a);
    
    P_t_min_1=P_t; 
    beta_s=Y_measure(:,t)-X_star(1:4,t)
    
end



yy=zeros(1,N);
yy0=zeros(1,N);
for t=1:N
    yy(1,t)=X_star(1,t)-X_true(1,t);
    yy0(1,t)=Y_measure(1,t)-temp(1,t);
end




%show the figure 
figure(1)
subplot(4,1,1)
plot(1:N,bool_value(1,:),':r')
% axis([1,N+1,-0.2,1.5])

subplot(4,1,2)
plot(1:N,bool_value(2,:),'-.')
% axis([1,N+1,-0.2,1.5])

subplot(4,1,3)
plot(1:N,bool_value(3,:),'--r')
% axis([1,N+1,-0.2,1.5])

subplot(4,1,4)
plot(1:N,bool_value(4,:))
% axis([1,N+1,-0.2,1.5])

figure(2)

subplot(4,1,1)
plot(1:N,delta_Rt(1,:),'.r')
% axis([1,N+1,-5,5])

subplot(4,1,2)
plot(1:N,delta_Rt(2,:),'.')
% axis([1,N+1,-5,5])
% 
subplot(4,1,3)
plot(1:N,delta_Rt(3,:),'.')
% axis([1,N+1,-5,5])

subplot(4,1,4)
plot(1:N,delta_Rt(4,:),'.')
% axis([1,N+1,-5,5])

figure(3)
plot(1:length(X_star(1,:)),X_star(1,:),'g*',1:length(Y_measure(1,:)),Y_measure(1,:),'b+',1:length(X_true(1,:)),X_true(1,:),'r-.')





figure(4)
plot(1:N,yy(1,:),'b-')
hold on 
plot(1:N,yy0(1,:),'r*')
hold off