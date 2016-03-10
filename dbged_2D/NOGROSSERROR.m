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


Tin_1=290;
Tin_2=350;
To_1=305;
To_2=330;

bo_1=0;
b0_2=0;
N=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the initiation
H=eye(2);
I_p=eye(4);
X_star=zeros(4,N);
delta_Rt=zeros(2,N);
Sigma_rt=zeros(2,2);
Z=zeros(2,N);
bool_value=zeros(2,N);
Z_tip_beta=zeros(2,N);
sqr_sigma_0 =zeros(2,N);%every row is the diag value
sqr_sigma_1 =zeros(2,N);%every row is the diag value
%==========================================================================
A=[1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)
      k*s*dt/(M2*cp2)     1-(m2*cp2+k*s)*dt/(M2*cp2)];

B=[m1*cp1*dt/(M1*cp1)   0
   0    m2*cp2*dt/(M2*cp2)];


A_star=[1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)              0  0;
        k*s*dt/(M2*cp2)             1-(m2*cp2+k*s)*dt/(M2*cp2)   0  0;
        0                           0                            1  0;
        0                           0                            0  1;];

B_star=[m1*cp1*dt/(M1*cp1)   0                        0       0;
                        0    m2*cp2*dt/(M2*cp2)       0       0;
                        0    0                        1       0;
                        0    0                        0       1;];

U_star=[Tin_1 Tin_2 0  0]';

H_star=[eye(2) diag([0 0])];
%==========================================================================

Q_state=0.02*eye(2);%where [290,350,305,330] is the true value
RY=diag([0.7,0.7]);%5*eye(4);the covariance of the measurement

MSE=[35 35]';
% Q_beta=0.002*eye(4);%control the covriance of random walk(the gross error)

%---------------------|
beta_s=[0 0]' ;  %| it also should be a matrix(4,1)
%---------------------|
%__________________________________________________________________________
[X_true,Y_measure,Beta]=cal_measeure_and_true_beta(N,Q_state,RY,MSE);

load('y_measure.mat')
load('x_true.mat')
load('beta.mat')
temp=Y_measure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_measure(1,20:50)=temp(1,20:50)+Beta(1,20:50);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
Active_GE_pos=zeros(1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the first kalman
X_star_0=[Y_measure(:,1)' Beta(:,1)']';
% P_t_min_1=0*eye(4);
P_t_min_1=diag([0.36,0.43,0.01,0.01]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
M_a=diag(Active_GE_pos);
%========
Kerr=1000; %=actually I dont know how to choose the Kerr
%========
RY_a=(Kerr-1)*M_a*RY+RY;

Q_a=diag([diag(Q_state)',MSE']);



%cal sigma_infinate
sigma_infinate=cal_sigma_infinate(H_star,A_star,RY_a,Q_a,P_t_min_1,N);
%==========================================================================
%--------------------------|
pi 	= [0.9 0.1;
       0.9 0.1]';%2×2
%--------------------------|
beta_c=1*ones(2,1);         
%--------------------------|
% P_z_t=zeros(4,4);
P_z_t=0.01*eye(2,2);
Q_z=diag([3 3]);


%==========================================================================
%cal Z_c
[V_z_c D_z_c]=eig(sigma_infinate);
temp_inv_D_z_c=sqrt(D_z_c)^-1;
Z_c=V_z_c*temp_inv_D_z_c/V_z_c*beta_c;
Z_c=0.15*ones(2,1);
%==========================================================================
%the main part of kalman filter!

temp_p_zt_ht=zeros(2,N);
temp_Q=zeros(2*N,2);
temp_Q_beta=zeros(2*N,2);
temp_pt=zeros(4,4*N);
temp_sigmart=zeros(2,2*N)
for t=1:N
    
    %predict
    X_tip_k_k_minus_1=A_star*X_star_0+B_star*U_star;
    P_k_k_minus_1=A_star*P_t_min_1*A_star'+Q_a
    
    %update
    y_tip_k=Y_measure(:,t)-H_star*X_tip_k_k_minus_1;
    S_k=H_star*P_k_k_minus_1*H_star'+RY_a;
    K_k=P_k_k_minus_1*H_star'/S_k;
    X_star(:,t)=X_tip_k_k_minus_1+K_k*y_tip_k;
%     S_t=H_star*(A_star*P_t_min_1*A_star'+Q_a)*H_star'+RY_a;
%    
%     K_a=((A_star*P_t_min_1*A_star'+Q_a)*H_star')/S_t;
%     
%     X_star(:,t)=A_star*X_star_0+K_a*(Y_measure(:,t)-H_star*X_star_0);
%     
    delta_Rt(:,t)=Y_measure(:,t)-beta_s-H_star*X_star_0;
  
    Sigma_rt=resiual_error_cov(H_star,A_star,RY_a,K_k,Q_a,P_t_min_1,Kerr);
    temp_pt(:,4*t-3:4*t)=P_k_k_minus_1;
    temp_sigmart(:,2*t-1:2*t)=Sigma_rt;
    P_t=(I_p-K_k*H_star)*P_k_k_minus_1;
    X_star_0=X_star(:,t);
    P_t_min_1=P_t; 
%     beta_s=(eye(2)-diag(Active_GE_pos))*(Y_measure(:,t)-X_star(1:2,t))
%      beat_s(2,1)=Y_measure(2,t)-X_star(2,t)
end


figure(1)
plot(1:length(X_star(1,:)),X_star(1,:),'g*',1:length(Y_measure(1,:)),Y_measure(1,:),'b+',1:length(X_true(1,:)),X_true(1,:),'r-.')
% 
% figure(2)
% subplot(2,1,1)
% plot(1:N,bool_value(1,:),':r')
% % axis([1,N+1,-0.2,1.5])
% 
% subplot(2,1,2)
% plot(1:N,bool_value(2,:),'-.')
% % axis([1,N+1,-0.2,1.5])
% 
figure(2)
plot(1:length(delta_Rt(2,:)),delta_Rt(2,:),'g*')


figure(3)
plot(1:length(X_star(2,:)),X_star(2,:),'g*',1:length(Y_measure(2,:)),Y_measure(2,:),'b+',1:length(X_true(2,:)),X_true(2,:),'r-.')




