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

% MSE=[10 10]';
MSE=[25 25]';
% Q_beta=0.002*eye(4);%control the covriance of random walk(the gross error)

%---------------------|
beta_s=[0 0]' ;  %| it also should be a matrix(4,1)
%---------------------|
%__________________________________________________________________________
[X_true,Y_measure,Beta]=cal_measeure_and_true_beta(N,Q_state,RY,MSE);

load('y_measure.mat')
load('x_true.mat')
% load('beta.mat')
temp=Y_measure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
one=ones(2,N);
Y_measure(1,20:N)=temp(1,20:N)+10*one(1,20:N);
%  Y_measure(1,20:40)=temp(1,20:40)+Beta(1,20:40);
%  Y_measure(2,70:80)=temp(2,70:80)+Beta(2,70:80);
% Y_measure(2,70:80)=temp(2,70:80)+6*one(2,70:80);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Active_GE_pos=zeros(1,2);

%==========================================================================
M_a=diag(Active_GE_pos);
%========
Kerr=1000; %=actually I dont know how to choose the Kerr
%========
RY_a=(Kerr-1)*M_a*RY+RY;
Q_a=diag([diag(Q_state)',MSE']);
X_star_0=[Y_measure(:,1)' Beta(:,1)']';
% P_t_min_1=cov(X_state');
P_t_min_1=diag([0.36,0.43,0.1,0.1])%0.01*eye(4);
%cal sigma_infinate
sigma_infinate=cal_sigma_infinate(H_star,A_star,RY_a,Q_a,P_t_min_1,N);
sigma_infinate=[0.650553818992836,-0.0230089260714441;-0.0230089260714441,0.647115929037386;];
%==========================================================================
%--------------------------|
pi 	= [0.9 0.9;
       0.1 0.1];%2�2
%--------------------------|
beta_c=2*ones(2,1);         
%--------------------------|
% P_z_t=zeros(4,4);
P_z_t=0.5*eye(2,2);
Q_z=diag([0.005 0.007]);


%==========================================================================
%cal Z_c
[V_z_c D_z_c]=eig(sigma_infinate);
temp_inv_D_z_c=sqrt(D_z_c)^-1;
Z_c=V_z_c*temp_inv_D_z_c/V_z_c*beta_c;
% Z_c=[2.524867344594904;
%      2.531450500598172];
%==========================================================================
%the main part of kalman filter!

temp_p_zt_ht=zeros(2,N);
temp_Q=zeros(2*N,2); 
temp_Q_beta=zeros(2*N,2);
temp_rt=zeros(2,2*N);
 temp_kk=zeros(4,2*N);
for t=1:N
    
    %predict
    X_tip_k_k_minus_1=A_star*X_star_0+B_star*U_star;
    P_k_k_minus_1=A_star*P_t_min_1*A_star'+Q_a;
    
    %update
    y_tip_k=Y_measure(:,t)-H_star*X_tip_k_k_minus_1;
    S_k=H_star*P_k_k_minus_1*H_star'+RY_a;
    K_k=P_k_k_minus_1*H_star'/S_k;
    X_star(:,t)=X_tip_k_k_minus_1+K_k*y_tip_k;
   
    delta_Rt(:,t)=Y_measure(:,t)-beta_s-H_star*X_star(:,t);
    for n=1:2
       if delta_Rt(n,t)<-7
           delta_Rt(n,t)=-delta_Rt(n,t);
       end
    end
    Sigma_rt=resiual_error_cov(H_star,A_star,RY_a,K_k,Q_a,P_t_min_1,Kerr);
    temp_rt(:,2*t-1:2*t)=Sigma_rt;
    [V D]=eig(Sigma_rt);
    temp_inv_D=sqrt(D)^-1;
    temp_kk(:,2*t-1:2*t)=K_k;
    Z(:,t)=V*temp_inv_D/V*delta_Rt(:,t);

    if (t==1)
        temp_Z_tip_beta=Z(:,t);
    end
   
    [Q_z_t,Q_beta_t]=cal_Qz_Qbeta(Q_z,temp_Z_tip_beta,sigma_infinate,Z_c);
    
    temp_Q(2*t-1:2*t,:)=Q_z_t;
    temp_Q_beta(2*t-1:2*t,:)=Q_beta_t;
%     [Z_tip_beta(:,t),temp_P_z_t]=EWMAF(Q_z_t,P_z_t,Z(:,t),temp_Z_tip_beta)
    [Z_tip_beta(:,t),temp_P_z_t]=EWMAF(Q_z,P_z_t,Z(:,t),temp_Z_tip_beta)

    sqr_sigma_0(:,t)=diag(temp_P_z_t);

    sqr_sigma_1(:,t)=Sigma_one(MSE,Sigma_rt,sqr_sigma_0(:,t));
    
    [bool_value(:,t),P_ht,p_zt_ht]=cal_HMM(sqr_sigma_0(:,t),sqr_sigma_1(:,t),pi,Z_tip_beta(:,t));
    temp_p_zt_ht(:,t)=p_zt_ht;
    pos=find(bool_value(:,t)==1);
    
    if (~isempty(pos))
        for j=1:length(pos)
            Active_GE_pos(pos(j))=1; 
        end  
    else
        Active_GE_pos=zeros(1,2);
    end
    %______________________________________________________________________
    
    Q_z=Q_z_t;
    Q_beta=Q_beta_t;
    Q_a=diag([diag(Q_state)',diag(Q_beta)']);
    %----------------------|
    M_a=diag(Active_GE_pos);
    RY_a=(Kerr-1)*M_a*RY+RY;
    
    %----------------------|
    H_star=[H,diag(Active_GE_pos)];

    
    X_star_0=X_star(:,t);
    
    pi=P_ht; %update the probability
    
    
    temp_Z_tip_beta=Z_tip_beta(:,t);
    
    
    P_z_t=temp_P_z_t;   %update the p_z_t in filtering the residual errors
    
    P_t=(I_p-K_k*H_star)*P_k_k_minus_1;
    
    P_t_min_1=P_t; 
%     beta_s=(eye(2)-diag(Active_GE_pos))*(Y_measure(:,t)-X_star(1:2,t))
%     beat_s(2,1)=Y_measure(2,t)-X_star(2,t)
end


figure(1)
plot(1:length(X_star(1,:)),X_star(1,:),'g*',1:length(Y_measure(1,:)),Y_measure(1,:),'b+',1:length(X_true(1,:)),X_true(1,:),'r-.')

figure(2)
subplot(2,1,1)
plot(1:N,bool_value(1,:),':r')
% axis([1,N+1,-0.2,1.5])

subplot(2,1,2)
plot(1:N,bool_value(2,:),'-.')
% axis([1,N+1,-0.2,1.5])



figure(3)
plot(1:length(X_star(2,:)),X_star(2,:),'g*',1:length(Y_measure(2,:)),Y_measure(2,:),'b+',1:length(X_true(2,:)),X_true(2,:),'r-.')

temp1=X_star(1,:)-X_true(1,:);
temp2=X_star(2,:)-X_true(2,:);
figure(4)
subplot(2,1,1)
plot(1:length(X_star(1,:)),temp1,'g*')
subplot(2,1,2)
plot(1:length(X_star(2,:)),temp2,'g*')


