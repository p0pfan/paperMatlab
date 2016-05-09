function [k_change,X_star]=onlyKalman()

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
I_p=eye(2);
X_star=zeros(2,N);
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


A_star=[1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)           ;
        k*s*dt/(M2*cp2)             1-(m2*cp2+k*s)*dt/(M2*cp2);];

B_star=[m1*cp1*dt/(M1*cp1)   0                 ;
                        0    m2*cp2*dt/(M2*cp2)];

U_star=[Tin_1 Tin_2]';

H_star=[eye(2)];
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
% load('beta.mat')
% load('beta_random_part.mat')
temp=Y_measure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_measure(1,20:50)=temp(1,20:50)+Beta(1,20:50);
one=ones(2,N);
% Y_measure(2,70:80)=temp(2,70:80)+6*one(2,70:80);
% Y_measure(1,50:N)=temp(1,50:N)+10*one(1,50:N);
%  Y_measure(1,50:100)=temp(1,50:100)+Beta(1,50:100);
Y_measure(2,70:80)=temp(2,70:80)+Beta(2,70:80);
% Y_measure(1,30:60)=temp(1,30:60)+Beta(1,30:60);
% Y_measure(2,70:80)=temp(2,70:80)+10*one(2,70:80);
% Y_measure(2,30:40)=temp(2,30:40)+Beta(2,30:40);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the first kalman
X_star_0=Y_measure(:,1);
% P_t_min_1=0*eye(4);
P_t_min_1=diag([0.36,0.43]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================

%========
Kerr=1000; %=actually I dont know how to choose the Kerr
%========
RY_a=RY;

Q_a=diag(diag(Q_state)');



temp_pt=zeros(2,2*N);



%****
k_change=zeros(2,2*N);
for t=1:N
    
    %predict
    X_tip_k_k_minus_1=A_star*X_star_0+B_star*U_star;
    P_k_k_minus_1=A_star*P_t_min_1*A_star'+Q_a
    
    %update
    y_tip_k=Y_measure(:,t)-H_star*X_tip_k_k_minus_1;
    S_k=H_star*P_k_k_minus_1*H_star'+RY_a;
    K_k=P_k_k_minus_1*H_star'/S_k;
    k_change(:,2*t-1:2*t)=K_k;
    X_star(:,t)=X_tip_k_k_minus_1+K_k*y_tip_k;
%     S_t=H_star*(A_star*P_t_min_1*A_star'+Q_a)*H_star'+RY_a;
%    
%     K_a=((A_star*P_t_min_1*A_star'+Q_a)*H_star')/S_t;
%     
%     X_star(:,t)=A_star*X_star_0+K_a*(Y_measure(:,t)-H_star*X_star_0);
%     
    delta_Rt(:,t)=Y_measure(:,t)-beta_s-H_star*X_star_0;
  
    temp_pt(:,2*t-1:2*t)=P_k_k_minus_1;

    P_t=(I_p-K_k*H_star)*P_k_k_minus_1;
    X_star_0=X_star(:,t);
    P_t_min_1=P_t; 
%     beta_s=(eye(2)-diag(Active_GE_pos))*(Y_measure(:,t)-X_star(1:2,t))
%      beat_s(2,1)=Y_measure(2,t)-X_star(2,t)

end
end

% figure(1)
% plot(1:length(X_star(1,:)),X_star(1,:),'g*',1:length(Y_measure(1,:)),Y_measure(1,:),'b+',1:length(X_true(1,:)),X_true(1,:),'r-.')
% % 
% % figure(2)
% % subplot(2,1,1)
% % plot(1:N,bool_value(1,:),':r')
% % % axis([1,N+1,-0.2,1.5])
% % 
% % subplot(2,1,2)
% % plot(1:N,bool_value(2,:),'-.')
% % % axis([1,N+1,-0.2,1.5])
% % 
% figure(2)
% plot(1:length(delta_Rt(2,:)),delta_Rt(2,:),'g*')
% 
% 
% figure(3)
% plot(1:length(X_star(2,:)),X_star(2,:),'g*',1:length(Y_measure(2,:)),Y_measure(2,:),'b+',1:length(X_true(2,:)),X_true(2,:),'r-.')
% 



