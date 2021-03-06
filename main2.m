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
H=eye(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the initiation
Active_GE_pos=zeros(1,4); 
count=0;
I_p=eye(8);
X_star=zeros(8,N);
delta_Rt=zeros(4,N);
Sigma_rt=zeros(4,4);
Z=zeros(4,N);
bool_value=zeros(4,N);
Z_tip_beta=zeros(4,N);
sqr_sigma_0 =zeros(4,N);%every row is the diag value
sqr_sigma_1 =zeros(4,N);%every row is the diag value
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

H_star=[eye(4),diag([0 0 0 0])];%this means no gross error.
%==========================================================================
%__________________________________________________________________________


Q_state=0.002*eye(4);%where [290,350,305,330] is the true value
RY=diag([0.007,0.007,0.007,0.007]);%5*eye(4);the covariance of the measurement

MSE=[35 35 35 35]';
% Q_beta=0.002*eye(4);%control the covriance of random walk(the gross error)

%---------------------|
beta_s=[0 0 0 0]' ;  %| it also should be a matrix(4,1)
%---------------------|
%__________________________________________________________________________
[X_true,Y_measure,Beta]=cal_measeure_and_true_beta(N,Q_state,RY,MSE);

% load('y_measure.mat')
% load('x_true.mat')

temp=Y_measure;




%==========================================================================
M_a=diag(Active_GE_pos);
%========
Kerr=1000; %=actually I dont know how to choose the Kerr
%========
RY_a=(Kerr-1)*M_a*RY+RY;
Q_a=diag([diag(Q_state)',MSE']);
X_star_0=[Y_measure(:,1)' Beta(:,1)']';
% P_t_min_1=cov(X_state');
P_t_min_1=0.01*eye(8);
%cal sigma_infinate
sigma_infinate=cal_sigma_infinate(H_star,A_star,RY_a,Q_a,P_t_min_1,N);
%==========================================================================


%--------------------------|
pi 	= [0.9 0.1;
       0.9 0.1;
       0.9 0.1;
       0.9 0.1]';%4�2
%--------------------------|
beta_c=1*ones(4,1);         
%--------------------------|
% P_z_t=zeros(4,4);
P_z_t=0.1*eye(4,4);
Q_z=0*eye(4);


%==========================================================================
%cal Z_c
[V_z_c D_z_c]=eig(sigma_infinate);
temp_inv_D_z_c=sqrt(D_z_c)^-1;
Z_c=V_z_c*temp_inv_D_z_c/V_z_c*beta_c;
Z_c=0.3*ones(4,1)
%==========================================================================
%the main part of kalman filter!
for t=1:N
 t
    S_t=H_star*(A_star*P_t_min_1*A_star'+Q_a)*H_star'+RY_a;
   
    K_a=((A_star*P_t_min_1*A_star'+Q_a)*H_star')/S_t;
    
    X_star(:,t)=A_star*X_star_0+K_a*(Y_measure(:,t)-H_star*A_star*X_star_0);

    delta_Rt(:,t)=Y_measure(:,t)-beta_s-H_star*X_star(:,t);
  
    Sigma_rt=resiual_error_cov(H_star,A_star,RY_a,K_a,Q_a,P_t_min_1,Kerr);
  
    [V D]=eig(Sigma_rt);
    temp_inv_D=sqrt(D)^-1;
    Z(:,t)=V*temp_inv_D/V*delta_Rt(:,t);
  
    if (t==1)
        temp_Z_tip_beta=Z(:,t);
    end
   
    [Q_z_t,Q_beta_t]=cal_Qz_Qbeta(Q_z,temp_Z_tip_beta,sigma_infinate,Z_c);
    
    [Z_tip_beta(:,t),temp_P_z_t]=EWMAF(Q_z_t,P_z_t,Z(:,t),temp_Z_tip_beta);
    
    sqr_sigma_0(:,t)=diag(temp_P_z_t);

    sqr_sigma_1(:,t)=Sigma_one(MSE,Sigma_rt,sqr_sigma_0(:,t));
    
    [bool_value(:,t),P_ht]=cal_HMM(sqr_sigma_0(:,t),sqr_sigma_1(:,t),pi,Z_tip_beta(:,t));

    pos=find(bool_value(:,t)==1);
    
    if (~isempty(pos))
        for j=1:length(pos)
            Active_GE_pos(pos(j))=1; 
        end  
    else
        Active_GE_pos=zeros(1,4);
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
    
    P_t=(I_p-K_a*H_star)*(A_star*P_t_min_1*A_star'+Q_a);
    
    P_t_min_1=P_t; 
%     beta_s=Y_measure(:,t)-(eye(4)-diag(Active_GE_pos))*X_star(1:4,t)
    
end



yy=zeros(4,N);
yy0=zeros(4,N);
yy1=zeros(4,N);
for j=1:4
for t=1:N
    yy(j,t)=X_star(j,t)-X_true(j,t);
    yy0(j,t)=Y_measure(j,t)-temp(j,t);
    yy1(j,t)=Y_measure(j,t)-X_true(j,t);
end

end


%show the figure 
figure(1)
plot(1:N,yy(1,:),'b-')
hold on 
plot(1:N,yy0(1,:),'r*')
hold on
plot(1:N,yy1(1,:),'g.')
hold off

figure(2)

subplot(4,1,1)
plot(1:N,yy(1,:),'.r')
% axis([1,N+1,-5,5])

subplot(4,1,2)
plot(1:N,yy(2,:),'.')
% axis([1,N+1,-5,5])
% 
subplot(4,1,3)
plot(1:N,yy(3,:),'.')
% axis([1,N+1,-5,5])

subplot(4,1,4)
plot(1:N,yy(4,:),'.')
% axis([1,N+1,-5,5])

figure(3)
plot(1:length(X_star(1,:)),X_star(1,:),'g*',1:length(Y_measure(1,:)),Y_measure(1,:),'b+',1:length(X_true(1,:)),X_true(1,:),'r-.')




figure(4)
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