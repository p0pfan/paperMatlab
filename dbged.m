clear;clc

M1=10;
M2=20;
m1=0.947;
cp1=4.18;
m2=1.25;
cp2=1.9;
k=1/(1/0.85+1/1.7);
k=0.6;
s=11.511;

dt=1;
N=50;

%====================================================================
A=[1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)
      k*s*dt/(M2*cp2)     1-(m2*cp2+k*s)*dt/(M2*cp2)];
  
B=[m1*cp1*dt/(M1*cp1)   0
      0    m2*cp2*dt/(M2*cp2)];
  
x(:,1)=[305
      330];
u=[290
   350];

H=[1 0
   0 1];

I=[1 0
   0 1];  
 
O = zeros(2,2);

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

A_star=[1      0                         0                     0                            0  0  0  0;
        0      1                         0                     0                            0  0  0  0;
        m1*cp1 0                   1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)              0  0  0  0;
        0      m2*cp2*dt/(M2*cp2)  k*s*dt/(M2*cp2)             1-(m2*cp2+k*s)*dt/(M2*cp2)   0  0  0  0;
        0      0                         0                     0                            1  0  0  0;
        0      0                         0                     0                            0  1  0  0;
        0      0                         0                     0                            0  0  1  0;
        0      0                         0                     0                            0  0  0  1
        ];

    
    
%=================================================================    
Active_GE_pos=zeros(1,4);                                   
%chose which instrrument has the active gross error         

%----------|
pos=[];   %|
%----------|

if (~isempty(pos))
    for j=1:length(pos)
       Active_GE_pos(j)=1; 
    end  
end
%=================================================================

H_star=[eye(4),diag(Active_GE_pos)];%this means no gross error.


Q_state=0.1*eye(4);%where [290,350,305,330] is the true value

w_state=sqrt(Q_state)*randn(4,N);

RY=0.5*eye(4);
V=sqrt(RY)*randn(4,N);
%===========get the random walk data=======================
%w_betat~N(0,Q_beta)
%suppose that N=50 and the covariance is supposed 3sigma~10sigma gross
%error

%-----------------
Q_beta=5*Q_state;
w_beta=sqrt(Q_beta)*randn(4,N);
%-----------------
W=[w_state;w_beta];
%-----------------
X_state=zeros(8,N);
Y_measure=zeros(4,N);
%the initial value
Tin_1=290;
Tin_2=350;
To_1=305;
To_2=330;
beta_in1=5;
beta_in2=5;
beta_ot1=5;
beta_ot2=5;

X_state(:,1)=[Tin_1; Tin_2; To_1; To_2; beta_in1; beta_in2; beta_ot1; beta_ot2];
Y_measure(:,1)=[Tin_1; Tin_2;To_1; To_2];

%get the measure of every sample time.
for i=2:N
    X_state(:,i)=A_star*X_state(:,i-1)+W(i-1);
    Y_measure(:,i)=H_star*X_state(:,i)+V(i);
end


%|-------------------------------------------------|
%|Y_measure also should be -betas(the static gross)|
%|Y_mea_mins_staG                                  |
%|-------------------------------------------------|


%===========================================================================
%the process of Kalman Filter.

M_a=diag(Active_GE_pos);
%============================
Kerr=100;                     %=actually I dont know how to choose the Kerr
%============================
RY_a=(Kerr-1)*M_a*RY+RY;

X_star=zeros(8,N);
K_a=zeros(8,4*(N+1));

I_p=eye(8);
Q_a=diag([diag(Q_state)',diag(Q_beta)']);

X_star(:,1)=X_state(:,1);
P_t_min_1=eye(8);
%==========================================================================
%the main part of kalman filter!
for t=2:N
    %start from the 2nd
    
    
    S_t=H_star*(A_star*P_t_min_1*A_star'+Q_a)*H_star'+RY_a;
   
    K_a(:,4*t-3:4*t  )=(A_star*P_t_min_1*A_star'+Q_a)*H_star'*inv(S_t);
    
    X_star(:,t)=A_star*X_star(:,t-1)+K_a(:,4*t-3:4*t )*(Y_measure(:,t-1)-H_star*A_star*X_star(:,t-1));
    P_t=(I_p-K_a(:,4*t-3:4*t  )*H_star)*(A_star*P_t_min_1*A_star'+Q_a)
    
    P_t_min_1=P_t;
    P_t_min_1
    
end
    %=================================================
    K_a(:,4*N+1:4*(N+1))=(A_star*P_t_min_1*A_star'+Q_a)*H_star'*inv(S_t);
%==========================================================================
%show the figure after Kalman Filter.
plot(X_star(1,:),'*r')
hold on 
plot(Y_measure(1,:),'+')
hold on
plot(X_state(1,:),'k:')
hold off
axis([1,51,280,300])

%===========================================================================
%after get the //STANDARD RESIDUAL ERROR//  


%---------------------|
beta_s=0 ;           %| it also should be a matrix(4,1)
%---------------------|

delta_Rt=zeros(4,N);
Sigma_rt=zeros(4,4*N);
Z=zeros(4,N);
for t=1:N
    delta_Rt(:,t)=Y_measure(:,t)-beta_s-H_star*X_star(:,t);
    
%============================================================================ 
K_a(:,4*t+1:4*(t+1));
    %get the covriance matrix of delta_Rt
    K_R=eye(4)-H_star*K_a(:,4*t+1:4*(t+1));
    K_Q=H_star*(eye(8)-K_a(:,4*t+1:4*(t+1))*H_star);
    Sigma_rt(:,4*t-3:4*t)= K_R*RY_a*K_R'+K_Q*(A_star*P_t_min_1*A_star'+Q_a)*K_Q';
    Z(:,t)=inv(Sigma_rt(:,4*t-3:4*t)^0.5)*delta_Rt(:,t);
end




















