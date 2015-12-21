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
      k*s*dt/(M2*cp2)     1-(m2*cp2+k*s)*dt/(M2*cp2)]
  
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
        ]

    
    
%=================================================================    
Active_GE_pos=zeros(1,4);                                   %=====
%chose which instrrument has the active gross error         
pos=[];
if (~isempty(pos))
    for j=1:length(pos)
       Active_GE_pos(j)=1; 
    end  
end
%=================================================================

H_star=[eye(4,4),diag(Active_GE_pos)];%this means no gross error.


Q_state=0.01*eyes(4);%where [290,350,305,330] is the true value

w_state=sqrt(Q_state)*randn(4,N);

RY=0.25*eyes(4);
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
X_true=zeros(8,N);
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

X_true(:,1)=[Tin_1; Tin_2; To_1; To_2; beta_in; beta_in; beta_ot; beta_ot];
Y_measure(:,1)=[Tin_1; Tin_2;To_1; To_2;];

%get the measure of every sample time.
for i=2:N
    X_true(:,i)=A_star*X_true(:,i-1)+W(i-1);
    Y_measure(:,i)=H_star*X_true(:,i)+V(i);
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

%P0 and Q0 should be gievn!


%the modified P and Q
%|-------------------------|
%|should use a cycle!!!!!!!|
%|-------------------------|
%|-------------------------|
%|modify                   |
%|-------------------------|


I_p=eye(8);
P_t=[]?
Q_a=[]?
X_star(0)=?_
Y_mea_mins_staG
P_t_min_1
%==========================================================================
%the main part of kalman filter!
S_t=H_star*[A_star*P_t*A_star'+Q_t]*H_star'+RY_a;
K_a=[A_star*P_t*A_star'+Q_t]*H_star'*inv(S_t);
X_star(t)=A_star*X_star(t-1)+K_a*(Y_mea_mins_staG-H_star*A_star*X_star(t-1));
P_t=[I_p-K_a*H_star]*[A_star*P_t_min_1*A_star'+Q_a];
%==========================================================================











%============================================
%this is used to cal the residual error
P_m_t=[I,O;O,1/Kerr*I]*P_t*[I,O;O,1/Kerr*I];
if(isempty(pos))
    Q_m_t=[I,O;O,O]*Q_t*[I,O;O,O];
end
Q_m_t=[I,O;I,O]*Q_t*[I,O;I,O];
%============================================

%==========================================================================



















