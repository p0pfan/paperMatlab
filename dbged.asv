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
Tin_1=290;
Tin_2=350;
T1(1)=305;
T2(1)=330;
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
        
H_star=[eye(4,4),zeros(4,4)];%this means no gross error.


Q_state=0.01*eyes(4);%where [290,350,305,330] is the true value

w_state=sqrt(Q_state)*randn(4,N);


%===========get the random walk data=======================
%w_betat~N(0,Q_beta)
%suppose that N=50 and the covariance is supposed 3sigma~10sigma gross
%error

%-----------------
Q_beta=5*Q_state;
w_beta=sqrt(Q_beta)*randn(4,N);
%-----------------




























%==========================================================================
%get the true value of the 
for i=1:N
    T1(i+1)=T1(i)+ (   m1*cp1*(Tin_1-T1(i))-k*s*(T1(i)-T2(i))   )/(M1*cp1)*dt;
    T2(i+1)=T2(i)+ (   m2*cp2*(Tin_2-T2(i))+k*s*(T1(i)-T2(i))   )/(M2*cp2)*dt;
% x(:,i+1)=A*x(:,i)+B*u;
end


w_x=normrnd(0,305*0.002,N,1)%get the measurement error of the output
%supposed that the

To_in_1=290;
To_in_2=350;
N=50;
















