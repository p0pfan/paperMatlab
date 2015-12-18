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
  
x=zeros(2,100);
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
%the augmented state matrix
A_a=[A,O;
     O,I];

B_a=[B;O];
%==========================================================================



%===========get the random walk data=======================
%w_betat~N(0,Q_beta)
%suppose that N=50 and the covariance is supposed 3sigma~10sigma gross
%error

%-----------------
w_beta=zeros(2,N);
Q_beta=36;
w_beta(1,:)=normrnd(0,Q_beta,1,50);
w_beta(2,:)=normrnd(0,Q_beta,1,50);
%-----------------
Beta=zeros(2,50);
Beta(:,1)=[10;5];%the initial value
betai=Beta(:,i);
for j=2:N
    Beta(:,j)=betai+w_beta(:,j-1);
    betai=Beta(:,j);
end
%=========================================================

%==========================================================================
%get the true value of the 
for i=1:N
    T1(i+1)=T1(i)+ (   m1*cp1*(Tin_1-T1(i))-k*s*(T1(i)-T2(i))   )/(M1*cp1)*dt;
    T2(i+1)=T2(i)+ (   m2*cp2*(Tin_2-T2(i))+k*s*(T1(i)-T2(i))   )/(M2*cp2)*dt;
% x(:,i+1)=A*x(:,i)+B*u;
end

%  m1*cp1*(T1(100)-Tin_1);
%  m2*cp2*(T2(100)-Tin_2);
%   k*s*(T2(100)-T1(100));


To_in_1=290;
To_in_2=350;
N=50;
















