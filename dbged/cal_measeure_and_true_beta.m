function [X_true,Y_measure,beta]=cal_measeure_and_true_beta(N,Q_state,RY,MSE)
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
b_in1=0;
b_in2=0;
bo_1=0;
b0_2=0;
X_true=zeros(4,N);
Y_measure=zeros(4,N);
beta=zeros(4,N);
w_state=zeros(4,N);
v=zeros(4,N);
w_beta=zeros(4,N);

for j=1:4
    w_state(j,:)=normrnd(0,sqrt(Q_state(j,j)),1,N);
    v(j,:)=normrnd(0,sqrt(RY(j,j)),1,N);
    w_beta=normrnd(0,sqrt(MSE(j)),1,N);
end
A=[1                    0                         0                     0                          ;
  0                    1                         0                     0                          ;
  m1*cp1*dt/(M1*cp1)   0                   1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)            ;
  0                    m2*cp2*dt/(M2*cp2)  k*s*dt/(M2*cp2)             1-(m2*cp2+k*s)*dt/(M2*cp2) ;];

X_0=[Tin_1 Tin_2 To_1 To_2]';

for i=1:N
    X_true(:,i)     = A*X_0;
    Y_measure(:,i)  = A*X_0+w_state(:,i)+v(:,i);
    X_0=X_true(:,i);
end

beta_0=[b_in1 b_in2 bo_1 b0_2]';
for i=1:N
    beta(:,i)=beta_0+w_beta(:,i);
    beta_0=beta(:,i);
end




end
% A=[1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)
%       k*s*dt/(M2*cp2)     1-(m2*cp2+k*s)*dt/(M2*cp2)];
%   
% B=[m1*cp1*dt/(M1*cp1)   0
%       0    m2*cp2*dt/(M2*cp2)];
%   
