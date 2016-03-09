function [Y_true,Y_measure,beta]=cal_measeure_and_true_beta(N,Q_state,RY,MSE)
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

A=[1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)
      k*s*dt/(M2*cp2)     1-(m2*cp2+k*s)*dt/(M2*cp2)];
  
B=[m1*cp1*dt/(M1*cp1)   0;
   0    m2*cp2*dt/(M2*cp2)];


Tin_1=290;
Tin_2=350;
To_1=305;
To_2=330;

bo_1=0;
b0_2=0;

Y_true=zeros(2,N);
Y_measure=zeros(2,N);
beta=zeros(2,N);
w_state=zeros(2,N);
v=zeros(2,N);
w_beta=zeros(2,N);

for j=1:2
    w_state(j,:)=normrnd(0,sqrt(Q_state(j,j)),1,N);
    v(j,:)=normrnd(0,sqrt(RY(j,j)),1,N);
    w_beta(j,:)=normrnd(0,sqrt(MSE(j)),1,N);
end


X_0=[To_1 To_2]';
U=[Tin_1 Tin_2]';
Y_true_0=X_0;


for i=1:N
    Y_true(:,i)= A*Y_true_0+B*U;
    X_state=A*X_0+B*U+w_state(:,i);
    Y_measure(:,i)= X_state+v(:,i);
    X_0=X_state;
    Y_true_0=Y_true(:,i);
end
beta_0=[bo_1 b0_2]';
for i=1:N
    beta(:,i)=beta_0+w_beta(:,i);
    beta_0=beta(:,i);
end


end