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
q=0.002;
for i=1:N
    real1(i)=To_in_1;
    real2(i)=To_in_2;
    Tin_2(i)=normrnd(To_in_2,To_in_2*q);
    Tin_1(i)=normrnd(To_in_1,To_in_1*q);
end

P0=Tin_1(1)*q;
x(1)=Tin_1(1);
y1(1)=Tin_1(1);
x(2)=Tin_1(1);
% Q=0.0001;
Q=0;
Pk=P0+Q;
R=Tin_1(2)*q;
Kk=Pk*inv(Pk+R);
x(2)=x(2)+Kk*(Tin_1(2)-x(2));
y1(2)=x(2);
Pk=Pk-Kk*Pk;

for i=2:N-1
x(i+1)=x(i);
Pk=Pk+Q;
R=Tin_1(i+1)*q;
Kk=Pk*inv(Pk+R);
x(i+1)=x(i+1)+Kk*(Tin_1(i+1)-x(i+1));
y1(i+1)=x(i+1);
Pk=Pk-Kk*Pk;
end



P0=Tin_2(1)*q;
x(1)=Tin_2(1);
y2(1)=Tin_2(1);
x(2)=x(1);
% Q=0.0001;
Q=0;
Pk=P0+Q;
R=Tin_2(2)*q;
Kk=Pk*inv(Pk+R);
x(2)=x(2)+Kk*(Tin_2(2)-x(2));
y2(2)=x(2);
Pk=Pk-Kk*Pk;

for i=2:N-1
x(i+1)=x(i);
Pk=Pk+Q;
R=Tin_2(i+1)*q;
Kk=Pk*inv(Pk+R);
x(i+1)=x(i+1)+Kk*(Tin_2(i+1)-x(i+1));
y2(i+1)=x(i+1);
Pk=Pk-Kk*Pk;
end
%入口温度校正






%==========================================================================
q=0.002;
% T11=normrnd(T1,T1*q);
% T22=normrnd(T2,T2*q);
T11= load('T11.txt')';
T22= load('T22.txt')';

%==========================================================================
%??????
L=20;
T11(L)=T1(30)*(1-20*q);
T22(L)=T2(L)*(1+20*q);

%==========================================================================
%w_betat~N(0,Q_beta)
%suppose that N=50 and the covariance=2
Q_beta=2;
w_beta=normrnd(0,Q_beta,1,1);%it should be changed every sample time
beta_0 = 0.5;



x=[T11
   T22];
y=zeros(2,N+1);


% u=[y1(1) 
%      y2(1)];

u=[290 
    350];
y(:,1)=x(:,1);
y(:,2)=A*x(:,1)+B*u;

P_k1=[(x(1,1)*q)^2   0
              0      (x(2,1)*q)^2];
Q=[0.02 0
    0  0.02];

P_k=A*P_k1*A'+Q;
R=[(x(1,2)*q)^2   0
    0      (x(2,2)*q)^2];
K_k=P_k*H'*inv(H*P_k*H'+R);
y(:,2)=y(:,2)+K_k*(x(:,2)-y(:,2));
P_k=(I-K_k*H)*P_k;

for i=2:N
  u=[y1(i) 
     y2(i)];
    y(:,i+1)=A*y(:,i)+B*u;
   P_k=A*P_k*A'+Q;

   r(1)=abs(  (y(1,i+1)-x(1,i+1))./(y(1,i+1)*q ) );
   r(2)=abs(  (y(2,i+1)-x(2,i+1))./(y(2,i+1)*q ) );
   c=4;

    if r(1)<c
       R(1,1)=(x(1,i+1)*q)^2;
   else
      R(1,1)=(x(1,i+1)*q)^2*r(1)/c;
   end
   
      if r(2)<c
       R(2,2)=(x(2,i+1)*q)^2;
   else
      R(2,2)=(x(2,i+1)*q)^2*r(2)/c;
   end
           R=[(x(1,i+1)*q)^2   0
     0      (x(2,i+1)*q)^2];

  R=[R(1,1)  0
      0    R(2,2)];
 
%       R=[(x(1,i+1)*q2)^2*(c^2+r(1)^2)/(2*c^2)   0
%      0      (x(2,i+1)*q1)^2*(c^2+r(2)^2)/(2*c^2)];     
%         R=[(x(1,i+1)*q)^2*(c^2+r(1)^2)/(2*c^2)   0
%        0      (x(2,i+1)*q)^2*(c^2+r(2)^2)/(2*c^2)];   
K_k=P_k*H'*inv(H*P_k*H'+R);
y(:,i+1)=y(:,i+1)+K_k*(x(:,i+1)-y(:,i+1));
P_k=(I-K_k*H)*P_k;
end
t1=load('shui.txt');
plot(1:length(T1),T1,'k:',1:length(x(1,:)),x(1,:),'k*',1:length(y(1,:)),y(1,:),'k-',1:length(t1),t1,'k.')
%  title('冷却水出口温度校正曲线')  
axis([1,51,296,312])  
xlabel('Time/s')
  ylabel('Output temperature of cooling water/K')
  legend('True value','Measurements','Reconciled value 1','Reconciled value 2')
SSRE= sum(  abs(       (y(1,:)-T1)/(T1*q)       ).^2   )

figure
t2=load('ben.txt');
plot(1:length(T2),T2,'k:',1:length(x(2,:)),x(2,:),'k*',1:length(y(2,:)),y(2,:),'k-',1:length(t2),t2,'k.')
%  title('苯出口温度校正曲线')  
axis([1,51,318,334])
  xlabel('Time/s')
  ylabel('Output temperature of benzene/K')
legend('True value','Measurements','Reconciled value 1','Reconciled value 2')
  SSRE= sum(  abs(       (y(2,:)-T2)/(T2*q)       ).^2   )
