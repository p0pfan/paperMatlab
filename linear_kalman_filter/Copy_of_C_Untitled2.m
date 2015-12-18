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
      k*s*dt/(M2*cp2)     1-(m2*cp2+k*s)*dt/(M2*cp2)];
  
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

end







q=0.002;
% T11=normrnd(T1,T1*q);
%  T22=normrnd(T2,T2*q);
T11= load('T11.txt')';
T22= load('T22.txt')';
L=20;
T11(L)=T1(30)*(1-20*q);
T22(L)=T2(L)*(1+20*q);

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
  u=[290 
    350];
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
%  R=[R(1,1)  0
%      0    R(2,2)];
%       R=[(x(1,i+1)*q2)^2*(c^2+r(1)^2)/(2*c^2)   0
%      0      (x(2,i+1)*q1)^2*(c^2+r(2)^2)/(2*c^2)];     
        R=[(x(1,i+1)*q)^2*(c^2+r(1)^2)/(2*c^2)   0
       0      (x(2,i+1)*q)^2*(c^2+r(2)^2)/(2*c^2)];   
K_k=P_k*H'*inv(H*P_k*H'+R);
y(:,i+1)=y(:,i+1)+K_k*(x(:,i+1)-y(:,i+1));
P_k=(I-K_k*H)*P_k;
end
plot(1:length(T1),T1,'--',1:length(x(1,:)),x(1,:),'*',1:length(y(1,:)),y(1,:),':')
  title('冷却水出口温度校正曲线')  
  xlabel('时间t/s')
  ylabel('冷却水出口温度{T_1}/K')
  legend('真实值','测量值','校正值')
SSRE= sum(  abs(       (y(1,:)-T1)./(T1*q)       ).^2   )

figure
plot(1:length(T2),T2,'--',1:length(x(2,:)),x(2,:),'*',1:length(y(2,:)),y(2,:),':')
  title('苯出口温度校正曲线')  
  xlabel('时间t/s')
  ylabel('苯出口温度{T_2}/K')
  legend('真实值','测量值','校正值')
  SSRE= sum(  abs(       (y(2,:)-T2)./(T2*q)       ).^2   )
