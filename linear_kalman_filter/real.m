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