clear;clc
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

% P0=Tin_1(1)*q;
% x(1)=Tin_1(1);
% y1(1)=Tin_1(1);
% x(2)=Tin_1(1);
% % Q=0.0001;
% Q=0;
% Pk=P0+Q;
% R=Tin_1(2)*q;
% Kk=Pk*inv(Pk+R);
% x(2)=x(2)+Kk*(Tin_1(2)-x(2));
% y1(2)=x(2);
% Pk=Pk-Kk*Pk;
% 
% for i=2:N-1
% x(i+1)=x(i);
% Pk=Pk+Q;
% R=Tin_1(i+1)*q;
% Kk=Pk*inv(Pk+R);
% x(i+1)=x(i+1)+Kk*(Tin_1(i+1)-x(i+1));
% y1(i+1)=x(i+1);
% Pk=Pk-Kk*Pk;
% end

% plot(1:length(real1),real1,'--',1:length(Tin_1),Tin_1,'*',1:length(y1),y1,':')
%   title('��ȴˮ����¶�У������')  
%   xlabel('ʱ��t/s')
%   ylabel('��ȴˮ����¶�{T_{in,1}/K}')
%   legend('��ʵֵ','����ֵ','У��ֵ')

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

  
  plot(1:length(real2),real2,'--',1:length(Tin_2),Tin_2,'*',1:length(y2),y2,':')
  title('������¶�У������')  
  xlabel('ʱ��t/s')
  ylabel('������¶�{T_{in,2}/K}')
  legend('��ʵֵ','����ֵ','У��ֵ')
