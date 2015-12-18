
clear all;clc
T0=[300
      340];
  [t,T]=ode45(@FF,[0:1:100],T0);
plot(t,T(:,1))
figure
plot(t,T(:,2))
