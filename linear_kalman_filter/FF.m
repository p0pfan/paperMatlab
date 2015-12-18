  function dydx=FF(x,T)
M1=10;
M2=20;
m1=0.947;
cp1=4.18;
m2=1.25;
cp2=1.9;
k=0.6;
s=11.511;
Tin_1=290;
Tin_2=350;
f1= ( m1*cp1*(Tin_1-T(1))-k*s*(T(1)-T(2)) )/(M1*cp1) ;
f2=(m2*cp2*(Tin_2-T(2))+k*s*(T(1)-T(2)))/(M2*cp2)  ;
dydx=[f1
          f2]; 

