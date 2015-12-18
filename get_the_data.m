%1. suppose the state follows Normal distribution
%2. 


I=eye(3,3);
X_error=normrnd(0,1,3,1);
X_t=[1000,500,0]'+[300 0 0;150 40 0;0 0 120]*X_error;

%get the value of state
point_number=10;

data=zeros(3,point_number);
data(:,1)=X_t;
for m=2:point_number
   
    X_t1=0.7*I*X_t+0.3*[1000,500,0]'+[300 0 0;150 40 0;0 0 120]*normrnd(0,1,3,1);
    data(:,m)=X_t1;
    X_t=X_t1;
end

%get the value of measurement
Y=zeros(7,point_number);
for x=1:point_number
    Y(1,x)=data(1,x)+data(3,x)+200*randn;
    Y(2,x)=data(2,x)+data(3,x)+90*randn;
    Y(3,x)=data(3)+130*randn;
    Y(4,x)=data(1,x)+80*randn;
    Y(5,x)=data(1,x)/2.1+data(2,x)+60*randn;
    Y(6,x)=(2.1*data(1,x)+data(2,x))/(data(1,x)+data(2,x))+0.2*randn;
    Y(7,x)=data(2,x)+50*randn;
end
data
Y

%get the matrix of A and C
A=0.7*eye(3)
ith_point=1
C_i=[1 0 1;1 0 1;0 0 1;1 0 0;1/2.1 1 0;2.1/(data(1,ith_point)+data(2,ith_point)),1/(data(1,ith_point)+data(2,ith_point)),0;0 1 0]


 %calculate the kalman filter

 
 data0=data(:,1)
 X_ba=[1000,500,0]'
 ppp=(data0-X_ba)*(data0-X_ba)'
 P0=cov(data0)
 %just only test in one instrument
 
 %when t=0
 w_beta=0;
 beta0=0
 beta_s=0
 %If the active gross error exists
 
 
 