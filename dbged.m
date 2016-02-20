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

dt=1;
N=100;

%====================================================================
A=[1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)
      k*s*dt/(M2*cp2)     1-(m2*cp2+k*s)*dt/(M2*cp2)];
  
B=[m1*cp1*dt/(M1*cp1)   0
      0    m2*cp2*dt/(M2*cp2)];
  
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
%get the state matrix
%[Uw_k      1      0                         0                     0                            0  0  0  0
% Ub_k   =  0      1                         0                     0                            0  0  0  0
% Tow_k     m1*cp1 0                   1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)              0  0  0  0
% Tob_k     0      m2*cp2*dt/(M2*cp2)  k*s*dt/(M2*cp2)             1-(m2*cp2+k*s)*dt/(M2*cp2)   0  0  0  0
% beta_wk   0      0                         0                     0                            1  0  0  0
% beta_bk   0      0                         0                     0                            0  1  0  0
% beta_o1k  0      0                         0                     0                            0  0  1  0
% beta_02k  0      0                         0                     0                            0  0  0  1
%]

A_star=[1                    0                         0                     0                            0  0  0  0;
        0                    1                         0                     0                            0  0  0  0;
        m1*cp1*dt/(M1*cp1)   0                   1-(m1*cp1+k*s)*dt/(M1*cp1)  k*s*dt/(M1*cp1)              0  0  0  0;
        0                    m2*cp2*dt/(M2*cp2)  k*s*dt/(M2*cp2)             1-(m2*cp2+k*s)*dt/(M2*cp2)   0  0  0  0;
        0                    0                         0                     0                            1  0  0  0;
        0                    0                         0                     0                            0  1  0  0;
        0                    0                         0                     0                            0  0  1  0;
        0                    0                         0                     0                            0  0  0  1
        ];

    
    
%=================================================================    
Active_GE_pos=zeros(1,4);                                   
%chose which instrrument has the active gross error         

%----------|
pos=[];   %|
%----------|

if (~isempty(pos))
    for j=1:length(pos)
       Active_GE_pos(j)=1; 
    end  
end
%=================================================================

H_star=[eye(4),diag(Active_GE_pos)];%this means no gross error.


Q_state=0.002*eye(4);%where [290,350,305,330] is the true value

w_state=sqrt(Q_state)*randn(4,N);

RY=diag([0.7,0.58,0.61,0.66])%5*eye(4);
V=sqrt(RY)*randn(4,N);
%===========get the random walk data=======================
%w_betat~N(0,Q_beta)
%suppose that N=50 and the covariance is supposed 3sigma~10sigma gross
%error

%-----------------
Q_beta=5*Q_state;
w_beta=sqrt(Q_beta)*randn(4,N);
%-----------------
W=[w_state;w_beta];
%-----------------
X_state=zeros(8,N+1);
Y_measure=zeros(4,N+1);
%the initial value
Tin_1=290;
Tin_2=350;
To_1=305;
To_2=330;
beta_in1=5;
beta_in2=5;
beta_ot1=5;
beta_ot2=5;

X_state(:,1)=[Tin_1; Tin_2; To_1; To_2; beta_in1; beta_in2; beta_ot1; beta_ot2];
Y_measure(:,1)=[Tin_1; Tin_2;To_1; To_2];
X_true=zeros(8,N+1);
%get the measure of every sample time.
for i=2:N+1
    X_true(:,i)=A_star*X_state(:,i-1);
    X_state(:,i)=A_star*X_state(:,i-1)+W(i-1);
    Y_measure(:,i)=H_star*X_state(:,i)+V(i);
end
X_true(:,1)=[];
X_state(:,1)=[];
Y_measure(:,1)=[];
%|-------------------------------------------------|
%|Y_measure also should be -betas(the static gross)|
%|Y_mea_mins_staG                                  |
%|-------------------------------------------------|


%===========================================================================
%the process of Kalman Filter.

M_a=diag(Active_GE_pos);
%============================
Kerr=100;                     %=actually I dont know how to choose the Kerr
%============================
RY_a=(Kerr-1)*M_a*RY+RY;

X_star=zeros(8,N);
K_a=zeros(8,4*(N+1));

I_p=eye(8);
Q_a=diag([diag(Q_state)',diag(Q_beta)']);

X_star(:,1)=X_state(:,1);
P_t_min_1=cov(X_state');
%==========================================================================
%the main part of kalman filter!
for t=2:N
    %start from the 2nd

    S_t=H_star*(A_star*P_t_min_1*A_star'+Q_a)*H_star'+RY_a;
   
    K_a(:,4*t-3:4*t  )=(A_star*P_t_min_1*A_star'+Q_a)*H_star'*inv(S_t);
    
    X_star(:,t)=A_star*X_star(:,t-1)+K_a(:,4*t-3:4*t )*(Y_measure(:,t-1)-H_star*A_star*X_star(:,t-1));
    P_t=(I_p-K_a(:,4*t-3:4*t  )*H_star)*(A_star*P_t_min_1*A_star'+Q_a);
    
    P_t_min_1=P_t; 
end
    %=================================================
    K_a(:,4*N+1:4*(N+1))=(A_star*P_t_min_1*A_star'+Q_a)*H_star'*inv(S_t);
%==========================================================================
%show the figure after Kalman Filter.
plot(1:N,X_star(3,:),'*r')
hold on 
plot(1:N,Y_measure(3,:),'+')
hold on
% plot(1:N,X_state(3,:),'k:')
% hold on
plot(1:N,X_true(3,:),'g:')
hold off
axis([1,101,280,340])

%===========================================================================
%after get the //STANDARD RESIDUAL ERROR//  


%---------------------|
beta_s=0 ;           %| it also should be a matrix(4,1)
%---------------------|

delta_Rt=zeros(4,N);
Sigma_rt=zeros(4,4*N);
Z=zeros(4,N);
for t=1:N
    delta_Rt(:,t)=Y_measure(:,t)-beta_s-H_star*X_star(:,t);
    
%============================================================================ 
    K_a(:,4*t+1:4*(t+1));
    %get the covriance matrix of delta_Rt
    K_R=eye(4)-H_star*K_a(:,4*t+1:4*(t+1));
    K_Q=H_star*(eye(8)-K_a(:,4*t+1:4*(t+1))*H_star);
    Sigma_rt(:,4*t-3:4*t)= K_R*RY_a*K_R'+K_Q*(A_star*P_t_min_1*A_star'+Q_a)*K_Q';
    [U S V]=svd(Sigma_rt(:,4*t-3:4*t));
    Z(:,t)=(U*S^(0.5)*V)\delta_Rt(:,t);
end



%==========================================================================
%==========================================================================
%DATA COLLECTION BY FILTERING RESIDUAL ERROR

D_ii=zeros(4,N);
for i=1:N
    D_ii(:,i)=diag(Sigma_rt(:,4*i-3:4*i));

end

%=================================================
%%%%%%%%%%%%EWMAF
Z_tip_beta=zeros(4,N);

P_z_t=zeros(4,4*N+4);
K_EWMAF=zeros(4,4*N+4);
I=eye(4);   

%=============================
%the initial
Z_tip_beta(:,1)=Z(:,1);
P_z_t(:,5:8)=0.5*eye(4);
P_z=cov(Z');
Q_z=0.05*eye(4);
for i=2:N
    
    K_EWMAF(:,4*i-3:4*i )=(P_z_t(:,4*i-3:4*i )+Q_z)/(P_z_t(:,4*i-3:4*i )+Q_z+I);
    
    Z_tip_beta(:,i)=Z_tip_beta(:,i-1)+ K_EWMAF(:,4*i-3:4*i )*(Z(:,i)-Z_tip_beta(:,i-1));
    
    P_z_t(:,4*i+1:4*i+4)=(I-K_EWMAF(:,4*i-3:4*i ))*(P_z_t(:,4*i-3:4*i )+Q_z);
   
end

%=====================================================
%get univariable statistics
sqr_sigma_0 =zeros(4,N);%every row is the diag value

for i=1:N
    sqr_sigma_0(:,i)=diag(P_z_t(:,4*i+1:4*i+4));%the row represent one variable N:the time
end

%get the sigma_one
%1.cal the MSE of every instrument
tempsum=zeros(4,1);
MSE=zeros(4,1);
prediected=H_star*X_star;
for i=1:4

    for j=1:N
        tempsum(i,1)=((prediected(i,j)-Y_measure(i,j))^0.5);
    end
    MSE(i,1)=tempsum(i,1)/N;
end

%2.cal sigma_beta_zero_s vector
sigma_beta_zero_s_vector=zeros(4,N)
for t=1:N
    [U S V]=svd(Sigma_rt(:,4*t-3:4*t));
    sigma_beta_zero_s_vector(:,t)=(U*S^(0.5)*V)\(diag(MSE))^0.5*eye(4,1);
    %readme:
    %the column of sigma_beta_zero_s_vector is the value of every sample
    %time.
    %in one cloumn, every row means one instrument
end


%3.cal sigma_one
square_sigma_one=zeros(4,N);
for i=1:N
    for j=1:4
        a=(sigma_beta_zero_s_vector(j,i))^2;
        b=sqr_sigma_0(j,i);
        square_sigma_one(j,i)=a+b;
    end
end

%===========================================================|
%this is used to see whether the simga_rt is positive or not|
count=0;                                                   %|
for j=1:N                                                  %|
    [t ccc]=zddc(Sigma_rt(:,4*j-3:4*j));                   %|
    count=count+ccc;                                       %|
end                                                        %|
%===========================================================|


%=====================================================================================
%get the cut-off Z-value

%1. get the cutoff gross error: beta_cutoff

%=================================
%so define the beta_cutoff value!|
infinte=100;                    %|
beta_cut_off=5*ones(4,1);       %| 
%=================================
Z_cutoff=Sigma_rt(:,infinte);


%========================================================================================

%the part of HMM

%1.SUPPOSED THE STATE TRANSITION 
A_CP=[0.95 0.05;0.05 0.95];

%2.cal the value of observe matrix
pi 	= [0.95;0.05];

state_probability=zeros(2,N);
likelihood_value=zeros(8,N);
P_ht=zeros(2,N);
%cal the likelihood_value
for i=1:4
for j=1:N
    likelihood_value(2*i-1,j)=1/sqrt(2*3.14*sqr_sigma_0(1,j))*exp(-0.5*Z_tip_beta(i,j)/sqr_sigma_0(1,j));
    likelihood_value(2*i,j)=1/sqrt(2*3.14*square_sigma_one(1,j))*exp(-0.5*Z_tip_beta(i,j)/square_sigma_one(1,j));
end
end
for i=1:N
    C_CP=[likelihood_value(1,i),0;0,likelihood_value(2,i)];
    P_ht(:,i)=A_CP*pi;
    P_zt_ht=C_CP*P_ht(:,i);
    state_probability(:,i)=P_zt_ht/(P_zt_ht(1,1)+P_zt_ht(2,1));
    pi=P_ht(:,i);
end

bool_value=zeros(1,N);
for t=1:N
    if(state_probability(1,t)>state_probability(2,t))
        bool_value(1,t)=0;
    else
        bool_value(1,t)=1;
    end
end

 plot(bool_value(1,1:100));
 axis([1,102,-0.5,1.5])
