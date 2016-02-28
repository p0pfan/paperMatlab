function sigma_infinate=cal_sigma_infinate(H_star,A_star,RY_a,Q_a,P_t_min_1,N)
    Sigma_rt=zeros(4,4*N);
for i=1:N
    S_t=H_star*(A_star*P_t_min_1*A_star'+Q_a)*H_star'+RY_a;
    K_a=(A_star*P_t_min_1*A_star'+Q_a)*H_star'/S_t;
    K_R=eye(4)-H_star*K_a;
    K_Q=H_star*(eye(8)-K_a*H_star);
    Sigma_rt(:,4*i-3:4*i)= K_R*RY_a*K_R'+K_Q*(A_star*P_t_min_1*A_star'+Q_a)*K_Q';
    P_t=(eye(8)-K_a*H_star)*(A_star*P_t_min_1*A_star'+Q_a);
    P_t_min_1=P_t;
end
    sigma_infinate=Sigma_rt(:,4*N-3:4*N);
end