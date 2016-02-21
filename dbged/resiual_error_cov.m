
function Sigma_rt=resiual_error_cov(C_a,A_star,RY_a,K_a,Q_a,P_t_min_1)
    %get the covriance matrix of delta_Rt
    K_R=eye(4)-C_a*K_a;
    K_Q=C_a*(eye(8)-K_a*C_a);
    Sigma_rt= K_R*RY_a*K_R'+K_Q*(A_star*P_t_min_1*A_star'+Q_a)*K_Q';
end