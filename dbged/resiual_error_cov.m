
function Sigma_rt=resiual_error_cov(C_a,A_star,RY_a,K_a,Q_a,P_t_min_1,Kerr)
    P_t_ba= diag([diag(eye(4)) ;diag(Kerr*eye(4))])*P_t_min_1*diag([diag(eye(4)) ;diag(Kerr*eye(4))]);
    %get the covriance matrix of delta_Rt
    K_R=eye(4)-C_a*K_a;
    K_Q=C_a*(eye(8)-K_a*C_a);
    Sigma_rt= K_R*RY_a*K_R'+K_Q*(A_star*P_t_ba*A_star'+Q_a)*K_Q';
end