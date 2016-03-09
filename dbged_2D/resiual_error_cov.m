function Sigma_rt=resiual_error_cov(C_a,A_star,RY_a,K_a,Q_a,P_t_min_1,Kerr)
    P_t_bar= diag([diag(eye(2)) ;diag(Kerr\eye(2))])*P_t_min_1*diag([diag(eye(2)) ;diag(Kerr\eye(2))]);
    Q_a_bar=diag([1 1  0 0])*Q_a*diag([1 1 0 0]);
    %get the covriance matrix of delta_Rt
    K_R=eye(2)-C_a*K_a;
    K_Q=C_a*(eye(4)-K_a*C_a);
    Sigma_rt= K_R*RY_a*K_R'+K_Q*(A_star*P_t_bar*A_star'+Q_a_bar)*K_Q';
end