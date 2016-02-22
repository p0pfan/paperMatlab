function[Q_z_t,Q_beta_t]=cal_Qz_Qbeta(Q_z_0,Z_tip_beta,temp_P_z_t,sigma_infinate)
    Q_z_adjust=diag(diag(Z_tip_beta*Z_tip_beta'-temp_P_z_t));
    Q_z_t=max(Q_z_0,0.2*Q_z_adjust);
    Q_beta_t=diag(diag(sigma_infinate))*max(Q_z_0,10*Q_z_adjust);
end