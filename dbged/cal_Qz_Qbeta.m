function[Q_z_t,Q_beta_t]=cal_Qz_Qbeta(Q_z_0,Z_tip_beta,sigma_infinate,z_c)
    K_p=Q_z_0*(diag(diag(z_c*z_c')))^-1;
    Q_z_adjust=K_p*diag(diag(Z_tip_beta*Z_tip_beta'));
    Q_z_t=max(Q_z_0,0.2*Q_z_adjust);
    Q_beta_t=diag(diag(sigma_infinate))*max(Q_z_0,10*Q_z_adjust);
 
%     Q_z_t=Q_z_0;
%     Q_beta_t=100*diag(diag(sigma_infinate))*Q_z_0
end