function squre_sigma_1=Sigma_one(MSE,sigma_rt,sqr_sigma_0)

    [V D]=eig(sigma_rt);
    temp_inv_D=sqrt(D)^-1;
    sigma_beta_zero_s_vector=V*temp_inv_D/V*sqrt(diag(MSE))*ones(4,1);
    a=sigma_beta_zero_s_vector.^2;
    squre_sigma_1=sqr_sigma_0+a;

end