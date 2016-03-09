function [bool_value,P_ht,P_zt_ht]=cal_HMM(sqr_sigma_0,sqr_sigma_1,pi,Z_tip_beta)
    N=size(sqr_sigma_0,1); %or N=size(sqr_sigma_1,1)
    bool_value=zeros(N,1);
    P_ht=zeros(2,N);
    A_CP=[0.999 0.001;0.001 0.999];
    
    for i=1:N
        h0_likelihood_value=1/sqrt(2*3.14*sqr_sigma_0(i,1))*exp(-0.5*Z_tip_beta(i,1)/sqr_sigma_0(i,1));
        h1_likelihood_value=1/sqrt(2*3.14*sqr_sigma_1(i,1))*exp(-0.5*Z_tip_beta(i,1)/sqr_sigma_1(i,1));
        C_CP=[h0_likelihood_value,0;0,h1_likelihood_value]
        
        P_ht(:,i)=A_CP*pi(:,i);
        P_zt_ht=C_CP*P_ht(:,i);
        state_probability=P_zt_ht/(P_zt_ht(1,1)+P_zt_ht(2,1));
        
        if(state_probability(1,1)>state_probability(2,1))
            bool_value(i,1)=0;
        else
            bool_value(i,1)=1;
        end  
    end
    

end