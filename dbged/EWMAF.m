function [next_Z_tip_beta,temp_P_z_t]=EWMAF(Q_z,P_z_t,Z,temp_Z_tip_beta)

    I=eye(4);   
 
    K_EWMAF=(P_z_t+Q_z)/(P_z_t+Q_z+I);
    
    next_Z_tip_beta=temp_Z_tip_beta+ K_EWMAF*(Z-temp_Z_tip_beta);
    
    temp_P_z_t=(I-K_EWMAF)*(P_z_t+Q_z);
   
end