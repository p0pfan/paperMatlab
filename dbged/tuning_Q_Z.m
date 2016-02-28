function [q_z,pz_infinate]=tuning_Q_Z(Q_z0,P_z_t,N)
    pz=zeros(4,4*N);
for i=1:N
    k=(P_z_t+Q_z0)/(P_z_t+Q_z0+eye(4));
    pz(:,4*i-3:4*i)=(eye(4)-k)*(P_z_t+Q_z0);
    P_z_t=pz(:,4*i-3:4*i);
    
end
    pz_infinate=pz(:,4*N-3:4*N);
    q_z=(eye(4)-pz_infinate)\(pz_infinate*pz_infinate);
end
