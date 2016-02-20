%==========================Calculate $\alpha$===========================
%Parameters:
%	T, N, pi, A, B, O, 
%	scale: if scale = 1, scale the calculation
%Return:
%	alpha: alpha vector
%	prob: MLE probability
%	scaleM: $\sum_{i=1}^{N}\alpha_t(i)$ to be used for scaling beta
%Algorithm:
%Not scale:
%	$\alpha_1(i) = \pi_ib_i(o_1)$\\
%	$\alpha_{t+1}(j) := [\sum_{i=1}^{N} \alpha_t(i)a_{i,j}]b_j(o_{t+1}), 1 \leq t\leq T-1$\\
%Scale:
%	$\hat{\alpha}_t(i) = \frac{\alpha_t(i)}{\sum_{i=1}^{N}\alpha_t(i)}$
function alpha=Forward(T, N, pi, A, B, O)
	alpha = zeros(T, N);
% 	scaleM = zeros(T, 1);
	for(i=1:N)
		alpha(1, i) = pi(i) * B(i, i);
        %according to the number of the state number
        
% 		scaleM(1) += alpha(1, i);
	end
	
% 	if(scale==1)
% 		alpha(1, :) /= scaleM(1);
% 	end
	for(t=1:T-1)
		for(j=1:N)
			suma = 0.0;
			for(i=1:N)
				suma =suma+ alpha(t, i) * A(i, j);
			end
			alpha(t+1, j) = suma * B(j, O(t+1));
% 			
		end
% 		if(scale==1)
% 			alpha(t+1, :) /= scaleM(t+1);
% 		end
	end
end