%==========================Calculate $\beta$===========================
%Parameters:
%	T, N, A, B, O, scale
%	scaleM: scale matrix calculated in Forward function
%Return:
%	beta: beta vector
%	prob: MLE probability
%Algorithm:
%Not scale
%	$\beta := \sum_{j=1}^{N} a_{i,j}b_j(o_{t+1})\beta_{t+1}(j), 1 \leq i \leq N, t = T-1,...,1 $ \\
%Scale
%	$\hat{\beta} = c_t\beta_t(i)$\\
%	$c_t = \frac{1}{\sum_{i=1}^{N}\alpha_t(i)}$\\
function beta=Backward(T, N, A, B, O)
	beta = ones(T, N);
	
	
	for(t=T-1:-1:1)
		for(i =1:N)
			sumb = 0.0;
			for(j=1:N)%use sum
				sumb =sumb+ A(i,j) * B(j, O(t+1)) * beta(t+1, j);
			end
			beta(t, i) = sumb;
% 			if(scale==1)
% 				beta(t, i) /= scaleM(t,1);
% 			end
		end
	end
	
% 	v = beta(1,:);
% 	prob = sum(v);
end