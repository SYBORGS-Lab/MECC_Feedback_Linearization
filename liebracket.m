function ad_fng = liebracket(f,g,x,n) 
% This function calcualtes the desired Lie Bracket of nth order given
% matrices f and g (of the same size), symbolic variable(s) x, and the
% order, n. 

% Lie Bracket Format: 
% [f,g](x)   = partial(g)/partial(x) * f(x) - partial(f)/partial(x) * g(x) 
% adf^k_g(x) = [f,adf^k-1_g](x) 

% Inputs: 
% f = symbolic matrix 
% g = symbolic matrix 
% x = symbolic vector of variables 
% n = order of the Lie Bracket 

% Output: 
% ad_fng = [ g   [f,g]  [f,[f,g]], ...]
%      n =   0     1        2      ...
% output = g     adf_g    adf2_g   ...

% Note: MATLAB's jacobian function calculates the partial derivatives of a
% matrix with respect to a specific variable 

% set up empty adjoint matrix to fill in 
ad_fng = sym(zeros(length(f),n+1)); 

% set the first column to be g 
ad_fng(:,1) = g; 


% fill in the rest of the adjoint based on n 
if n > 0 
    for k = 2:n+1 % starting at the second column 
        ad_fng(:,k) = jacobian(ad_fng(:,k-1),x) * f - jacobian(f,x) * ad_fng(:,k-1); 
    end 
end 

ad_fng = expand(ad_fng); 

end 