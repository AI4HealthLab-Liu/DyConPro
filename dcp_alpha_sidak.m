function alpha_sidak=dcp_alpha_sidak(alpha,m)
 
% 
% Inputs:
% 1. alpha is regular alpha value, i.e., .05
% 2. m is number tests for which to correct
% 

alpha_sidak=1-((1-alpha).^(1/m));

end
