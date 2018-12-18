function dfc=dcp_dfc_rand(td,x)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% DYnamic CONnectivity PROcessing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Inputs:
% 1. td is duration of dfc
% 2. x is number of cols or rows (symmetric)

dfc=dcp_mat2tens(dcp_ten2mat(randn(td,x,x)));

end