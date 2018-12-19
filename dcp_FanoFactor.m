function FF=dcp_FanoFactor(x)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% DYnamic CONnectivity PROcessing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Inputs:
% 1. x is a time series of neural spiking/activity
% 
% Outputs:
% 1. FF is the Fano Factor
% 

vx=var(x);
mx=mean(x);
FF=vx./mx;

end