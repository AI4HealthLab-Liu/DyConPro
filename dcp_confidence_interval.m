function [CI,ts,SEM]=dcp_confidence_interval(QW)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% DYnamic CONnectivity PROcessing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% QW is vector of DCS values (or other dFC values)

SEM = std(QW)/sqrt(length(QW));               % Standard Error
ts = tinv([0.025  0.975],length(QW)-1);      % T-Score
CI = mean(QW) + ts*SEM;   % Confidence Intervals

end
