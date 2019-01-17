function xrescaled=dcp_median_rescale(X,type)


% Inputs:
% 1. X is time x channel matrix
% 2. type is for median or mean scaling
%       type==1 is median scaling
%       type==2 is mean scaling
% Output:
% 1. xrescaled is data rescaled to central tendency of 1000

% clear,clc
% x=randi([375 800],150,35);

% xrescaled=X+(1000-median(X,2));

if isempty(type)
    type=2;
end

if type==1
    xrescaled=(X./median(X,2)).*1000;
end

if type==2
    xrescaled=(X./(mean(X,2))).*1000;
end

end