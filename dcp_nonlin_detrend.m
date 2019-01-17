function nldtdat=dcp_nonlin_detrend(X,dtord)

% 
% Linear and nonlinear quadratic, cubic detrending
% 
% Inputs:
% 1. X is time x channel matrix
% 2. dtord is order (polynomial) of detrending
%   1=linear only; 2=linear and quadratic; 
%   3=linear, quadrtic and cubic; 4=all + quartic; 5=all + quintic;
% 
% Output:
% 1. nldtdat is non-linear detrended data
% 
% NOTES:
% 1. probably should always just set dtord=2
% 

lengtht=size(X,1);

% linear detrending
if dtord>=1
    X1=detrend(X);
%     subplot(dtord+1,1,1);plot(X1-X)
    X=X1;
end

% quadratic detrending
if dtord>=2
    t=(1:lengtht)';
    XX=[ones(lengtht,1) t t.^2];
    b=XX\X;
    nltrends=XX*b;
%     subplot(dtord+1,1,2);plot(nltrends)
    X=X-nltrends;
end

% cubic detrending
if dtord>=3
    XX=[ones(lengtht,1) t t.^2 t.^3];
    b=XX\X;
    nltrends=XX*b;
%     subplot(dtord+1,1,3);plot(nltrends)
    X=X-nltrends;
end

% quartic detrending
if dtord>=4
    XX=[ones(lengtht,1) t t.^2 t.^3 t.^4];
    b=XX\X;
    nltrends=XX*b;
%     subplot(dtord+1,1,4);plot(nltrends)
    X=X-nltrends;
end

% quintic detrending
if dtord>=5
    XX=[ones(lengtht,1) t t.^2 t.^3 t.^4 t.^5];
    b=XX\X;
    nltrends=XX*b;
%     subplot(dtord+1,1,5);plot(nltrends)
    X=X-nltrends;
end

nldtdat=X;
% subplot(dtord+1,1,6);plot(nldtdat)

end