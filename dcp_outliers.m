function [outrem,outremstruct]=dcp_outliers(X,fmeth,ometh,outliersin)

% 
% Inputs:
% 1. X is time x channel matrix; this is required
% 2. fmeth is string indicating preferred interpolation method; not required
%   NOTE: fmeth can be alot of things, but default is pchip
% 3. ometh is outlier detection method; also a string variable; not required
%   NOTE: ometh can be: median, mean, quartiles, grubbs, or gesd
% 4. outliersin is a vector of zeros and ones indicating known outliers; 
%     converted to a matrix internally if X is a matrix
% 
% Output:
% 1. outrem is data set x with outliers interpolated
% 2. outremstruct is structure with useful info about outliers and stuff
% 

if isempty(fmeth)
    fmeth='pchip';
end
if isempty(ometh)
    ometh='median';
end
% mometh=[];

if isempty(outliersin) || isequal(sum(outliersin),0)
    [outrem,TF,lowerbound,upperbound,centerbound]=filloutliers(X,fmeth,ometh);
end

if ~isempty(outliersin) && ~isequal(sum(outliersin),0)
    [outrem,TF,lowerbound,upperbound,centerbound]=filloutliers(X,fmeth,'OutLierLocations',logical(repmat(outliersin,1,size(X,2))));
    outremstruct.outliersin=find(outliersin==1);
end

outremstruct.fmeth=fmeth;
outremstruct.ometh=ometh;
outremstruct.outdetect=sum(TF);
outremstruct.TF=TF;
outremstruct.lowerbound=lowerbound;
outremstruct.upperbound=upperbound;
outremstruct.centerbound=centerbound;

end
