function [rzr,rzrt,rzrp]=dcp_rzr(corrval,dofused,doforig)

% 
% 
% Inputs:
% 1. corrval is pearson correlation value to be corrected
% 2. dofused is # of df removed from data (regressors used in denoising)
%     it may be necessary to enter a two value vector with dof from x and dof from y
%     for example [19 6] which would indicate that there were 19 regressors
%     for x time series denoising and 6 regressors for y time series denoising
%         Notenote: 19 and 6, in this example, may correspond to ROI mask
%         averages of the voxelwise df within each roi such that on
%         average, the # of noise regressors was 19 in roi x but actually there were 5 voxels with a varying # of df used;
%         if you want you could also use the max df of the voxels in roi x, etc., as you wish, and, as you can defend
% 3. doforig is # of df in original data (N of time series)
% 
% Notes:
%     1. DoF in this usage is probably not actually Dof as in statistics.
%     When computing the Pearson correlation, the DoF are never computed,
%     instead the N, or length of the time series is computed and used as
%     the denominator
%     2. The correction is as follows: convert to r to a fisher zcorr, then
%     calculate the test statistic with DoF consumed by denoising, then reconvert
%     to correlation value (un-fisher z transform) using the original DoF
% 

if length(dofused)>1
    dofused=mean(dofused);
end

DFnew=doforig-dofused;
zfc=atanh(corrval);
zfcT=zfc/(1/(sqrt(DFnew-3)));
corcorrected=zfcT/(sqrt(doforig-3));
rzr=tanh(corcorrected);

rzrt=(rzr.*sqrt(DFnew-2))./sqrt(1-rzr.^2);
rzrp=2*tcdf(abs(rzrt),(DFnew-2),'upper');


end
