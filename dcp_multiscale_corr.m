function ms_corr=dcp_multiscale_corr(x,y,scales)

% 
% Inputs:
% 1. x is time series
% 2. y is time series
% 3. scales is number of integration scales to investigate
% 
% clear,clc
% % x=randn(2500,1);y=randn(2500,1);
% xy=dcp_GenCorrTS(900,2,.35);
% x=xy(:,1);
% y=xy(:,2);

if ~exist('scales','var') || isempty(scales)
    scales=10;
end

ms_corr=zeros(1,scales);
ms_pval=zeros(1,scales);
endfix=zeros(1,scales);

for loop1=1:scales
        
    endfix(loop1)=loop1*floor(size(x,1)/loop1);   
    xs=mean(reshape(x(1:endfix(loop1),:),loop1,endfix(loop1)/loop1),1);
    ys=mean(reshape(y(1:endfix(loop1),:),loop1,endfix(loop1)/loop1),1);
    [rr,pp]=corr(xs',ys');
    ms_corr(loop1)=rr;
    ms_pval(loop1)=pp;
        
%     figure(1);subplot(5,scales/5,loop1)
%     plot(xs);hold on;plot(ys)
%     xlim([0 length(xs)])

end
    
% figure();bar(ms_corr);ylim([-1 1])
% figure();bar(1-ms_pval);ylim([.9 1])
    
end
