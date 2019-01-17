function swfc=dcp_swfc(x,winlength)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% X is matrix with time X region
% winlength is size of sliding window

if mod(winlength,2)~=0 %if winlength is odd number then add 1 to make it even
    winlength=winlength+1;
end

x=[flipud(x(1:winlength/2,:)) ; x ; flipud(x(end-(winlength/2)+1:end,:))];

[row,col]=size(x);
wl=winlength-1;
windows=row-winlength;
swfc=zeros(windows,col,col);
for win=1:windows    
    swfc(win,:,:)=corr(x(win:win+wl,:));   
end
swfc=dcp_mat2tens(dcp_ten2mat(swfc));


end