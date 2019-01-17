function dfc=dcp_dfc_rand(td,x,type)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% DYnamic CONnectivity PROcessing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Inputs:
% 1. td is duration of dfc
% 2. x is number of cols or rows (symmetric)
% 3. type isempty or ~exist for randn dfc, or vector indicating how many networks
%     and correlation level of each, e.g., [3 .5] for 3 networks with a
%     intranetwork correlation of .5
% 
% NOTES:
% 1. you can use fc_snapshot to inspect your dfc and scrool through time;
%     just type fc_snapshot at the command line, then click on 'Get FC' to load
%     the dfc data from the workspace; then scroll through time!; also can type
%     mdfc=squeeze(mean(dfc,1)); at command line to inspect simulated time averaged dfc
%     networks also in fc_snapshot (so crude, but so useful)
% 

if ~exist('type','var') || isempty(type)
    dfc=dcp_mat2tens(dcp_ten2mat(randn(td,x,x)));
end

if exist('type','var') && ~isempty(type) && length(type)==2
    nets=[];
    for loop1=1:type(1)
        sigs=dcp_GenCorrTS(td,x,type(2));
        nets=[nets sigs];
    end
    dfc=dcp_dcs(nets);
end
    

end