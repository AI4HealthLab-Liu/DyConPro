function sigs_filt=dcp_cfc_filt(data,style,Fs)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% DYnamic CONnectivity PROcessing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Inputs:
% 1. data is a time x channel matrix
% 2. style is a string of either EB (equal bandwidth) or BB (Buzsaki bandwidths)
% 
% Outputs:
% 1. sigs_filt is a concatenated matrix of time x (channel x fbands)
% 
% Notes:
% 1.

% These are based on Buzsaki 2004 fband bins
if strcmp(style,'BB')
    slow2=[.198 .24];
    slow3=[.073 .198];
    slow3a=[.135 .198];
    slow3b=[.073 .135];
    slow4=[.027 .073];
    slow5=[.01 .027];
    buzbands=[slow3b;slow4;slow5];
    BB=[];
    for loop1=1:size(buzbands,1)
        F=dcp_buttfilt(data,4,buzbands(loop1,1),buzbands(loop1,2),Fs);
        BB=[BB F];
    end
    sigs_filt=BB;
end

if strcmp(style,'EB')
% These are equally spaced bandwidth bins
    band1=[.01 .05];
    band2=[.06 .1];
    band3=[.11 .15];
    band4=[.16 .2];
    band5=[.21 .24];
    eqbands=[band3;band2;band1];
    EB=[];
    for loop1=1:size(eqbands,1)
        F=dcp_buttfilt(data,4,eqbands(loop1,1),eqbands(loop1,2),Fs);
        EB=[EB F];
    end
    sigs_filt=EB;
end

end









