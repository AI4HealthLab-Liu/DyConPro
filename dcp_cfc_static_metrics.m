function CFC=dcp_cfc_static_metrics(lfo,hfo)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% DYnamic CONnectivity PROcessing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Inputs:
% 1. lfo is low frequency oscillation
% 2. hfo is high frequency osicllation
% 
% Outputs:
% 1. CFC is a structure with multiple static CFC metrics
% 
% Notes:
% 1. histogram based metrics use 20 phase bins

phibins=-pi:(pi--pi)/20:pi;
phi_lfo=angle(hilbert(lfo));
amp_hfo=abs(hilbert(hfo));
amp_hfo_phi=angle(hilbert(zscore(amp_hfo)));
phi_hfo=angle(hilbert(hfo));

cplx=[amp_hfo(:) phi_lfo(:)];
sc=sortrows(cplx,2);

for loop1=1:length(phibins)-1
    binned(loop1)=mean(amp_hfo(phi_lfo>=phibins(loop1) & phi_lfo<phibins(loop1+1)));
    binsum(loop1)=sum(amp_hfo(phi_lfo>=phibins(loop1) & phi_lfo<phibins(loop1+1)));
    bincounts(loop1)=length(amp_hfo(phi_lfo>=phibins(loop1) & phi_lfo<phibins(loop1+1)));
    this(loop1,:)=[phibins(loop1) phibins(loop1+1)];
end
binval=binned./sum(binned);
Z=zscore(amp_hfo).*exp(j.*phi_lfo); % complex pac signal; plot(Z) to inspect for PAC
HP=-sum(binval.*log(binval));
% Order=1-(HP./log(length(binval))); % this is 1-percentropy; higher values are higher cfc strength

pli=(abs(sum(sign(phi_lfo-amp_hfo_phi))))/length(phi_lfo);
plv=(abs(sum(cos(phi_lfo-amp_hfo_phi))))/length(phi_lfo);

% Miller measures nested multi-freq signals, not AMPAC
MI_canolty=abs(mean(Z));
MI_schnitzler=(MI_canolty/sqrt(mean(cplx(:,1).^2))).*(1/sqrt(length(lfo)));
MI_cohen=abs(mean(cos(phi_lfo-amp_hfo_phi)));
MI_tort=(log(length(phibins)-1)-HP)./log(length(phibins)-1); % Same as Order
MI_miller=((abs(hilbert(lfo))+amp_hfo)/2).*cos(phi_lfo-phi_hfo);
W=1-abs((abs(hilbert(lfo))-amp_hfo)./(abs(hilbert(lfo))+amp_hfo));
MI_tobia=W.*cos(phi_lfo-amp_hfo_phi);
ue=dcp_getspline(MI_miller',[]);le=-dcp_getspline(-MI_miller',[]);
MI_miller_env=(ue+le)/2;
AM_ratio=(max(binval)-min(binval))/max(binval);
AM2_ratio=(max(binval)-min(binval))/(max(binval)+min(binval));

CFC.MI_canolty=MI_canolty;
CFC.MI_shcnitzler=MI_schnitzler;
CFC.MI_cohen=MI_cohen;
CFC.MI_tort=MI_tort;
CFC.MI_miller=mean(MI_miller);
CFC.MI_tobia=mean(MI_tobia);
CFC.AM_ratio=AM_ratio;
CFC.AM2_ratio=AM2_ratio;
CFC.plv=plv;
CFC.pli=pli;
CFC.dMiller=MI_miller;
CFC.dMillerE=MI_miller_env';
CFC.pbins=mean(this,2);
CFC.binval=binval;


end