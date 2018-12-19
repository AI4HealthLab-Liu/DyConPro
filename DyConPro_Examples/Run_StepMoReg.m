

clear,clc

sub=1:1:148; 

% Motion is: roll pitch yaw Zd Yd Xd

% Unused Opts
% Opts.morexs=; %additional params to regress out from data; path to textfile name
% Opts.swthr=0; %apply bonferonni to stepwise threshold; 0=no, 1=yes

for loop1=1 %:148 %length(sub)
    
    Opts.prefix=['subj_',int2str(sub(loop1)),'_nostep_mo6_fwhmxyz']; %name of output BRIK file
    Opts.data=['/Users/mtobia/Desktop/COBRE/subj_',int2str(sub(loop1)),'/subj_',int2str(sub(loop1)),'.volreg+orig.BRIK'];
    Opts.mask=['/Users/mtobia/Desktop/COBRE/subj_',int2str(sub(loop1)),'/subj_',int2str(sub(loop1)),'.epimaskdi+orig.BRIK'];
    Opts.fwhm_xyz=['/Users/mtobia/Desktop/COBRE/subj_',int2str(sub(loop1)),'/subj_',int2str(sub(loop1)),'.localfwhm_xyz+orig.BRIK'];
    Opts.moparams=['/Users/mtobia/Desktop/COBRE/subj_',int2str(sub(loop1)),'/subj_',int2str(sub(loop1)),'.motion.1D'];
%     Opts.tissuemask=['/Users/mtobia/Desktop/COBRE/subj_',int2str(sub(loop1)),'/subj_',int2str(sub(loop1)),'.tissuemasks+orig.BRIK'];
%     Opts.edgemask=['/Users/mtobia/Desktop/COBRE/subj_',int2str(sub(loop1)),'/subj_',int2str(sub(loop1)),'.EDGEer+orig.BRIK']; %path to edgemask file
%     Opts.edgenum=24; %this is number of comps to retain as edgemask regs
    Opts.brikout=1; %write out a BRIK file? 0=no, 1=yes
    Opts.framewised=1; %compute framewise displacement? 0=no, 1=yes
    Opts.F24=0; %compute Friston24? 0=no, 1=yes
    Opts.fd_thr=.35; %threshold for framewise displacement scrubbing
    Opts.nodenoise=0; %Only output motion metrics & Xvars? 0=no, 1=yes
    Opts.view='+orig'; %can be +tlrc or +orig
    Opts.gsr=0; %global signal regression? 0=no, 1=yes
    Opts.remean=0; %remean denoised data? 0=no, 1=yes
    Opts.nostep=1; %perform regular regression instead of stepwise; 0=no, 1=yes
    Opts.voxcorr=0; %calculate global correlation values for all voxels; 0=no, 1=yes
    Opts.compregs=[]; %reduce the full set of all regressors to just n princomps; []=no, n=yes
%     Opts.Fs=.5; %use this Fs for bandpass filtering the nuisregs; exist=yes
%     Opts.fbands=[.01 .1]; %this is the bandpass filter for nuisregs; exist=yes    
%     Opts.xvarsnorm=1; %this norms the nuisregs prior to regression; 0=no, 1=yes (does not norm the fwhm_xyz regressors)

    cd(['/Users/mtobia/Desktop/COBRE/subj_',int2str(sub(loop1))]);
    
    if exist(Opts.data,'file')
        [BRIK_den,outstrc]=dcp_stepmoreg(Opts);
    end   
    
    save(['subj_',int2str(loop1),'_nostep_mo6_fwhmxyz.mat'],'outstrc','-mat')
    
end


