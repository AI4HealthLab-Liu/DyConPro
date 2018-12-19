
clear,clc

wdir='/Users/mtobia/Desktop/COBRE/';

hemiseed={'R','L'};
suffix={'_HF','_MF','_LF'};
freqfs=[.11 .15 .5;.06 .1 .5;.01 .05 .5];

for hemiloop=1:2
    for floop=1:3
        for loop1=1:25

            cd([wdir,'subj_',int2str(loop1)]);

            S=([wdir,'piri-',hemiseed{hemiloop},'-seed-mask+tlrc.BRIK']);
            X=([wdir,'subj_',int2str(loop1),'/subj_',int2str(loop1),'_step_gsr_f24_wmcsf_edge24_comp36_denoised_mni_6mm.nii']);
            M=([wdir,'MNI_2009_3mm_mask+tlrc.BRIK']);

            prefix=['subj_',int2str(loop1),'_',hemiseed{hemiloop},'_piri'];
            outview='+tlrc';

            [ifc,ifc_std,ifc_map,ifc_dfcmap,outstrc]=dcp_fast_s2b_ifc(S,X,M,freqfs(floop,:),prefix,outview,suffix{floop},1);

            save([prefix,suffix{floop},'_ifc_s2b_confidence_intervals.mat'],'outstrc','-mat')

        end
    end  
end
    
    