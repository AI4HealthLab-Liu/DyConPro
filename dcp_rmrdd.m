% function gRMRDD=dcp_rmrdd(X)


% Residual Motion Related Distance Dependence
% Currently in script format, but soon to be funcitonal
% motion is loaded via .mat files, or however you want if you code it
% MRI data is loaded via afni_matlab BrikLoad function

clear,clc

atlas_data='/scratch/mtobia/COBRE_derivatives/fmriprep/Power_seed_masks+tlrc.BRIK';
[~,ATLAS,~,~]=BrikLoad(atlas_data);
power_coords=load('/scratch/mtobia/COBRE_derivatives/fmriprep/Power_236_Coords.1D');

[ar,ac]=size(power_coords);
for loop1=1:ar
    for loop2=1:ar
        distances(loop1,loop2)=norm(power_coords(loop1,1:3)-power_coords(loop2,1:3));
    end
end
D=dcp_2Dto1D(distances);

wdir='/scratch/mtobia/COBRE_derivatives/fmriprep/sub-';
sns={'nostep','step','noproc'};

for snsloop=1:size(sns,2)

    for sub=1:148
        cd([wdir,int2str(sub),'/func'])
        pwd
        if snsloop==1
            load(['sub-',int2str(sub),'_',sns{snsloop},'_p36_outliers.mat']);
            gMo(sub)=outstrc.fd_avg;
            gMo(sub)=outstrc.fd_avg;
            data=([wdir,int2str(sub),'/func/sub-',int2str(sub),'_',sns{snsloop},'_p36_outliers_denoised+tlrc.BRIK']);
            roi_tcs=dcp_roi_timecourse_extract(data,atlas_data,1,[],[]);
            [fc,pp]=corr(roi_tcs);
            gFC(:,sub)=atanh(dcp_2Dto1D(fc));
        end
        if snsloop==2
            load(['sub-',int2str(sub),'_',sns{snsloop},'_p36_compreg_outliers.mat']);
            gMo(sub)=outstrc.fd_avg;
            data=([wdir,int2str(sub),'/func/sub-',int2str(sub),'_',sns{snsloop},'_p36_compreg_outliers_denoised+tlrc.BRIK']);
            roi_tcs=dcp_roi_timecourse_extract(data,atlas_data,1,[],[]);
            [fc,pp]=corr(roi_tcs);
            gFC(:,sub)=atanh(dcp_2Dto1D(fc));
        end
        if snsloop==3
            data=([wdir,int2str(sub),'/func/sub-',int2str(sub),'_scaled+tlrc.BRIK']);                        
            roi_tcs=dcp_roi_timecourse_extract(data,atlas_data,1,[],1);
            [fc,pp]=corr(roi_tcs);
            gFC(:,sub)=atanh(dcp_2Dto1D(fc));
        end
    end
    
    xx=size(gFC,1);
    for loop1=1:xx
        [mocorr(loop1,snsloop),mopval(loop1,snsloop)]=corr(real(gFC(loop1,:))',gMo');
    end

%     mocorr=atanh(mocorr);
    [fit1{snsloop},gof1{snsloop},obs1{snsloop}]=fit(D,mocorr(:,snsloop),'poly1','Robust','Off','Normalize','Off');
    [fit2{snsloop},gof2{snsloop},obs2{snsloop}]=fit(D,mocorr(:,snsloop),'exp2','Robust','Off','Normalize','Off');
    [fit3{snsloop},gof3{snsloop},obs3{snsloop}]=fit(D,mocorr(:,snsloop),'poly2','Robust','Off','Normalize','Off');
    figure();subplot(1,2,1);
    scatplot(D,mocorr(:,snsloop));colorbar('off')
    ylim([-1 1])
    hold on;plot(fit1{snsloop},'r')
    hold on;plot(fit2{snsloop},'g')
    hold on;plot(fit3{snsloop},'y')
    title(sns{snsloop})
    subplot(1,2,2);
    histogram(mocorr(:,snsloop),45)
    title(sns{snsloop})
    

end

[rr1,pp1]=corr(D,mocorr(:,1)) 
[rr2,pp2]=corr(D,mocorr(:,2)) 
[rr3,pp3]=corr(D,mocorr(:,3))


