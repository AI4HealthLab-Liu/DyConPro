function snic=dcp_snic(x1,x2,fps)
%
%
% 
% USAGE: find changes in power spectrum due to denoising (or some other
%   manipulation of the time series)
% 
% Inputs:
% 1. x1 is original time series
% 2. x2 is denoised time series
% 3. fps=[lf,hf,fs]
% 
% Outputs:
% 1. snic is structure with info about change in power spectrum
%     snic.snic=total percent new power injected due to denoising
%     snic.corr=correlation of x1 and x2 spectra
%     snic.percpow=[total power x1 total power x2 x2/x1]
%     snic.powerup=frequencies that increased power after denoising
%     snic.powerdown=frequencies that decreased power after denoising
% 
% NOTE:
%     1. snic = Spectral Noise Injection Coefficient
%     2. power spectral density of x2 is normailized by that of x1 to obtain
%         percent power at each frequency; then pow2-pow1 to get percent
%         power change; then negative changes are zeroed
%

[z1,pts1]=dcp_mkfreq(x1-mean(x1),fps(1),fps(2),fps(3));
[~,pts2]=dcp_mkfreq(x2-mean(x2),fps(1),fps(2),fps(3));

ptsdif=(pts2'./sum(pts1))-(pts1'./sum(pts1));
powerdown=z1(ptsdif<0);
ptsdif(ptsdif<0)=0;

% figure();
% subplot(3,1,1);plot(z1,pts1./sum(pts1)) %plot noisy spectrum
% subplot(3,1,2);plot(z1',[pts1'./sum(pts1) pts2'./sum(pts1)]) %plot denoised spectrum
% subplot(3,1,3);plot(z1,ptsdif) % plot denoised spectrum percentage of noisy spectrum

snic.snic=sum(ptsdif);
snic.corr=corr(pts1',pts2');
snic.percpow=[sum(pts1) sum(pts2) sum(pts2)/sum(pts1) 1-sum(pts2)/sum(pts1)];
snic.powerup=z1(ptsdif>0);
snic.powerdown=powerdown;
snic.ptsdif=ptsdif;
snic.z=z1';

end
