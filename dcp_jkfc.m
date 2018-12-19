function jkfc=dcp_jkfc(x)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% DYnamic CONnectivity PROcessing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 

    jkfc=zeros(size(x,1),size(x,2),size(x,2));
    for loop1=1:size(x,1)
        dat=x;
        dat(loop1,:)=[];
        jkfc(loop1,:,:)=corr(dat);
    end

end