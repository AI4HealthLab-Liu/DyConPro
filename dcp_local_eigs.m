function eigloc=dcp_local_eigs(x,nhood,ckeep,MASK)

% 
% Inputs:
% 1. x is 4D data (numbers, not BRIK/HEAD); you should already load it with BrikLoad.m
% 2. nhood is radius around center (not spherical)
% 3. ckeep is comps to retain
% 4. MASK is a binary mask and is required for the function
% 
% Output:
% 1. eigloc is reconstituted 4D volume with local svd filtering
% 

[xd,yd,zd,td]=size(x);
eigloc=zeros(xd,yd,zd,td);

for loop1=1+nhood:xd-nhood
    for loop2=1+nhood:yd-nhood
        for loop3=1+nhood:zd-nhood
            if MASK(loop1,loop2,loop3)~=0
                local_data=x(loop1-nhood:loop1+nhood,loop2-nhood:loop2+nhood,loop3-nhood:loop3+nhood,:);
                [u,~,~]=svd(local_data);
                eigloc(loop1,loop2,loop3,:)=sum(u(:,1:ckeep),2);      
            end
        end
    end
end


end