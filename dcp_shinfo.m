function shinfo=dcp_shinfo(x)

[row,c1,c2]=size(x);
xnew=zeros(row,c1,c2);        
dx=[zeros(1,c1,c2);diff(x)];
xnew(x>0 & dx>0)=2;
xnew(x>0 & dx<0)=1;
xnew(x<0 & dx>0)=-1;
xnew(x<0 & dx<0)=-2;

shinfo=zeros(1,row);
for win=1:row
    f=reshape(triu(squeeze(xnew(win,:,:))),1,c1*c2);
    f(f==0)=[];
    e1=length(find(f==-2))/length(f);
    e2=length(find(f==-1))/length(f);
    e3=length(find(f==1))/length(f);
    e4=length(find(f==2))/length(f);
    ee1=e1*log2(e1);ee1(isnan(ee1)==1)=0;
    ee2=e2*log2(e2);ee2(isnan(ee2)==1)=0;
    ee3=e3*log2(e3);ee3(isnan(ee3)==1)=0;
    ee4=e4*log2(e4);ee4(isnan(ee4)==1)=0;
    shinfo(win)=-sum([ee1 ee2 ee3 ee4]);
end
shinfo(1)=mean(shinfo);

end