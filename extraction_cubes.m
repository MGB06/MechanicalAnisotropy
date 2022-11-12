close all; clc; clear; pack
 loadiing_dir='U:\input_directory';
extracton_dir='U:\output_directory';
cd(loadiing_dir)
% Loading bone and mask volume
load MaskTest.mat; load uCTTest.mat 
% Specify image resolution in mm
resolution=0.02981; 
% Insert cube length side
cubic_length= 5; 
p=cubic_length/resolution; p=ceil(p);
[Vrows,Vcols,Vnfiles]=size(Vm);
ngridrow=floor(Vrows/p); ngridcol=floor(Vcols/p); 
nfacesh=ngridcol*ngridrow; nfacesv=floor(Vnfiles/p); 
ncubes=nfacesh*nfacesv;
% Identification of bone regions from mask segmentation
w=1;
for u=1:nfacesv   
s=1+(p*(u-1)); e=p*u;
esse=num2str(s);  eeee=num2str(e);  
       k=1; 
          for i=1:ngridrow-1
            for j=1:ngridcol-1
                    Vsub=Vm(1+(p*i-p):i*p,1+(p*j-p):j*p,s:e);
                    [a,b,c]=size(Vsub); Vnall(w)=a*b*c;
                    Vnbone(w)=max(size(find(Vsub)));
                    BVTV(w)=Vnbone(w)/Vnall(w);
              k=k+1; w=w+1;
            end 
          end
end
mBVTV=zeros(nfacesh,nfacesv,p); mVnall=zeros(nfacesh,nfacesv,p);
mVnbone=zeros(nfacesh,nfacesv,p); nlBVTV=max(size(mBVTV));
i=1;
for k=1:nfacesv
    for m=1:nfacesh
        for n=1:nfacesv
            mBVTV(m,n,k)=BVTV(i); mVnall(m,n,k)=Vnall(i); mVnbone(m,n,k)=Vnbone(i);         
        end
    end
end
%cube extraction algorithm  
w=1;
for u=1:nfacesv   
s=1+(p*(u-1)); e=p*u;
 
if (s >= 0) && (s <= 9)
    esse=['000',num2str(s)];  
elseif  (s >= 10) && (s <= 99)
    esse=['00',num2str(s)];  
elseif (s >= 100) && (s <= 999)
    esse=['0',num2str(s)];  
else   
    esse=num2str(s);  
end 
 
if (e >= 0) && (e <= 9)
    eeee=['000',num2str(e)];  
elseif (e >= 10) && (e <= 99)
    eeee=['00',num2str(e)];  
elseif (e >= 100) && (e <= 999)
    eeee=['0',num2str(e)];  
else   
    eeee=num2str(e);  
end
  k=1; 
    for i=1:ngridrow-1
        for j=1:ngridcol-1
              if BVTV(w)>0.01
                if (i >= 0) && (i <= 9)
                     iiii=['0',num2str(i)];  
                else
                     iiii=num2str(i);
                end
                 if (j >= 0) && (j <= 9)
                            jjjj=['0',num2str(j)];
                        else
                            jjjj=num2str(j);
                        end 
namefolder=[esse,'_',eeee,'__',iiii,'_',jjjj,];
            mkdir([extracton_dir,namefolder])
            cd([extracton_dir,namefolder]) 
            Vsub=Vb(1+(p*i-p):i*p,1+(p*j-p):j*p,s:e);
                            for t=1:p
                                    I(:,:)=Vsub(:,:,t);
                            	if (t >= 0) && (t <= 9)
                                    tttt=['0',num2str(t)];  
                              else
                                    tttt=num2str(t);
                              end
                                namef=['cube',tttt,'.tif']; 
 		                    imwrite(uint8(I), namef)
                                k=k+1; t=t-100;
                            end 
                    end
                    w=w+1;
            end
        end
end
