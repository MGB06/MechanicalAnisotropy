close all; clc; clear; pack
% Insert the cordinate in pixel of ring artifact center
xc=2195; yc=1845;
% Setting the filter characteristc
decNum=4; wname='db2'; sigma=1; 
 projectdir = 'C:\original_microCT';  
dinfo = dir(fullfile(projectdir)); dinfo([dinfo.isdir]) = []; nfiles = length(dinfo);
for j=1:nfiles
  filename = fullfile(projectdir, dinfo(j).name);
  Im=double(imread(filename));
    % Image pre-processing
  [row_Im,col_Im]=size(Im); I=Im;
   deltah=min(abs(col_Im-xc),xc); deltav=min(abs(row_Im-yc),yc); deltal=min(deltah,deltav);
  fin_quad_dx=I(:,xc+deltal:col_Im); dx_lim=col_Im-(xc+deltal)+1
  fin_quad_sx=I(:,1:xc-deltal);   sx_lim=xc-deltal;
  fin_quad_dow=I(yc+deltal:row_Im,:);  dow_lim=row_Im-(yc+deltal)+1
  fin_quad_up=I(1:yc-deltal,:);  up_lim=yc-deltal;
  I(:,xc+deltal:col_Im)=[];   I(:,1:xc-deltal)=[];
  I(yc+deltal:row_Im,:)=[];   I(1:yc-deltal,:)=[];
  [row_try,col_try]=size(I);
%% Polar Transformation
 [row,col]=size(I);  cc=round(length(I)/2);
  Ip=polartrans(I,row,col,cc,cc);  I_pol=Ip';  
%% Filter application
  I_pol_f=real(xRemoveStripesVertical(I_pol,decNum,wname,sigma));
  [a,b]=size(I_pol_f);  I_pol_f(:,b)=[];  I_pol_f(a,:)=[];
   delta=I_pol-I_pol_f;   I_pol=I_pol-delta;
%% Cartesia Transfomration
  Ic=invpolar2(I_pol',size(I,2),size(I,1));
  [ra,ca]=size(I_pol_f);   [rb,cb]=size(Ic);
  r=abs(ra-rb);   c=abs(ca-cb);   rdelta=r/2;   cdelta=c/2;
  A=NaN(rdelta,rb);   B=NaN(ra,cdelta);  Ic=[A; Ic; A];  Ic=[B Ic B];
 for k=1:row
    for t=1:col
        TF =isnan(Ic(k,t));
        if TF == 1;
         Ic(k,t)=I(k,t);
        end
    end
end
  reco_sx=NaN(row_try,sx_lim);  reco_dx=NaN(row_try,dx_lim);  I_back=[reco_sx,Ic,reco_dx];
  [row_back,col_back]=size(I_back);  reco_dow=NaN(dow_lim,col_back);  reco_up=NaN(up_lim,col_back);
  I_back2=[reco_up;I_back;reco_dow];   diff_verifica=Im-I_back2;
for k=1:row_Im
    for t=1:col_Im
        TF =isnan(I_back2(k,t));
        if TF == 1;
         I_back2(k,t)=Im(k,t);
        end
    end
end 
  namef{j}= ['f_',dinfo(j).name];  %filtrata croppata
  namediff{j}= ['f_diff_',dinfo(j).name];  %filtrata croppata
  mshow(diff_verifica,[])
  imwrite(uint8(I_back2), namef{j},'bmp')
  imwrite(uint8(diff_verifica), namediff{j},'bmp')
end
 
%%% Functions Copyright (c) 2002 Peter Kovesi - http://www.csse.uwa.edu.au/
function pim = polartrans(im, nrad, ntheta, cx, cy, linlog, shape)
[rows, cols] = size(im);
if nargin==3        
    cx = cols/2+.5;  
    cy = rows/2+.5;
end
if nargin < 7, shape = 'full'; end
if nargin < 6, linlog = 'linear'; end
if strcmp(shape,'full')         
    dx = max([cx-1, cols-cx]);
    dy = max([cy-1, rows-cy]);
    rmax = sqrt(dx^2+dy^2);
elseif strcmp(shape,'valid')   
    rmax = min([cx-1, cols-cx, cy-1, rows-cy]);
else
    error('Invalid shape specification');
end
deltatheta = 2*pi/ntheta;
if strcmp(linlog,'linear')
    deltarad = rmax/(nrad-1);
    [theta, radius] = meshgrid([0:ntheta-1]*deltatheta, [0:nrad-1]*deltarad);    
elseif strcmp(linlog,'log')
    maxlogr = log(rmax);
    deltalogr = maxlogr/(nrad-1);    
    [theta, radius] = meshgrid([0:ntheta-1]*deltatheta, exp([0:nrad-1]*deltalogr));
else
    error('Invalid radial transformtion (must be linear or log)');
end
xi = radius.*cos(theta) + cx;  yi = radius.*sin(theta) + cy;  
[x,y] = meshgrid([1:cols],[1:rows]);
pim = interp2(x, y, double(im), xi, yi);
end
 
function im = invpolar2(pim,s1,s2)
r1 = sqrt(2*(s1-1)^2); r2 = sqrt(2*(s2-1)^2); rmax = max(r1,r2);
ntheta = size(pim,2); nrad = size(pim,1);
if mod(s1,2)==0
    rad0 = 0.5;
else
    rad0 = 0;
end
rad = linspace(rad0,rmax,nrad);
theta = linspace(0,2*pi,ntheta+1);
[R,T] = meshgrid(rad, theta);
[Xmat,Ymat] = meshgrid(linspace(1,(rmax*sqrt(2)),s1),linspace((rmax*sqrt(2)),1,s2));
Xcart = Xmat - max(Xmat(:)/2+0.5); Ycart = Ymat - max(Ymat(:)/2+0.5);
[invT,invR] = cart2pol(Xcart,Ycart);
ext_pim = ([pim, pim(:,end)])'; im_rot = interp2(R,T,ext_pim,invR,invT+pi,'*linear');
im = flipdim(im_rot,2); im( all( isnan( im ), 2 ), : ) = []; im( :, all( isnan( im ), 1 ) ) = []; 
end
 
function [nima]=xRemoveStripesVertical(ima,decNum,wname,sigma)
for ii=1:decNum
[ima,Ch{ii},Cv{ii},Cd{ii}]=dwt2(ima,wname);
end
for ii=1:decNum
fCv=fftshift(fft(Cv{ii}));
[my,mx]=size(fCv);
damp=1-exp(-[-floor(my/2):-floor(my/2)+my-1].^2/(2*sigma^2));
fCv=fCv.*repmat(damp',1,mx);
Cv{ii}=ifft(ifftshift(fCv));
end
nima=ima;
for ii=decNum:-1:1
nima=nima(1:size(Ch{ii},1),1:size(Ch{ii},2));
nima=idwt2(nima,Ch{ii},Cv{ii},Cd{ii},wname);
end
return
