clc; clear; close all
projectdir = 'U:\input_files'; cd(projectdir)
dinfo_mainfolder = dir(fullfile(projectdir));
 
volumeRVE=0.000125;
dir='U:\MATLAB_codes\'; cd(dir)
[loadingfile] = inputfile_stringsV7;
 
dirout='U:\output_files’’
 
for i=3:length(dinfo_mainfolder)
    try
uFEM.name{i-2}=string(dinfo_mainfolder(i).name);
 
resultsdir = [dinfo_mainfolder(i).folder,'\',dinfo_mainfolder(i).name,'\']; cd (resultsdir);
geometric_parameters.lo=0.005; geometric_parameters.Ao=0.000025;
 
%% Loading numerical outcomes for each of six BCs simulted
load Compression_X; load Compression_Y; load Compression_Z;
load Shear_XY; load_Shear_XZ; load_Shear_YZ;
 
%% Numerical outcomes processing
EE11=(sum((CompressionX.EPTOX).*CompressionX.volume*1000))/volumeRVE;
EE21=(sum((CompressionX.EPTOY).*CompressionX.volume*1000))/volumeRVE;
EE31=(sum((CompressionX.EPTOZ).*CompressionX.volume*1000))/volumeRVE;
EE41=(sum((CompressionX.EPTOXY).*CompressionX.volume*1000))/volumeRVE;
EE51=(sum((CompressionX.EPTOYZ).*CompressionX.volume*1000))/volumeRVE;
EE61=(sum((CompressionX.EPTOXZ).*CompressionX.volume*1000))/volumeRVE;
EE12=(sum((CompressionY.EPTOX).*CompressionY.volume*1000))/volumeRVE;
EE22=(sum((CompressionY.EPTOY).*CompressionY.volume*1000))/volumeRVE;
EE32=(sum((CompressionY.EPTOZ).*CompressionY.volume*1000))/volumeRVE ;
EE42=(sum((CompressionY.EPTOXY).*CompressionY.volume*1000))/volumeRVE ;
EE52=(sum((CompressionY.EPTOYZ).*CompressionY.volume*1000))/volumeRVE ;
EE62=(sum((CompressionY.EPTOXZ).*CompressionY.volume*1000))/volumeRVE ;
EE13=(sum((CompressionZ.EPTOX).*CompressionZ.volume*1000))/volumeRVE ;
EE23=(sum((CompressionZ.EPTOY).*CompressionZ.volume*1000))/volumeRVE ;
EE33=(sum((CompressionZ.EPTOZ).*CompressionZ.volume*1000))/volumeRVE ;
EE43=(sum((CompressionZ.EPTOXY).*CompressionZ.volume*1000))/volumeRVE ;
EE53=(sum((CompressionZ.EPTOYZ).*CompressionZ.volume*1000))/volumeRVE ;
EE63=(sum((CompressionZ.EPTOXZ).*CompressionZ.volume*1000))/volumeRVE ;
EE14=(sum((ShearXY.EPTOX).*ShearXY.volume*1000))/volumeRVE ;
EE24=(sum((ShearXY.EPTOY).*ShearXY.volume*1000))/volumeRVE ;
EE34=(sum((ShearXY.EPTOZ).*ShearXY.volume*1000))/volumeRVE ;
EE44=(sum((ShearXY.EPTOXY).*ShearXY.volume*1000))/volumeRVE; 
EE54=(sum((ShearXY.EPTOYZ).*ShearXY.volume*1000))/volumeRVE ;
EE64=(sum((ShearXY.EPTOXZ).*ShearXY.volume*1000))/volumeRVE ;
EE15=(sum((ShearYZ.EPTOX).*ShearYZ.volume*1000))/volumeRVE ;
EE25=(sum((ShearYZ.EPTOY).*ShearYZ.volume*1000))/volumeRVE ;
EE35=(sum((ShearYZ.EPTOZ).*ShearYZ.volume*1000))/volumeRVE ;
EE45=(sum((ShearYZ.EPTOXY).*ShearYZ.volume*1000))/volumeRVE ;
EE55=(sum((ShearYZ.EPTOYZ).*ShearYZ.volume*1000))/volumeRVE ;
EE65=(sum((ShearYZ.EPTOXZ).*ShearYZ.volume*1000))/volumeRVE ;
EE16=(sum((ShearXZ.EPTOX).*ShearXZ.volume*1000))/volumeRVE ;
EE26=(sum((ShearXZ.EPTOY).*ShearXZ.volume*1000))/volumeRVE ;
EE36=(sum((ShearXZ.EPTOZ).*ShearXZ.volume*1000))/volumeRVE;
EE46=(sum((ShearXZ.EPTOXY).*ShearXZ.volume*1000))/volumeRVE;
EE56=(sum((ShearXZ.EPTOYZ).*ShearXZ.volume*1000))/volumeRVE;
EE66=(sum((ShearXZ.EPTOXZ).*ShearXZ.volume*1000))/volumeRVE;
 
SS11=((sum((CompressionX.SX).*CompressionX.volume*1000))/volumeRVE)/10^6;
SS21=((sum((CompressionX.SY).*CompressionX.volume*1000))/volumeRVE)/10^6;
SS31=((sum((CompressionX.SZ).*CompressionX.volume*1000))/volumeRVE)/10^6;
SS41=((sum((CompressionX.SXY).*CompressionX.volume*1000))/volumeRVE)/10^6;
SS51=((sum((CompressionX.SYZ).*CompressionX.volume*1000))/volumeRVE)/10^6;
SS61=((sum((CompressionX.SXZ).*CompressionX.volume*1000))/volumeRVE)/10^6;
SS12=((sum((CompressionY.SX).*CompressionY.volume*1000))/volumeRVE)/10^6;
SS22=((sum((CompressionY.SY).*CompressionY.volume*1000))/volumeRVE)/10^6;
SS32=((sum((CompressionY.SZ).*CompressionY.volume*1000))/volumeRVE)/10^6;
SS42=((sum((CompressionY.SXY).*CompressionY.volume*1000))/volumeRVE)/10^6;
SS52=((sum((CompressionY.SYZ).*CompressionY.volume*1000))/volumeRVE)/10^6;
SS62=((sum((CompressionY.SXZ).*CompressionY.volume*1000))/volumeRVE)/10^6;
SS13=((sum((CompressionZ.SX).*CompressionZ.volume*1000))/volumeRVE)/10^6;
SS23=((sum((CompressionZ.SY).*CompressionZ.volume*1000))/volumeRVE)/10^6;
SS33=((sum((CompressionZ.SZ).*CompressionZ.volume*1000))/volumeRVE)/10^6;
SS43=((sum((CompressionZ.SXY).*CompressionZ.volume*1000))/volumeRVE)/10^6;
SS53=((sum((CompressionZ.SYZ).*CompressionZ.volume*1000))/volumeRVE)/10^6;
SS63=((sum((CompressionZ.SXZ).*CompressionZ.volume*1000))/volumeRVE)/10^6;
 
SS14=((sum((ShearXY.SX).*ShearXY.volume*1000))/volumeRVE)/10^6;
SS24=((sum((ShearXY.SY).*ShearXY.volume*1000))/volumeRVE)/10^6;
SS34=((sum((ShearXY.SZ).*ShearXY.volume*1000))/volumeRVE)/10^6;
SS44=((sum((ShearXY.SXY).*ShearXY.volume*1000))/volumeRVE)/10^6;
SS54=((sum((ShearXY.SYZ).*ShearXY.volume*1000))/volumeRVE)/10^6;
SS64=((sum((ShearXY.SXZ).*ShearXY.volume*1000))/volumeRVE)/10^6;
SS15=((sum((ShearYZ.SX).*ShearYZ.volume*1000))/volumeRVE)/10^6;
SS25=((sum((ShearYZ.SY).*ShearYZ.volume*1000))/volumeRVE)/10^6;
SS35=((sum((ShearYZ.SZ).*ShearYZ.volume*1000))/volumeRVE)/10^6;
SS45=((sum((ShearYZ.SXY).*ShearYZ.volume*1000))/volumeRVE)/10^6;
SS55=((sum((ShearYZ.SYZ).*ShearYZ.volume*1000))/volumeRVE)/10^6;
SS65=((sum((ShearYZ.SXZ).*ShearYZ.volume*1000))/volumeRVE)/10^6;
SS16=((sum((ShearXZ.SX).*ShearXZ.volume*1000))/volumeRVE)/10^6;
SS26=((sum((ShearXZ.SY).*ShearXZ.volume*1000))/volumeRVE)/10^6;
SS36=((sum((ShearXZ.SZ).*ShearXZ.volume*1000))/volumeRVE)/10^6;
SS46=((sum((ShearXZ.SXY).*ShearXZ.volume*1000))/volumeRVE)/10^6;
SS56=((sum((ShearXZ.SYZ).*ShearXZ.volume*1000))/volumeRVE)/10^6;
SS66=((sum((ShearXZ.SXZ).*ShearXZ.volume*1000))/volumeRVE)/10^6;
 
SS=abs([SS11 SS12 SS13 SS14 SS15 SS16;
        SS21 SS22 SS23 SS24 SS25 SS26;
        SS31 SS32 SS33 SS34 SS35 SS36;
        SS41 SS42 SS43 SS44 SS45 SS46;
        SS51 SS52 SS53 SS54 SS55 SS56;
        SS61 SS62 SS63 SS64 SS65 SS66;]);
 
EE=abs([EE11 EE12 EE13 EE14 EE15 EE16;
          EE21 EE22 EE23 EE24 EE25 EE26;
          EE31 EE32 EE33 EE34 EE35 EE36;
          EE41 EE42 EE43 EE44 EE45 EE46;
          EE51 EE52 EE53 EE54 EE55 EE56;
          EE61 EE62 EE63 EE64 EE65 EE66;]);
 
SS_theory= abs([SS(1,1) SS(1,2) SS(1,3) SS(1,5) SS(1,6) SS(1,4);
               SS(2,1) SS(2,2) SS(2,3) SS(2,5) SS(2,6) SS(2,4);
               SS(3,1) SS(3,2) SS(3,3) SS(3,5) SS(3,6) SS(3,4);
               SS(5,1) SS(5,2) SS(5,3) SS(5,5) SS(5,6) SS(5,4);
               SS(6,1) SS(6,2) SS(6,3) SS(6,5) SS(6,6) SS(6,4);
               SS(4,1) SS(4,2) SS(4,3) SS(4,5) SS(4,6) SS(4,4);]); 
 
EE_theory= abs([EE(1,1) EE(1,2) EE(1,3) EE(1,5) EE(1,6) EE(1,4);
               EE(2,1) EE(2,2) EE(2,3) EE(2,5) EE(2,6) EE(2,4);
               EE(3,1) EE(3,2) EE(3,3) EE(3,5) EE(3,6) EE(3,4);
               EE(5,1) EE(5,2) EE(5,3) EE(5,5) EE(5,6) EE(5,4);
               EE(6,1) EE(6,2) EE(6,3) EE(6,5) EE(6,6) EE(6,4);
               EE(4,1) EE(4,2) EE(4,3) EE(4,5) EE(4,6) EE(4,4);]);
 
%% Implementation of homogenization theory
TissueStiffnessTensor= SS*inv(EE);
ApparentStrainTensor = eye(6)*0.01;
ApparentStrainTensor(1,1)=-ApparentStrainTensor(1,1);
ApparentStrainTensor(2,2)=-ApparentStrainTensor(2,2);
ApparentStrainTensor(3,3)=-ApparentStrainTensor(3,3);
LocalStructureTensor=EE*inv(ApparentStrainTensor);
stiffness=(TissueStiffnessTensor*LocalStructureTensor);
TissueStiffnessTensor_theory= SS_theory*inv(EE_theory);
ApparentStrainTensor_theory = eye(6)*0.01;
ApparentStrainTensor_theory(1,1)=-ApparentStrainTensor_theory(1,1);
ApparentStrainTensor_theory(2,2)=-ApparentStrainTensor_theory(2,2);
ApparentStrainTensor_theory(3,3)=-ApparentStrainTensor_theory(3,3);
LocalStructureTensor_theory=EE_theory*inv(ApparentStrainTensor_theory);
Stiffness_theory=(TissueStiffnessTensor_theory*LocalStructureTensor_theory);
Theory_Stiffness(:,:,i-2)=Stiffness_theory;
Theory_Stress(:,:,i-2)=SS_theory;
Theory_Strain(:,:,i-2)=EE_theory;
 
uFE_Stiffness_AnsysInputContinuum(:,:,i-2)=abs([stiffness(1,1) stiffness(2,1) stiffness(3,1) stiffness(4,1)  stiffness(5,1) stiffness(6,1); stiffness(2,1) stiffness(2,2) stiffness(3,2) stiffness(4,2) stiffness(5,2) stiffness(6,2);                                              stiffness(3,1) stiffness(3,2) stiffness(3,3) stiffness(4,3) stiffness(5,3) stiffness(6,3);  stiffness(4,1) stiffness(4,2) stiffness(4,3) stiffness(4,4) stiffness(5,4) stiffness(6,4); stiffness(5,1) stiffness(5,2) stiffness(5,3) stiffness(5,4) stiffness(5,5) stiffness(6,5); stiffness(6,1) stiffness(6,2) stiffness(6,3) stiffness(6,4) stiffness(6,5), stiffness(6,6);]);
uFE_Stiffness_AnsysCorrection(:,:,i-2)=abs([stiffness(1,1) stiffness(2,1) stiffness(3,1) stiffness(4,1)/2 stiffness(5,1)/2 stiffness(6,1)/2; stiffness(2,1) stiffness(2,2) stiffness(3,2) stiffness(4,2)/2 stiffness(5,2)/2 stiffness(6,2)/2; stiffness(3,1) stiffness(3,2) stiffness(3,3) stiffness(4,3)/2 stiffness(5,3)/2 stiffness(6,3)/2;                                          stiffness(4,1)/2 stiffness(4,2)/2 stiffness(4,3)/2 stiffness(4,4) stiffness(5,4)/2 stiffness(6,4)/2;                                         stiffness(5,1)/2 stiffness(5,2)/2 stiffness(5,3)/2 stiffness(5,4)/2 stiffness(5,5) stiffness(6,5)/2;                                         stiffness(6,1)/2 stiffness(6,2)/2 stiffness(6,3)/2 stiffness(6,4)/2 stiffness(6,5)/2 stiffness(6,6);]);
 
uFE_Stiffness_TheoryInputContinuum(:,:,i-2)=abs([Stiffness_theory(1,1) Stiffness_theory(2,1) Stiffness_theory(3,1) Stiffness_theory(4,1) Stiffness_theory(5,1) Stiffness_theory(6,1);
Stiffness_theory(2,1) Stiffness_theory(2,2) Stiffness_theory(3,2) Stiffness_theory(4,2) Stiffness_theory(5,2) Stiffness_theory(6,2);Stiffness_theory(3,1) Stiffness_theory(3,2) Stiffness_theory(3,3) Stiffness_theory(4,3) Stiffness_theory(5,3) Stiffness_theory(6,3); Stiffness_theory(4,1) Stiffness_theory(4,2) Stiffness_theory(4,3) Stiffness_theory(4,4) Stiffness_theory(5,4) Stiffness_theory(6,4); Stiffness_theory(5,1) Stiffness_theory(5,2) Stiffness_theory(5,3) Stiffness_theory(5,4) Stiffness_theory(5,5) Stiffness_theory(6,5); Stiffness_theory(6,1) Stiffness_theory(6,2) Stiffness_theory(6,3) Stiffness_theory(6,4) Stiffness_theory(6,5) Stiffness_theory(6,6);]);

uFE_Stiffness_TheoryCorrection(:,:,i-2)=abs([Stiffness_theory(1,1) Stiffness_theory(2,1) Stiffness_theory(3,1) Stiffness_theory(4,1)/2 Stiffness_theory(5,1)/2 Stiffness_theory(6,1)/2; Stiffness_theory(2,1) Stiffness_theory(2,2) Stiffness_theory(3,2) Stiffness_theory(4,2)/2 Stiffness_theory(5,2)/2 Stiffness_theory(6,2)/2; Stiffness_theory(3,1) Stiffness_theory(3,2) Stiffness_theory(3,3) Stiffness_theory(4,3)/2 Stiffness_theory(5,3)/2 Stiffness_theory(6,3)/2;   Stiffness_theory(4,1)/2 Stiffness_theory(4,2)/2 Stiffness_theory(4,3)/2 Stiffness_theory(4,4) Stiffness_theory(5,4)/2 Stiffness_theory(6,4)/2; Stiffness_theory(5,1)/2 Stiffness_theory(5,2)/2 Stiffness_theory(5,3)/2 Stiffness_theory(5,4)/2 Stiffness_theory(5,5) Stiffness_theory(6,5)/2; Stiffness_theory(6,1)/2 Stiffness_theory(6,2)/2 Stiffness_theory(6,3)/2 Stiffness_theory(6,4)/2 Stiffness_theory(6,5)/2 Stiffness_theory(6,6);]);
 end 
end 
 
[a,b,c]=size(uFE_Stiffness_TheoryInputContinuum);
 
for i=1:c
    try
  fprintf('\n uFE PostProcessing: Cubes#%d --> %s  loading ...  %s\n\n',i,uFEM.name{i},datestr(now,'mmmm dd, yyyy HH:MM:SS'));
 
  uFE_Stiffness_TIContinuum=uFE_Stiffness_TheoryInputContinuum(:,:,i);
  uFE_Stiffness_TCorrection=uFE_Stiffness_TheoryCorrection(:,:,i);
 
    if (uFE_Stiffness_TIContinuum(1,1)>=uFE_Stiffness_TIContinuum(2,2) && uFE_Stiffness_TIContinuum(2,2)>=uFE_Stiffness_TIContinuum(3,3))
        l1=1; l2=0; l3=0; m1=0; m2=1; m3=0; n1=0; n2=0; n3=1; 
    elseif (uFE_Stiffness_TIContinuum(1,1)>=uFE_Stiffness_TIContinuum(3,3) && uFE_Stiffness_TIContinuum(3,3)>=uFE_Stiffness_TIContinuum(2,2))
        l1=1; l2=0; l3=0; m1=0; m2=0; m3=1; n1=0; n2=1; n3=0;    
    elseif (uFE_Stiffness_TIContinuum(2,2)>=uFE_Stiffness_TIContinuum(1,1) && uFE_Stiffness_TIContinuum(1,1)>=uFE_Stiffness_TIContinuum(3,3))
        l1=0; l2=1; l3=0; m1=1; m2=0; m3=0; n1=0; n2=0; n3=1;        
    elseif (uFE_Stiffness_TIContinuum(2,2)>=uFE_Stiffness_TIContinuum(3,3) && uFE_Stiffness_TIContinuum(3,3)>=uFE_Stiffness_TIContinuum(1,1))
        l1=0; l2=0; l3=1; m1=1; m2=0; m3=0; n1=0; n2=1; n3=0;   
    elseif (uFE_Stiffness_TIContinuum(3,3)>=uFE_Stiffness_TIContinuum(1,1) && uFE_Stiffness_TIContinuum(1,1)>=uFE_Stiffness_TIContinuum(2,2))
        l1=0; l2=1; l3=0; m1=0; m2=0; m3=1; n1=1; n2=0; n3=0;      
    elseif (uFE_Stiffness_TIContinuum(3,3)>=uFE_Stiffness_TIContinuum(2,2) && uFE_Stiffness_TIContinuum(2,2)>=uFE_Stiffness_TIContinuum(1,1))
        l1=0; l2=0; l3=1; m1=0; m2=1; m3=0; n1=1; n2=0; n3=0;  
    end 
    
   Evet1= [l1,l2,l3]; Evet2= [m1,m2,m3];  Evet3= [n1,n2,n3];
   
   rotm = [Evet1; Evet2; Evet3]; eulZYX = rotm2eul(rotm)
   degX_magnitude(i)=radtodeg(eulZYX(3)); degY_magnitude(i)=radtodeg(eulZYX(2));    degZ_magnitude(i)=radtodeg(eulZYX(1)); 
   
    Q(:,:,i) =  [l1^2        m1^2        n1^2         2*m1*n1              2*n1*l1             2*l1*m1;
                 l2^2        m2^2        n2^2         2*m2*n2              2*n2*l2             2*l2*m2;
                 l3^2        m3^2        n3^2         2*m3*n3              2*n3*l3             2*l3*m3;
                 l2*l3       m2*m3       n2*n3        m2*n3+m3*n2          n2*l3+n3*l2         m2*l3+m3*l2;
                 l3*l1       m3*m1       n3*n1        m3*n1+m1*n3          n3*l1+n1*l3         m3*l1+m1*l3;
                 l1*l2       m1*m2       n1*n2        m1*n2+m2*n1          n1*l2+n2*l1         m1*l2+m2*l1] ;
 
    uFE_Stiffness_TIC_Q(:,:,i)=Q(:,:,i)*uFE_Stiffness_TIContinuum*Q(:,:,i)'; 
    uFE_Stiffness_TC_Q(:,:,i)=Q(:,:,i)*uFE_Stiffness_TCorrection*Q(:,:,i)'; 
        
    T=uFE_Stiffness_TIC_Q(:,:,i); 
    uFE_Stiffness_AIC_Q(:,:,i)=[T(1,1) T(1,2) T(1,3) T(1,6) T(1,4) T(1,5);
                            T(2,1) T(2,2) T(2,3) T(2,6) T(2,4) T(2,5);
                            T(3,1) T(3,2) T(3,3) T(3,6) T(3,4) T(3,5);
                            T(6,1) T(6,2) T(6,3) T(6,6) T(6,4) T(6,5);
                            T(4,1) T(4,2) T(4,3) T(4,6) T(4,4) T(4,5);
                            T(5,1) T(5,2) T(5,3) T(5,6) T(5,4) T(5,5);];
 
    TT=uFE_Stiffness_TC_Q(:,:,i);
    uFE_Stiffness_AC_Q(:,:,i)=[TT(1,1) TT(1,2) TT(1,3) TT(1,6) TT(1,4) TT(1,5);
                           TT(2,1) TT(2,2) TT(2,3) TT(2,6) TT(2,4) TT(2,5);
                           TT(3,1) TT(3,2) TT(3,3) TT(3,6) TT(3,4) TT(3,5);
                           TT(6,1) TT(6,2) TT(6,3) TT(6,6) TT(6,4) TT(6,5);
                           TT(4,1) TT(4,2) TT(4,3) TT(4,6) TT(4,4) TT(4,5);
                           TT(5,1) TT(5,2) TT(5,3) TT(5,6) TT(5,4) TT(5,5);];
end
end
 
for i=1:c
    try
  fprintf('\n uFE Anisotropy Determination: Cubes#%d --> %s  loading ...  %s\n\n',i,uFEM.name{i},datestr(now,'mmmm dd, yyyy HH:MM:SS'));
Stiffness_ani=uFE_Stiffness_TIC_Q(:,:,i); 
cd('D:\June2021\QUT_February2021')             
opt = optimoptions('fmincon', 'StepTolerance', 1e-50, 'Display', 'iter', 'FiniteDifferenceStepSize', 1e-9);
fun=@(x)(t(Stiffness_ani))
       nonlcon = @CowinObjFun;
    x0 = [1 0 0 0 1 0 0 0 1];  A = []; b = []; Aeq = []; beq = []; lb = []; ub = [];
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,opt);
 
l1=x(1);l2=x(2);l3=x(3);m1=x(4);m2=x(5);m3=x(6);n1=x(7);n2=x(8);n3=x(9);
v1=(x(1))*(x(4))+(x(2))*(x(5))+(x(3))*(x(6));
v2=(x(1))*(x(7))+(x(2))*(x(8))+(x(3))*(x(9));
v3=(x(4))*(x(7))+(x(5))*(x(8))+(x(6))*(x(9));
v4=rad2deg(atan2(norm(cross([x(1),x(2),x(3)],[x(4),x(5),x(6)])), dot([x(1),x(2),x(3)],[x(4),x(5),x(6)])));
v5=rad2deg(atan2(norm(cross([x(1),x(2),x(3)],[x(7),x(8),x(9)])), dot([x(1),x(2),x(3)],[x(7),x(8),x(9)])));
v6=rad2deg(atan2(norm(cross([x(7),x(8),x(9)],[x(4),x(5),x(6)])), dot([x(7),x(8),x(9)],[x(4),x(5),x(6)])));
v7= norm([x(1),x(2),x(3)]);
v8= norm([x(4),x(5),x(6)]);
v9= norm([x(7),x(8),x(9)]);
 
Q(:,:,i) = [l1*l1      m1*m1      n1*n1       2*m1*n1          2*n1*l1         2*l1*m1;
           l2*l2      m2*m2      n2*n2       2*m2*n2          2*n2*l2         2*l2*m2;
           l3*l3      m3*m3      n3*n3       2*m3*n3          2*n3*l3         2*l3*m3;
           l2*l3      m2*m3      n2*n3       m2*n3+m3*n2      n2*l3+n3*l2     m2*l3+m3*l2;
           l3*l1      m3*m1      n3*n1       m3*n1+m1*n3      n3*l1+n1*l3     m3*l1+m1*l3;
           l1*l2      m1*m2      n1*n2       m1*n2+m2*n1      n1*l2+n2*l1     m1*l2+m2*l1] ;
  
   uFE_Stiffness_TheoryOrth_Q(:,:,i)=Q(:,:,i)*Stiffness_ani*Q(:,:,i)';
   T=uFE_Stiffness_TheoryOrth_Q(:,:,i); 
 
   uFE_Stiffness_AnsysOrth_Q(:,:,i)=abs([T(1,1) T(1,2) T(1,3) T(1,6) T(1,4) T(1,5);
                                    T(2,1) T(2,2) T(2,3) T(2,6) T(2,4) T(2,5);
                                    T(3,1) T(3,2) T(3,3) T(3,6) T(3,4) T(3,5);
                                    T(6,1) T(6,2) T(6,3) T(6,6) T(6,4) T(6,5);
                                    T(4,1) T(4,2) T(4,3) T(4,6) T(4,4) T(4,5);
                                    T(5,1) T(5,2) T(5,3) T(5,6) T(5,4) T(5,5);]);
                           
   Evet1= [l1,l2,l3];   Evet2= [m1,m2,m3];   Evet3= [n1,n2,n3];  rotm = [Evet1; Evet2; Evet3];
   eulZYX = rotm2eul(rotm);
       
   uFEM_Vect1x(i,1)=Evet1(3);  uFEM_Vect1y(i,1)=-Evet1(1);  uFEM_Vect1z(i,1)=-Evet1(2);
   uFEM_Vect2x(i,1)=Evet2(3);  uFEM_Vect2y(i,1)=-Evet2(1);  uFEM_Vect2z(i,1)=-Evet2(2);
   uFEM_Vect3x(i,1)=Evet3(3);  uFEM_Vect3y(i,1)=-Evet3(1);  uFEM_Vect3z(i,1)=-Evet3(2);  
    
   degX_orth(i,1)=radtodeg(eulZYX(3))+degX_magnitude(i);
   degY_orth(i,1)=radtodeg(eulZYX(2))+degY_magnitude(i);
   degZ_orth(i,1)=radtodeg(eulZYX(1))+degZ_magnitude(i);
   
   Sani_AniOrth=uFE_Stiffness_AnsysOrth_Q(:,:,i)*10^6;   Sort_Orh=Sani_AniOrth;
   Sort_Orh(4,1)=0;   Sort_Orh(5,1)=0;   Sort_Orh(6,1)=0;   Sort_Orh(4,2)=0;   Sort_Orh(5,2)=0;   Sort_Orh(6,2)=0;    Sort_Orh(4,3)=0;   Sort_Orh(5,3)=0;   Sort_Orh(6,3)=0;   Sort_Orh(1,4)=0;   Sort_Orh(1,5)=0;   Sort_Orh(1,6)=0;     Sort_Orh(2,4)=0;   Sort_Orh(2,5)=0;   Sort_Orh(2,6)=0;   Sort_Orh(3,4)=0;   Sort_Orh(3,5)=0;   Sort_Orh(3,6)=0;    Sort_Orh(5,4)=0;   Sort_Orh(6,4)=0;   Sort_Orh(4,5)=0;   Sort_Orh(6,5)=0;   Sort_Orh(4,5)=0;   Sort_Orh(4,6)=0;   Sort_Orh(5,4)=0;   Sort_Orh(5,6)=0;    
     
    end
end
 
    

