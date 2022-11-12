clc; clear; close all
 
projectdir_input = 'D:\input_folder';
cd(projectdir_input); dinfo = dir(fullfile(projectdir_input));
resultsdir='D:\output_folder';
mkdir(resultsdir);
 
%Iteration over all the selected cubes
for t=3:length(dinfo)
cd ([resultsdir,dinfo(t).name]);
BCs{1}='Sxy'; BCs{2}='Sxz'; BCs{3}='Syz'; 
BCs{4}='Cx'; BCs{5}='Cy'; BCs{6}='Cz';
folder_model_BCs{1}=[model.mCT_location_output,model.fname,'\',BCs{1}];
folder_model_BCs{2}=[model.mCT_location_output,model.fname,'\',BCs{2}];
folder_model_BCs{3}=[model.mCT_location_output,model.fname,'\',BCs{3}];
folder_model_BCs{4}=[model.mCT_location_output,model.fname,'\',BCs{4}];
folder_model_BCs{5}=[model.mCT_location_output,model.fname,'\',BCs{5}];
folder_model_BCs{6}=[model.mCT_location_output,model.fname,'\',BCs{6}];
name_model_BCs{1}=[model.fname,'_',BCs{1}]; name_model_BCs{2}=[model.fname,'_',BCs{2}];
name_model_BCs{3}=[model.fname,'_',BCs{3}]; name_model_BCs{4}=[model.fname,'_',BCs{4}];
name_model_BCs{5}=[model.fname,'_',BCs{5}]; name_model_BCs{6}=[model.fname,'_',BCs{6}];
 
inputBCs{1}= 'disp_FE_shearXY'; inputBCs{2}= 'disp_FE_shearXZ'; 
inputBCs{3}= 'disp_FE_shearYZ'; inputBCs{4}= 'disp_FE_Xcompression'; 
inputBCs{5}= 'disp_FE_Ycompression'; inputBCs{6}= 'disp_FE_Zcompression';
 
for h=1:length(BCs)  
fprintf('Writing APDL - Load case %d: %s\n',h,BCs{h});
mkdir(sprintf('%s',BCs{h}));
 
cd (folder_model_BCs{h});
fid = fopen([[name_model_BCs{h}],'.ans'],'w');
fprintf(fid,'/CWD,''%s''\n',model.basedir);
fprintf(fid,'/TITLE,%s\n',name_model_BCs{h});
fprintf(fid,'/FILNAME,%s,0 \n',name_model_BCs{h});
fprintf(fid,'/PREP7\n');
fprintf(fid,'/INPUT,%s,nodes,\n',model.fname);
fprintf(fid,'/INPUT,%s,mlis,\n',model.fname);
fprintf(fid,'/INPUT,%s,elem,\n',model.fname);
fprintf(fid,'/INPUT,%s,txt,\n',inputBCs{h});
fprintf(fid,'/CWD,''%s''\n',folder_model_BCs{h});
fprintf(fid,'/SOLU\n');
fprintf(fid,'NROPT,full,,off,\n');
fprintf(fid,'solve \n');
fprintf(fid,'/POST1\n');
 
fprintf(fid,'!COLLECT NODAL DISPLACEMENTS\n');
fprintf(fid,'*GET,NNUMMAX,NODE, ,COUNT\n');
fprintf(fid,'*DIM,DX,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,DY,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,DZ,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*VGET,DX,NODE, ,U,X, , ,2\n');
fprintf(fid,'*VGET,DY,NODE, ,U,Y, , ,2\n');
fprintf(fid,'*VGET,DZ,NODE, ,U,Z, , ,2\n');
fprintf(fid,'*CFOPEN,%s_nodal_displacement.txt\n',inputBCs{h});
fprintf(fid,'*VWRITE,DX(1,1),DY(1,1),DZ(1,1)\n');
fprintf(fid,'(E15.7,E15.7,E15.7)\n');
fprintf(fid,'*CFCLOSE\n');
fprintf(fid,'*DEL,DX, ,\n');
fprintf(fid,'*DEL,DY, ,\n');
fprintf(fid,'*DEL,DZ, ,\n');
 
fprintf(fid,'!COLLECT NODAL REACTION FORCE\n');
fprintf(fid,'*GET,NNUMMAX,NODE, ,COUNT\n');
fprintf(fid,'*DIM,RFX,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,RFY,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,RFZ,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*VGET,RFX,NODE, ,RF,FX, , ,2\n');
fprintf(fid,'*VGET,RFY,NODE, ,RF,FY, , ,2\n');
fprintf(fid,'*VGET,RFZ,NODE, ,RF,FZ, , ,2\n');
fprintf(fid,'*CFOPEN,%s_nodal_reactions.txt\n',inputBCs{h});
fprintf(fid,'*VWRITE,RFX(1,1),RFY(1,1),RFZ(1,1)\n');
fprintf(fid,'(E15.7,E15.7,E15.7)\n');
fprintf(fid,'*CFCLOSE\n');
fprintf(fid,'*DEL,RFX, ,\n');
fprintf(fid,'*DEL,RFY, ,\n');
fprintf(fid,'*DEL,RFZ, ,\n');
 
fprintf(fid,'/post1\n');
fprintf(fid,'SET,FIRST\n');
fprintf(fid,'outres,all,all\n');
fprintf(fid,'mxpand,,,,yes,,yes\n');
fprintf(fid,'/post1\n');
fprintf(fid,'/output,%s_nodal_forces,txt\n',inputBCs{h});
fprintf(fid,'nforce,all\n');
fprintf(fid,'/output\n');
 
fprintf(fid,'FILE,''%s'',''rst'',''.'' \n',name_model_BCs{h});
fprintf(fid,'SET,LAST\n');
fprintf(fid,'*GET,NNUMMAX,NODE, ,COUNT\n');
fprintf(fid,'*DIM,SX,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,SY,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,SZ,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,SXY,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,SXZ,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,SYZ,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*VGET,SX,NODE, ,S,X, , ,2\n');
fprintf(fid,'*VGET,SY,NODE, ,S,Y, , ,2\n');
fprintf(fid,'*VGET,SZ,NODE, ,S,Z, , ,2\n');
fprintf(fid,'*VGET,SXY,NODE, ,S,XY, , ,2\n');
fprintf(fid,'*VGET,SXZ,NODE, ,S,XZ, , ,2\n');
fprintf(fid,'*VGET,SYZ,NODE, ,S,YZ, , ,2\n');
fprintf(fid,'*CFOPEN,%s_nodal_stresses.txt\n',inputBCs{h});
fprintf(fid,'*VWRITE,SX(1,1),SY(1,1),SZ(1,1),SXY(1,1),SXZ(1,1),SYZ(1,1)\n');
fprintf(fid,'(E15.7,E15.7,E15.7,E15.7,E15.7,E15.7)\n');
fprintf(fid,'*CFCLOSE\n');
fprintf(fid,'*DEL,SX, ,\n');
fprintf(fid,'*DEL,SY, ,\n');
fprintf(fid,'*DEL,SZ, ,\n');
fprintf(fid,'*DEL,SXY, ,\n');
fprintf(fid,'*DEL,SXZ, ,\n');
fprintf(fid,'*DEL,SYZ, ,\n');
 
fprintf(fid,'FILE,''%s'',''rst'',''.'' \n',name_model_BCs{h});
fprintf(fid,'SET,LAST\n');
fprintf(fid,'*GET,NNUMMAX,NODE, ,COUNT\n');
fprintf(fid,'*DIM,EX,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,EY,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,EZ,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,EXY,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,EXZ,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*DIM,EYZ,ARRAY,NNUMMAX,1\n');
fprintf(fid,'*VGET,EX,NODE, ,EPTO,X, , ,2\n');
fprintf(fid,'*VGET,EY,NODE, ,EPTO,Y, , ,2\n');
fprintf(fid,'*VGET,EZ,NODE, ,EPTO,Z, , ,2\n');
fprintf(fid,'*VGET,EXY,NODE, ,EPTO,XY, , ,2\n');
fprintf(fid,'*VGET,EXZ,NODE, ,EPTO,XZ, , ,2\n');
fprintf(fid,'*VGET,EYZ,NODE, ,EPTO,YZ, , ,2\n');
fprintf(fid,'*CFOPEN,%s_nodal_strains.txt\n',inputBCs{h});
fprintf(fid,'*VWRITE,EX(1,1),EY(1,1),EZ(1,1),EXY(1,1),EXZ(1,1),EYZ(1,1)\n');
fprintf(fid,'(E15.7,E15.7,E15.7,E15.7,E15.7,E15.7)\n');
fprintf(fid,'*CFCLOSE\n');
fprintf(fid,'*DEL,EX, ,\n');
fprintf(fid,'*DEL,EY, ,\n');
fprintf(fid,'*DEL,EZ, ,\n');
fprintf(fid,'*DEL,EXY, ,\n');
fprintf(fid,'*DEL,EXZ, ,\n');
fprintf(fid,'*DEL,EYZ, ,\n');
 
fprintf(fid,'! COLLECT ELEMENT RESULTS\n');
fprintf(fid,'*GET,ENUMMAX,ELEM, ,COUNT\n');
fprintf(fid,'*DIM,ERES,ARRAY,ENUMMAX,6\n');
fprintf(fid,'ETABLE,SX,S,x\n');
fprintf(fid,'ETABLE,SY,S,Y\n');
fprintf(fid,'ETABLE,SZ,S,Z\n');
fprintf(fid,'ETABLE,SXY,S,XY\n');
fprintf(fid,'ETABLE,SYZ,S,YZ\n');
fprintf(fid,'ETABLE,SXZ,S,XZ\n');
fprintf(fid,'*VGET,ERES(1,1),ELEM,1,ETAB,SX\n');
fprintf(fid,'*VGET,ERES(1,2),ELEM,1,ETAB,SY\n');
fprintf(fid,'*VGET,ERES(1,3),ELEM,1,ETAB,SZ\n');
fprintf(fid,'*VGET,ERES(1,4),ELEM,1,ETAB,SXY\n');
fprintf(fid,'*VGET,ERES(1,5),ELEM,1,ETAB,SYZ\n');
fprintf(fid,'*VGET,ERES(1,6),ELEM,1,ETAB,SXZ\n');
fprintf(fid,'*CFOPEN,%s_tensor_stress.txt\n',inputBCs{h});
fprintf(fid,'*VWRITE,ERES(1,1),ERES(1,2),ERES(1,3),ERES(1,4),ERES(1,5),ERES(1,6)\n');
fprintf(fid,'(E15.7,E15.7,E15.7,E15.7,E15.7,E15.7)\n');
fprintf(fid,'*CFCLOSE\n');
fprintf(fid,'*DEL,SX, ,\n');
fprintf(fid,'*DEL,SY, ,\n');
fprintf(fid,'*DEL,SZ, ,\n');
fprintf(fid,'*DEL,SXY, ,\n');
fprintf(fid,'*DEL,SYZ, ,\n');
fprintf(fid,'*DEL,SXZ, ,\n');
 
fprintf(fid,'! COLLECT ELEMENT RESULTS\n');
fprintf(fid,'*GET,ENUMMAX,ELEM, ,COUNT\n');
fprintf(fid,'*DIM,ERES,ARRAY,ENUMMAX,6\n');
fprintf(fid,'ETABLE,EX,EPTO,x\n');
fprintf(fid,'ETABLE,EY,EPTO,Y\n');
fprintf(fid,'ETABLE,EZ,EPTO,Z\n');
fprintf(fid,'ETABLE,EXY,EPTO,XY\n');
fprintf(fid,'ETABLE,EYZ,EPTO,YZ\n');
fprintf(fid,'ETABLE,EXZ,EPTO,XZ\n');
fprintf(fid,'*VGET,ERES(1,1),ELEM,1,ETAB,EX\n');
fprintf(fid,'*VGET,ERES(1,2),ELEM,1,ETAB,EY\n');
fprintf(fid,'*VGET,ERES(1,3),ELEM,1,ETAB,EZ\n');
fprintf(fid,'*VGET,ERES(1,4),ELEM,1,ETAB,EXY\n');
fprintf(fid,'*VGET,ERES(1,5),ELEM,1,ETAB,EYZ\n');
fprintf(fid,'*VGET,ERES(1,6),ELEM,1,ETAB,EXZ\n');
fprintf(fid,'*CFOPEN,%s_tensor_strain.txt\n',inputBCs{h});
fprintf(fid,'*VWRITE,ERES(1,1),ERES(1,2),ERES(1,3),ERES(1,4),ERES(1,5),ERES(1,6)\n');
fprintf(fid,'(E15.7,E15.7,E15.7,E15.7,E15.7,E15.7)\n');
fprintf(fid,'*CFCLOSE\n');
fprintf(fid,'*DEL,EX, ,\n');
fprintf(fid,'*DEL,EY, ,\n');
fprintf(fid,'*DEL,EZ, ,\n');
fprintf(fid,'*DEL,EXY, ,\n');
fprintf(fid,'*DEL,EYZ, ,\n');
fprintf(fid,'*DEL,EXZ, ,\n');
 
fprintf(fid,'! COLLECT ELEMENT RESULTS\n');
fprintf(fid,'*GET,ENUMMAX,ELEM, ,COUNT\n')
fprintf(fid,'*DIM,ERES,ARRAY,ENUMMAX,1\n');
fprintf(fid,'ETABLE,volume,VOLU\n');
fprintf(fid,'*VGET,ERES(1,1),ELEM,1,ETAB,volume\n');
fprintf(fid,'*CFOPEN,%s_volume.txt\n',inputBCs{h});
fprintf(fid,'*VWRITE,ERES(1,1)\n');
fprintf(fid,'(E15.7)\n');
fprintf(fid,'*CFCLOSE\n');
fprintf(fid,'*DEL,volume, ,\n');
 
fprintf(fid,'cdwrite,all,''%s'',cdb,\n',name_model_BCs{h});
 
fprintf(fid,'FINISH\n');
fclose(fid);
 
fprintf('\nSimulation loading case %d: %s in progress...\n',h,BCs{h});
fprintf('%s\n\n',datestr(now,'mmmm dd, yyyy HH:MM:SS'));
    str =['"',model.ANSYS_dir_exe,'" -p ansys',...
        ' -dir "',folder_model_BCs{h},'"',...
        ' -j "',name_model_BCs{h},'"',...
        ' -s read -l en-us -b ',...
        ' -i "',name_model_BCs{h},'.ans"',...
        ' -o "',name_model_BCs{h},'.out"',...
        ' -np 1'];
    system(str);
        
fprintf('\nSimulation loading case %d complete!\n',h);
fprintf('%s\n\n',datestr(now,'mmmm dd, yyyy HH:MM:SS'));
cd (model.basedir);
end 
end 

