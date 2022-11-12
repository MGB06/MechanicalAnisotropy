clear all; clc; close all;
 
projectdir = 'C:\input_data';
dinfofolders = dir(fullfile(projectdir));
nfiles = length(dinfofolders); 
dinfofiles = dir(fullfile(projectdir));
 num=1:nfiles-2;
cd (projectdir)
fid = fopen('list_C122360.ctl','w'); 
fprintf(fid,'[Dataset list]\n');
 spacer='\'; fpath=strcat(dinfofolders(3).folder,spacer,dinfofolders(3).name,spacer);
cd (fpath) dinfofiles = dir(fullfile(fpath));
% the following string is specific for each CTAn 
infost='Info=0000000100000FFFF00000009F02603C8386AE3F00000000000000000400000000036';
pathdef=strcat(fpath,dinfofiles(3).name);
fprintf(fid,'Next=@0\n');
fprintf(fid,'[@0]\n');
fprintf(fid,'File=%s\n',pathdef);
fprintf(fid,'%s\n',infost);
 
for i=4:nfiles
spacer='\';
fpath=strcat(dinfofolders(i).folder,spacer,dinfofolders(i).name,spacer);
cd (fpath)
dinfofiles = dir(fullfile(fpath));
infost='Info=0000000100000FFFF00000009F02603C8386AE3F00000000000000000400000000036';
pathdef=strcat(fpath,dinfofiles(3).name);
fprintf(fid,'Next=@%d\n',num(i-3));
fprintf(fid,'[@%d] \n',num(i-3));
fprintf(fid,'File=%s\n',pathdef);
fprintf(fid,'%s\n',infost);
end
fclose(fid);
 
