close all; clear all;clc 
dir_saving='D:\output_models'
mkdir(dir_saving);
projectdir_input = 'D:\input';
dinfo = dir(fullfile(projectdir_input));
pixels=84;  
p=5/10000; %cube length size in mm 
 
%% Model Generation
 parfor t=3:length(dinfo)
     try
 
mkdir([dir_saving,dinfo(t).name])
cd([dir_saving,dinfo(t).name])
fprintf('Model generation %d - Cubes #%s\n',t-2,dinfo(t).name)
 
model_fname            = fullfile( dinfo(t).name);
model_basedir          = pwd;
model_results          = [model_basedir,'\',model_fname,'\','results'];
model_mCT_im_folder    = fullfile( dinfo(t).name);  
model_gray_thr         = [80 255];   
model_MFE_size         = [1, 1, 84; 1, 1, 84;  1, 1, 84];
model_MFE_softMatNum   = 5;  
model_check_thr_levels = [500 2500/2 2000];  
model_node_file        = [model_fname,'.nodes'];  
model_elem_file        = [model_fname,'.elem'];  
model_mat_file         = [model_fname,'.mlis']; 
model_MFE_Young_Modul  = [10 17000000000];     
model_ANSYS_dir_exe    = 'C:\Program Files\Ansys Inc\v201\ansys\bin\winx64\MAPDL.exe';
model_volume = loadImageFolder([model_mCT_location,model_mCT_im_folder],'BMP');
[mFE,mlist,inde,mater] = material_model(model_volume,...
    model_gray_thr(1),model_gray_thr(2),model_MFE_Young_Modul(1),model_MFE_Young_Modul(2));
    [elem,nodes,inde,indn] = gen_lattice_v2(mFE);
    nodes = double(nodes);
    nodes(:,2) = nodes(:,2) + repmat(model_MFE_size(2,1),size(nodes,1),1,1); % mesh
 
    e = zeros(size(elem));
    for i = 1:8
        [c,ia,ib] = intersect(indn,elem(:,i));  e(:,i) = ia;
    end
    eln = [1:length(elem)]';   m =  mFE(inde);     elem = cat(2,eln,e,double(m(:)));
 
    resolution=P/pixels;
    writeAnsysElementFileIsotropicIntensityMaterials(elem,model_elem_file);
    writeAnsysIsotropicMaterialFile(model_MFE_Young_Modul(2),model_mat_file);
    writeAnsysNodeFile([[1:length(nodes)]',nodes*resolution],model_node_file);  
 
    addpath(genpath(model_lattice_gen_path));
    file_nodes = [model_mCT_location_output,model_mCT_im_folder,'\',model_node_file];
    [testv3] = import_unodes(file_nodes);
 
    disp=-0.00005;     shear_disp=0.000025;
 
%% Generation of Loading condition file
vinc_x_elements=find(x==min(x));  
vinc_x_numbers=length(vinc_x_elements); vinc_x_coordinate_x=x(vinc_x_elements);       
vinc_x_coordinate_y=y(vinc_x_elements);  vinc_x_coordinate_z=z(vinc_x_elements);      
 
vinc_y_elements=find(y==min(y));   
vinc_y_numbers=length(vinc_y_elements); vinc_y_coordinate_x=x(vinc_y_elements);       
vinc_y_coordinate_y=y(vinc_y_elements);  vinc_y_coordinate_z=z(vinc_y_elements);       
 
vinc_z_elements=find(z==min(z));   
vinc_z_numbers=length(vinc_z_elements); vinc_z_coordinate_x=x(vinc_z_elements);       
vinc_z_coordinate_y=y(vinc_z_elements);  vinc_z_coordinate_z=z(vinc_z_elements);       
 
disp_x_elements=find(x==max(x));     
disp_x_numbers=length(disp_x_elements);  disp_x_coordinate_x=x(disp_x_elements);        
disp_x_coordinate_y=y(disp_x_elements);       
disp_x_coordinate_z=z(disp_x_elements);        
 
disp_y_elements=find(y==max(y));   
disp_y_numbers=length(disp_y_elements);  disp_y_coordinate_x=x(disp_y_elements);      
disp_y_coordinate_y=y(disp_y_elements);    disp_y_coordinate_z=z(disp_y_elements);         
 
disp_z_elements=find(z==max(z));    
disp_z_numbers=length(disp_z_elements);   disp_z_coordinate_x=x(disp_z_elements);     
disp_z_coordinate_y=y(disp_z_elements);    disp_z_coordinate_z=z(disp_z_elements);        
 
           %% -------------- Compression X -------------- 
    fid = fopen('disp_FE_Xcompression.txt','w');
    fprintf(fid,'!Defining constrains \n');
            for p=1:vinc_x_numbers
                    fprintf(fid,'D,%d, uy,  0\n',vinc_x_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',vinc_x_elements(p));
            end
            for p=1:disp_x_numbers
                    fprintf(fid,'D,%d, uy,  0\n',disp_x_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',disp_x_elements(p));
            end
            
            for p=1:vinc_y_numbers
                    fprintf(fid,'D,%d, ux,  0\n',vinc_y_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',vinc_y_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',vinc_y_elements(p));
            end    
            for p=1:disp_y_numbers
                    fprintf(fid,'D,%d, ux,  0\n',disp_y_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',disp_y_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',disp_y_elements(p));
            end
            for p=1:vinc_z_numbers
                    fprintf(fid,'D,%d, ux,  0\n',vinc_z_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',vinc_z_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',vinc_z_elements(p));
            end       
            for p=1:disp_z_numbers
                    fprintf(fid,'D,%d, ux,  0\n',disp_z_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',disp_z_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',disp_z_elements(p));
            end
    fprintf(fid,'\n!Defining displacements \n');
            for j=1:disp_x_numbers
                    fprintf(fid,'D,%d, ux,  %d\n',disp_x_elements(j),shear_disp);
            end
            for j=1:vinc_x_numbers
                    fprintf(fid,'D,%d, ux,  %d\n',vinc_x_elements(j),-shear_disp);
            end
    fclose(fid);
 
           %% -------------- Compression Y --------------  
    fid = fopen('disp_FE_Ycompression.txt','w');
    fprintf(fid,'!Defining constrains \n');
            for p=1:vinc_y_numbers
                    fprintf(fid,'D,%d, ux,  0\n',vinc_y_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',vinc_y_elements(p));
            end
            for p=1:disp_y_numbers
                    fprintf(fid,'D,%d, ux,  0\n',disp_y_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',disp_y_elements(p));
            end
            
            for p=1:vinc_x_numbers
                    fprintf(fid,'D,%d, ux,  0\n',vinc_x_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',vinc_x_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',vinc_x_elements(p));
            end    
            for p=1:disp_x_numbers
                    fprintf(fid,'D,%d, ux,  0\n',disp_x_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',disp_x_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',disp_x_elements(p));
            end
            for p=1:vinc_z_numbers
                    fprintf(fid,'D,%d, ux,  0\n',vinc_z_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',vinc_z_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',vinc_z_elements(p));
            end       
            for p=1:disp_z_numbers
                    fprintf(fid,'D,%d, ux,  0\n',disp_z_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',disp_z_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',disp_z_elements(p));
            end
    fprintf(fid,'\n!Defining displacements \n');
            for j=1:disp_y_numbers
                    fprintf(fid,'D,%d, uy,  %d\n',disp_y_elements(j),shear_disp);
            end
            for j=1:vinc_y_numbers
                    fprintf(fid,'D,%d, uy,  %d\n',vinc_y_elements(j),-shear_disp);
            end
    fclose(fid);
 
 
           %% -------------- Compression Z --------------  
    fid = fopen('disp_FE_Zcompression.txt','w');
    fprintf(fid,'!Defining constrains \n');
            for p=1:disp_z_numbers
                    fprintf(fid,'D,%d, ux,  0\n',disp_z_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',disp_z_elements(p));
            end
            for p=1:vinc_z_numbers
                    fprintf(fid,'D,%d, ux,  0\n',vinc_z_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',vinc_z_elements(p));
            end
            
            for p=1:vinc_x_numbers
                    fprintf(fid,'D,%d, ux,  0\n',vinc_x_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',vinc_x_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',vinc_x_elements(p));
            end    
            for p=1:disp_x_numbers
                    fprintf(fid,'D,%d, ux,  0\n',disp_x_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',disp_x_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',disp_x_elements(p));
            end
            for p=1:vinc_y_numbers
                    fprintf(fid,'D,%d, ux,  0\n',vinc_y_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',vinc_y_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',vinc_y_elements(p));
            end       
            for p=1:disp_y_numbers
                    fprintf(fid,'D,%d, ux,  0\n',disp_y_elements(p));
                    fprintf(fid,'D,%d, uy,  0\n',disp_y_elements(p));
                    fprintf(fid,'D,%d, uz,  0\n',disp_y_elements(p));
            end
            
    fprintf(fid,'\n!Defining displacements \n');
            for j=1:disp_z_numbers
                    fprintf(fid,'D,%d, uz,  %d\n',disp_z_elements(j),shear_disp);
            end
            for j=1:vinc_z_numbers
                    fprintf(fid,'D,%d, uz,  %d\n',vinc_z_elements(j),-shear_disp);
            end
    fclose(fid);
 
    
           %% -------------- Shear XY --------------  
    fid = fopen('disp_FE_shearXY.txt','w');
 
    fprintf(fid,'!Defining constrains \n');
        for p=1:disp_x_numbers
                fprintf(fid,'D,%d, ux,  0\n',disp_x_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',disp_x_elements(p));
        end
        for p=1:vinc_x_numbers
                fprintf(fid,'D,%d, ux,  0\n',vinc_x_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',vinc_x_elements(p));
        end
        for p=1:disp_y_numbers
                fprintf(fid,'D,%d, uy,  0\n',disp_y_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',disp_y_elements(p));
        end
        for p=1:vinc_y_numbers
                fprintf(fid,'D,%d, uy,  0\n',vinc_y_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',vinc_y_elements(p));
        end
        for p=1:disp_z_numbers
                fprintf(fid,'D,%d, ux,  0\n',disp_z_elements(p));
                fprintf(fid,'D,%d, uy,  0\n',disp_z_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',disp_z_elements(p));
        end
        for p=1:vinc_z_numbers
                fprintf(fid,'D,%d, ux,  0\n',vinc_z_elements(p));
                fprintf(fid,'D,%d, uy,  0\n',vinc_z_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',vinc_z_elements(p));
        end
 
    fprintf(fid,'\n!Defining displacements \n');
         for j=1:disp_x_numbers
                 fprintf(fid,'D,%d, uy,  %d\n',disp_x_elements(j),shear_disp);
         end
         for j=1:vinc_x_numbers
                fprintf(fid,'D,%d, uy,  %d\n',vinc_x_elements(j),-shear_disp);
        end
         for j=1:disp_y_numbers
                 fprintf(fid,'D,%d, ux,  %d\n',disp_y_elements(j),shear_disp);
         end
         for j=1:vinc_y_numbers
                fprintf(fid,'D,%d, ux,  %d\n',vinc_y_elements(j),-shear_disp);
         end
 fclose(fid);  
 
         
     %% -------------- Shear XZ -------------- 
    fid = fopen('disp_FE_shearXZ.txt','w');
 
    fprintf(fid,'!Defining constrains \n');
        for p=1:disp_x_numbers
                fprintf(fid,'D,%d, ux,  0\n',disp_x_elements(p));
                fprintf(fid,'D,%d, uy,  0\n',disp_x_elements(p));
        end
        for p=1:vinc_x_numbers
                fprintf(fid,'D,%d, ux,  0\n',vinc_x_elements(p));
                fprintf(fid,'D,%d, uy,  0\n',vinc_x_elements(p));
        end
        
        for p=1:disp_y_numbers
                fprintf(fid,'D,%d, ux,  0\n',disp_y_elements(p));
                fprintf(fid,'D,%d, uy,  0\n',disp_y_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',disp_y_elements(p));
        end
        for p=1:vinc_y_numbers
                fprintf(fid,'D,%d, ux,  0\n',vinc_y_elements(p));
                fprintf(fid,'D,%d, uy,  0\n',vinc_y_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',vinc_y_elements(p));
        end
 
        for p=1:disp_z_numbers
                fprintf(fid,'D,%d, uy,  0\n',disp_z_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',disp_z_elements(p));
        end         
        for p=1:vinc_z_numbers
                fprintf(fid,'D,%d, uy,  0\n',vinc_z_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',vinc_z_elements(p));
        end
        
        fprintf(fid,'\n!Defining displacements \n');
        for j=1:disp_z_numbers
                fprintf(fid,'D,%d, ux,  %d\n',disp_z_elements(j),shear_disp);
        end
        for j=1:vinc_z_numbers
                fprintf(fid,'D,%d, ux,  %d\n',vinc_z_elements(j),-shear_disp);
        end
        for j=1:disp_x_numbers
                 fprintf(fid,'D,%d, uz,  %d\n',disp_x_elements(j),shear_disp);
        end
        for j=1:vinc_x_numbers
                 fprintf(fid,'D,%d, uz,  %d\n',vinc_x_elements(j),-shear_disp);
        end
        fclose(fid);    
 
           %% -------------- Shear YZ -------------- =
    fid = fopen('disp_FE_shearYZ.txt','w');
 
    fprintf(fid,'!Defining constrains \n');
        for p=1:disp_z_numbers
                fprintf(fid,'D,%d, ux,  0\n',disp_z_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',disp_z_elements(p));
        end
        for p=1:vinc_z_numbers
                fprintf(fid,'D,%d, ux,  0\n',vinc_z_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',vinc_z_elements(p));
        end
        
        for p=1:disp_x_numbers
                fprintf(fid,'D,%d, ux,  0\n',disp_x_elements(p));
                fprintf(fid,'D,%d, uy,  0\n',disp_x_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',disp_x_elements(p));
        end
        for p=1:vinc_x_numbers
                fprintf(fid,'D,%d, ux,  0\n',vinc_x_elements(p));
                fprintf(fid,'D,%d, uy,  0\n',vinc_x_elements(p));
                fprintf(fid,'D,%d, uz,  0\n',vinc_x_elements(p));
        end
 
        for p=1:disp_y_numbers
                fprintf(fid,'D,%d, ux,  0\n',disp_y_elements(p));
                fprintf(fid,'D,%d, uy,  0\n',disp_y_elements(p));
        end
        for p=1:vinc_y_numbers
                fprintf(fid,'D,%d, ux,  0\n',vinc_y_elements(p));
                fprintf(fid,'D,%d, uy,  0\n',vinc_y_elements(p));
        end
        
    fprintf(fid,'\n!Defining displacements \n');
        for j=1:disp_y_numbers
                fprintf(fid,'D,%d, uz,  %d\n',disp_y_elements(j),shear_disp);
        end
        for j=1:vinc_y_numbers
                fprintf(fid,'D,%d, uz,  %d\n',vinc_y_elements(j),-shear_disp);
        end
        for j=1:disp_z_numbers
                 fprintf(fid,'D,%d, uy,  %d\n',disp_z_elements(j),shear_disp);
        end
        for j=1:vinc_z_numbers
                 fprintf(fid,'D,%d, uy,  %d\n',vinc_z_elements(j),-shear_disp);
        end
        fclose(fid);
 
     end
    
end 
