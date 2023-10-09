%%%%% Master script to process PET data from LEADS %%%%%
%%%%% 04/08/2019 Leonardo.Iaccarino@ucsf.edu %%%%%
%%%%% Sep 2021: Code Cleanup %%%%%

%% Set useful paths 
addpath('/mnt/coredata/Projects/LEADS/script_f7p1/');
addpath('/mnt/coredata/Projects/LEADS/script_f7p1/service/');
addpath('/mnt/coredata/Projects/LEADS/script_f7p1/templates/');
addpath('/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/');
addpath(genpath('/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/'));
addpath(genpath('/mnt/coredata/Projects/Resources/scripts/leotools/'))

path_master='/mnt/coredata/Projects/LEADS/data_f7p1/';
path_newmris='/mnt/coredata/Projects/LEADS/data_f7p1/newdata/mri/';
path_newfbbs='/mnt/coredata/Projects/LEADS/data_f7p1/newdata/fbb/';
path_newftps='/mnt/coredata/Projects/LEADS/data_f7p1/newdata/ftp/';
path_newfdgs='/mnt/coredata/Projects/LEADS/data_f7p1/newdata/fdg/';
path_processed='/mnt/coredata/Projects/LEADS/data_f7p1/processed/';
path_freesurfer='/mnt/coredata/Projects/LEADS/data_f7p1/freesurfer_processing/';
path_mni='/mnt/coredata/Projects/LEADS/data_f7p1/mni/';
path_qccompwm='/shared/petcore/Projects/LEADS/data_f7p1/summary/mri/wm_fsf/';

path_dbs='/mnt/coredata/Projects/LEADS/data_f7p1/dbs/';

path_wscore_meta='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/meta/';

%% Ask for required action

 mmethod = input(['\n\n Choose what analysis to run today:' ,...
'\n     [0] Raw data handling: convert dicoms to nifti',...  
'\n     [1] MRI Update and Processing',...
'\n     [2] PET Update and Processing',...
'\n     [3] MRI Update and Processing + PET Update and Processing',...
'\n     [4] PET Quantification',...
'\n     [5] MRI and PET warping to MNI',...
'\n     [6] MRI and PET Wmap Estimation',...
'\n     [7] Full PET Pipeline (Proc,Quant,Warp,Wmap)',...
'\n     [8] Path Database generator',...
'\n     [9] Patient Lookup and Data Grabber',...
'\n     [10] Group Comparisons Generator',...
'\n     --> ']);

    if mmethod==0
        
        run LEADS_DicomtoNifti_Conversion.m
        
        cd (path_master);
        clear;
        
    elseif mmethod==1

        run LEADS_MRI_Processing.m

        % let's clean the corresponding newdata folder SAFELY

        cd (path_newmris); 
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1); 
        
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new MRI folder
            clear tobedeleted

        cd (path_master);
        clear;


    elseif mmethod==2

        run LEADS_PET_Processing.m

        % let's clean the corresponding newdata folder

        cd (path_newfbbs); 
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1); 
        
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FBB folder
            clear tobedeleted

        cd (path_newftps); 
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1); 
        
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FTP folder
            clear tobedeleted

        cd (path_newfdgs); 
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1); 
        
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FDG folder
            clear tobedeleted

        cd (path_master);
        clear;

    elseif mmethod==3

        run LEADS_MRI_Processing.m
        run LEADS_PET_Processing.m

        cd (path_newmris); 
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1); 
        
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new MRI folder
            clear tobedeleted

        cd (path_newfbbs); 
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1); 
        
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FBB folder
            clear tobedeleted

        cd (path_newftps); 
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1); 
        
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FTP folder
            clear tobedeleted

        cd (path_newfdgs); 
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1); 
        
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FDG folder
            clear tobedeleted

        cd (path_master);
        clear;

    elseif mmethod==4

        run LEADS_PET_Quantification.m
        clear;

    elseif mmethod==5

        run LEADS_PETMRI_MNI.m
        clear;

    elseif mmethod==6

        run LEADS_Wmapping.m
        clear;

    elseif mmethod==7

         run LEADS_PET_Processing.m

        % let's clean the corresponding newdata folder

        cd (path_newfbbs); 
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1); 
        
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FBB folder
            clear tobedeleted

        cd (path_newftps); 
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1); 
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FTP folder
            clear tobedeleted
        
        cd (path_newfdgs); 
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);  
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FDG folder
            clear tobedeleted

        cd (path_master);

        run LEADS_PET_Quantification.m
        clear;

        run LEADS_PETMRI_MNI.m
        clear;

        run LEADS_Wmapping.m
        clear;

    elseif mmethod==8

        run LEADS_DatabasePaths.m
        clear;

    elseif mmethod==9

        run LEADS_lookup.m
        clear;

    elseif mmethod==10

        run LEADS_GrpComp.m
        clear;

    end % end if condition required action


    fprintf(1,'Finished, thanks!\n');
    clear
