%%%%% Master script to process LEADS PET data %%%%%
%%%%% Original code by Leo Iaccarino in 2019 %%%%%
%%%%% Updated by Daniel Schonhaut in 2024 %%%%%

% Define directory paths
dirs = containers.Map;
dirs('proj') = '/mnt/coredata/data/leads';
dirs('code') = fullfile(dirs('proj'), 'code');
dirs('data') = fullfile(dirs('proj'), 'data');
dirs('extraction') = fullfile(dirs('data'), 'extraction');
dirs('freesurfer') = fullfile(dirs('data'), 'freesurfer');
dirs('links') = fullfile(dirs('data'), 'links');
dirs('newdata') = fullfile(dirs('data'), 'newdata');
dirs('processed') = fullfile(dirs('data'), 'processed');

% Get directories in newdata
newdata_dirs = dir(fullfile(dirs('newdata'), '*'));
newdata_dirs = newdata_dirs([newdata_dirs.isdir]);
newdata_dirs = newdata_dirs(~ismember({newdata_dirs.name}, {'.', '..'}));
newdata_dirs = fullfile(dirs('newdata'), {newdata_dirs.name});

% Add paths to other scripts we need to access
addpath(genpath(dirs('code')));
addpath(genpath('/mnt/coredata/Projects/Resources/scripts/leotools'));

% Ask the user what processing to perform
msg = [
    '\nWhat do you want to do?' ,...
    '\n  [1] Process MRIs through FreeSurfer',...
    '\n  [2] Process PET through MRI-based pipeline',...
    '\n  [3] Process MRIs and PET',...
    '\n  [4] Extract PET SUVR means from FreeSurfer ROIs',...
    '\n  [5] Warp MRI and PET to MNI space',...
    '\n  [6] Create MRI and PET W-score maps',...
    '\n  [7] Run the full pipeline (Process MRI+PET, Extract ROIs, Warp to MNI, W-score)',...
    '\n  --> '
]
 user_action = input(msg);
    if user_action == 1
        run LEADS_convert_dicoms.m
        run LEADS_MRI_Processing.m

        % let's clean the corresponding newdata folder SAFELY

        cd(path_newmris);
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);

            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new MRI folder
            clear tobedeleted

        cd(path_master);
        clear;

    elseif user_action == 2
        run LEADS_convert_dicoms.m
        run LEADS_PET_Processing.m

        % let's clean the corresponding newdata folder

        cd(path_newfbbs);
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);

            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FBB folder
            clear tobedeleted

        cd(path_newftps);
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);

            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FTP folder
            clear tobedeleted

        cd(path_newfdgs);
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);

            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FDG folder
            clear tobedeleted

        cd(path_master);
        clear;

    elseif user_action == 3
        run LEADS_convert_dicoms.m
        run LEADS_MRI_Processing.m
        run LEADS_PET_Processing.m

        cd(path_newmris);
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);

            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new MRI folder
            clear tobedeleted

        cd(path_newfbbs);
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);

            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FBB folder
            clear tobedeleted

        cd(path_newftps);
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);

            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FTP folder
            clear tobedeleted

        cd(path_newfdgs);
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);

            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FDG folder
            clear tobedeleted

        cd(path_master);
        clear;

    elseif user_action == 4

        run LEADS_PET_Quantification.m
        clear;

    elseif user_action == 5

        run LEADS_PETMRI_MNI.m
        clear;

    elseif user_action == 6

        run LEADS_Wmapping.m
        clear;

    elseif user_action == 7
        run LEADS_convert_dicoms.m
        run LEADS_MRI_Processing.m
         run LEADS_PET_Processing.m

        % let's clean the corresponding newdata folder

        cd(path_newfbbs);
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);

            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FBB folder
            clear tobedeleted

        cd(path_newftps);
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FTP folder
            clear tobedeleted

        cd(path_newfdgs);
        tobedeleted=dir('LDS*');
        tobedeleted=transpose(struct2cell(tobedeleted));
        tobedeleted=tobedeleted(:,1);
            for i=1:size(tobedeleted,1)
            rmdir(tobedeleted{i},'s')
            end % end for loop for each new FDG folder
            clear tobedeleted

        cd(path_master);

        run LEADS_PET_Quantification.m
        clear;

        run LEADS_PETMRI_MNI.m
        clear;

        run LEADS_Wmapping.m
        clear;

    end


    fprintf(1,'Finished, thanks!\n');
    clear
