%% Script to automatically grab images for a list of patients
%% Created June 4th 2020 
%% Last update 08/25/2021 by Leonardo Iaccarino - Leonardo.Iaccarino@ucsf.edu
%%% Sep 2021: Code Cleanup

% Get directory where to store data
outdir=uigetdir(pwd,'Select output directory');

% Get label for folder in which to put the data, useful for log, and make
% dir
labfold=input(['\n\n Enter label for the datasearch (e.g. AD_only):' ,...
        '\n     --> '],'s');
labfold2=strcat(labfold,'_',datestr(now,'mm-dd-yyyy_HH-MM-SS'));
mkdir(char(strcat(outdir,'/',labfold2,'/'))); 

% Ask User what Timepoint they are interested in
            
tpchc=input(['\n\n Enter Timepoint(s) of interest (if multiple, space separated):' ,...
        '\nExample: 1 2', ...
        '\n     --> '],'s');
tpchc=str2num(tpchc);

% Ask what kind of images

imgtypes = input(['\n\n What Image Type are you looking for? (if multiple, space separated):' ,...
                '\nExample: 1 2', ...
                '\n     [1] Native',...
                '\n     [2] Warped',...
                '\n     [3] Affine Warped',...
                '\n     [4] Wmaps',...
                '\n     [5] Multislices',...
                '\n     [6] Wmap 3D renderings',...
                '\n     --> '],'s');
imgtypes=strsplit(imgtypes,' ');
            
% Finally ask for the list of images

ids=input(['\n\n Enter Vector with LEADS IDs (space separated):' ,...
        '\nExample: LDS0730024 LDS0730044 LDS0730101', ...
        '\n     --> '],'s');
ids = strsplit(ids,' ');
            

% Time to finally copy the data

for t=1:size(tpchc,2)
    
    % First create directories
    
    mkdir(char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/'))); 
    mkdir(char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fbb/'))); 
    mkdir(char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/ftp/'))); 
    mkdir(char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fdg/'))); 
    mkdir(char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/mri/'))); 
    
    % now query needed data and copy it
    
    for i=1:size(ids,2)
        
        if exist(strcat(path_processed,ids{i}),'dir')==7

            if exist(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t))),'dir')==7

                for cc=1:size(imgtypes,2)
                
                    imgtype=str2num(imgtypes{cc});
                    
                    if imgtype==1 % Looking for native space data

                        % MRI data must be there

                        srcnu=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/MRI*/',ids{i},'*nu.nii'));
                        srcnu=strcat(srcnu.folder,'/',srcnu.name);

                        srcaparc=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/MRI*/',ids{i},'*raparc+aseg.nii'));
                        srcaparc=strcat(srcaparc.folder,'/',srcaparc.name);

                        srcbs=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/MRI*/',ids{i},'*rbrainstemSsLabels_v12_VoxSpace.nii'));
                        srcbs=strcat(srcbs.folder,'/',srcbs.name);

                        copyfile(srcnu,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/mri/')))
                        copyfile(srcaparc,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/mri/')))
                        copyfile(srcbs,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/mri/')))

                        clear srcnu srcaparc srcbs


                        % FBB data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*')),1)==1 % existence of FBB folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*/',ids{i},'*suvr_cbl.nii'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fbb/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FBB does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FBB

                        % FTP data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*')),1)==1 % existence of FTP folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*/',ids{i},'*suvr_infcblg.nii'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/ftp/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FTP does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FTP

                        % FDG data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*')),1)==1 % existence of FDG folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*/',ids{i},'*suvr_pons.nii'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fdg/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FDG does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FDG

                    elseif imgtype==2 % Looking for warped data

                        % MRI data must be there

                        srcnu=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/MRI*/s8iso_mwc1',ids{i},'*nu.nii'));
                        srcnu=strcat(srcnu.folder,'/',srcnu.name);

                        copyfile(srcnu,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/mri/')))

                        clear srcnu 

                        % FBB data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*')),1)==1 % existence of FBB folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*/w',ids{i},'*suvr_cbl.nii'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fbb/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FBB does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FBB

                        % FTP data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*')),1)==1 % existence of FTP folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*/w',ids{i},'*suvr_infcblg.nii'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/ftp/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FTP does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FTP

                        % FDG data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*')),1)==1 % existence of FDG folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*/w',ids{i},'*suvr_pons.nii'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fdg/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FDG does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FDG

                    elseif imgtype==3 % Looking for Affine Warped data

                     % MRI data must be there

                        srcnu=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/MRI*/w_affine',ids{i},'*nu.nii'));
                        srcnu=strcat(srcnu.folder,'/',srcnu.name);

                        copyfile(srcnu,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/mri/')))

                        clear srcnu 

                        % FBB data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*')),1)==1 % existence of FBB folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*/w_affine',ids{i},'*suvr_cbl.nii'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fbb/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FBB does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FBB

                        % FTP data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*')),1)==1 % existence of FTP folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*/w_affine',ids{i},'*suvr_infcblg.nii'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/ftp/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FTP does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FTP

                        % FDG data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*')),1)==1 % existence of FDG folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*/w_affine',ids{i},'*suvr_pons.nii'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fdg/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FDG does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FDG

                    elseif imgtype==4

                         % MRI data must be there

                        srcnu=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/MRI*/Wmap_s8iso_mwc1',ids{i},'*nu.nii'));
                        srcnu=strcat(srcnu.folder,'/',srcnu.name);

                        copyfile(srcnu,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/mri/')))

                        clear srcnu 

                        % FBB data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*')),1)==1 % existence of FBB folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*/Wmap_w',ids{i},'*suvr_cbl.nii'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fbb/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FBB does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FBB

                        % FTP data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*')),1)==1 % existence of FTP folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*/Wmap_w',ids{i},'*suvr_infcblg.nii'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/ftp/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FTP does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FTP

                        % FDG data, check if it exists first. Wmaps not
                        % available yet for FDG

    %                     if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*')),1)==1 % existence of FDG folder 
    % 
    %                         srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*/Wmap_w',ids{i},'*suvr_pons.nii'));
    %                         srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);
    % 
    %                         copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fdg/')))
    % 
    %                         clear srcsuvr
    % 
    %                     else
    %                        fprintf(2,'Warning! FDG does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
    %                     end % end if condition existence of FDG


                    elseif imgtype==5

                     % MRI data must be there

                        srcnu=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/MRI*/Multiaxial_',ids{i},'*nu.jpg'));
                        srcnu=strcat(srcnu.folder,'/',srcnu.name);

                        copyfile(srcnu,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/mri/')))

                        clear srcnu 

                        % FBB data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*')),1)==1 % existence of FBB folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*/Multiaxial_',ids{i},'*suvr_cbl.jpg'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fbb/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FBB does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FBB

                        % FTP data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*')),1)==1 % existence of FTP folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*/Multiaxial_',ids{i},'*suvr_infcblg.jpg'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/ftp/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FTP does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FTP

                        % FDG data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*')),1)==1 % existence of FDG folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*/Multiaxial_',ids{i},'*suvr_pons.jpg'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fdg/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FDG does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FDG


                    elseif imgtype==6


                       % MRI data must be there

                        srcnu=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/MRI*/3DRend_Wmap_s8iso_mwc1',ids{i},'*_nu_GM.jpg'));
                        srcnu=strcat(srcnu.folder,'/',srcnu.name);

                        copyfile(srcnu,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/mri/')))

                        clear srcnu 

                        % FBB data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*')),1)==1 % existence of FBB folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FBB*/3DRend_Wmap_w',ids{i},'*suvr_cbl_GM.jpg'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fbb/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FBB does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FBB

                        % FTP data, check if it exists first

                        if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*')),1)==1 % existence of FTP folder 

                            srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FTP*/3DRend_Wmap_w',ids{i},'*suvr_infcblg_GM.jpg'));
                            srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);

                            copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/ftp/')))

                            clear srcsuvr

                        else
                           fprintf(2,'Warning! FTP does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
                        end % end if condition existence of FTP

                        % FDG data, check if it exists first. Wmap do not exist
                        % for FDG yet

    %                     if size(dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*')),1)==1 % existence of FDG folder 
    % 
    %                         srcsuvr=dir(strcat(path_processed,ids{i},'/Timepoint',num2str(tpchc(t)),'/FDG*/3DRend_Wmap_w',ids{i},'*suvr_pons_GM.jpg'));
    %                         srcsuvr=strcat(srcsuvr.folder,'/',srcsuvr.name);
    % 
    %                         copyfile(srcsuvr,char(strcat(outdir,'/',labfold2,'/Timepoint',num2str(tpchc(t)),'/fdg/')))
    % 
    %                         clear srcsuvr
    % 
    %                     else
    %                        fprintf(2,'Warning! FDG does not exist for ID %s at Timepoint %d in the LEADS repository. I am skipping it!\n',ids{i},tpchc(t));
    %                     end % end if condition existence of FDG  


                    end % end if condition image types requested
                
                end % for each imgtype requested
                
            else
                fprintf(2,'Warning! Timepoint %d does not exist for ID %s in the LEADS repository. I am skipping it!\n',tpchc(t),ids{i});
            end % end if condition the requested timepoint exists for the given subject


        else
            fprintf(2,'Warning! The ID %s does not exist in the LEADS repository. I am skipping it!\n',ids{i});
        end % end if condition the subject exists
    
    end % end for loop for each ID requested
    
end % end for loop for each requested ID

fprintf(1,'Just completed, have fun!\n');
clear;


