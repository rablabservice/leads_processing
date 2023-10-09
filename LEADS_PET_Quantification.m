%%% Script to run quantification for LEADS processed data %%%
%%% Sep 2021: Code Cleanup %%%

path_extraction='/mnt/coredata/Projects/LEADS/data_f7p1/extraction/'; % define extraction folder

cd (path_processed)

%%% Create empty arrays to fill with images we already have processed and
%%% ready

tdb_fbb={};
tdb_ftp={};   
tdb_fdg={};

%%% Scan through the images to populate the arrays

% 1. FBB
fprintf(1,'Gathering info on available FBBs.....\n');

avfbbs=dir('*/*/FBB*');
avfbbs=strcat({avfbbs.folder}','/',{avfbbs.name}');

for i=1:size(avfbbs,1)
    
    srcfbbsuvr=dir(strcat(avfbbs{i,1},'/LDS*suvr_cbl.nii'));
    srcfbbaparc=dir(strcat(avfbbs{i,1},'/LDS*raparc+aseg.nii'));

    if size(srcfbbsuvr,1)==1 && size(srcfbbaparc,1)==1
        
    temp={strcat(srcfbbsuvr.folder,'/',srcfbbsuvr.name),strcat(srcfbbaparc.folder,'/',srcfbbaparc.name)};
    tdb_fbb=vertcat(tdb_fbb,temp);
        
    else
        
        fprintf(2,'Warning! In the FBB Folder %s I could not find both the SUVR and APARC files. That is not normal!\n', avfbbs{i,1});
    
    end % end if condition there are SUVR and APARC Files in the folder
    
    clear srcfbbsuvr srcfbbaparc temp
    
end % end for condition for each FBB folder

% 2. FTP
fprintf(1,'Gathering info on available FTPs.....\n');

avftps=dir('*/*/FTP*');
avftps=strcat({avftps.folder}','/',{avftps.name}');

for i=1:size(avftps,1)
    
    srcftpsuvr=dir(strcat(avftps{i,1},'/LDS*suvr_infcblg.nii'));
    srcftpaparc=dir(strcat(avftps{i,1},'/LDS*raparc+aseg.nii'));

    if size(srcftpsuvr,1)==1 && size(srcftpaparc,1)==1
        
    temp={strcat(srcftpsuvr.folder,'/',srcftpsuvr.name),strcat(srcftpaparc.folder,'/',srcftpaparc.name)};
    tdb_ftp=vertcat(tdb_ftp,temp);
        
    else
        
        fprintf(2,'Warning! In the FTP Folder %s I could not find both the SUVR and APARC files. That is not normal!\n', avftps{i,1});
    
    end % end if condition there are SUVR and APARC Files in the folder
    
    clear srcftpsuvr srcftpaparc temp
    
end % end for condition for each FTP folder

% 3. FDG
fprintf(1,'Gathering info on available FDGs.....\n');

avfdgs=dir('*/*/FDG*');
avfdgs=strcat({avfdgs.folder}','/',{avfdgs.name}');

for i=1:size(avfdgs,1)
    
    srcfdgsuvr=dir(strcat(avfdgs{i,1},'/LDS*suvr_pons.nii'));
    srcfdgaparc=dir(strcat(avfdgs{i,1},'/LDS*raparc+aseg.nii'));

    if size(srcfdgsuvr,1)==1 && size(srcfdgaparc,1)==1
        
    temp={strcat(srcfdgsuvr.folder,'/',srcfdgsuvr.name),strcat(srcfdgaparc.folder,'/',srcfdgaparc.name)};
    tdb_fdg=vertcat(tdb_fdg,temp);
        
    else
        
        fprintf(2,'Warning! In the FDG Folder %s I could not find both the SUVR and APARC files. That is not normal!\n', avfdgs{i,1});
    
    end % end if condition there are SUVR and APARC Files in the folder
    
    clear srcfdgsuvr srcfdgaparc temp
    
end % end for condition for each FDG folder

cd (path_processed);

%%%%%%%%%% Now we have the list of images we could run the freesurfer
%%%%%%%%%% extraction on.

%%%% Running extraction tool for PET data %%%%

if size(tdb_fbb,1)>0
 run LEADS_PET_FBB_ROIExtraction.m
else
    fprintf(2,'There are no processed FBB-PET images available to extract values from.\n');
end % end if condition there are FBB available

if size(tdb_ftp,1)>0
 run LEADS_PET_FTP_ROIExtraction.m
else
    fprintf(2,'There are no processed FTP-PET images available to extract values from.\n');
end % end if condition there are FTP available

if size(tdb_fdg,1)>0
 run LEADS_PET_FDG_ROIExtraction.m
else
    fprintf(2,'There are no processed FDG-PET images available to extract values from.\n');
end % end if condition there are FDG available
 
fprintf(1,'Finished all the extractions, have a good one!\n');
clear;
