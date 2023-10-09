%% All the new FBBs have been now moved and renamed, so it is time to process them
%% This script still scouts the whole LEADS database to try and process any scan in need. 
%% Normally this would include only the scans that were just added

% save the ids of regions that may have to be used to process new FBB scans

ind_accpcc=[1002 1010 1023 1026 2002 2010 2023 2026]';
ind_front=[1003 1032 1012 1014 1018 1019 1020 1027 1028 2003 2032 2012 2014 2018 2019 2020 2027 2028]';
ind_temp=[1015 1030 2015 2030]';
ind_pariet=[1008 1025 1029 1031 2008 2025 2029 2031]';
ind_cbl=[7 8 46 47]';

%

cd (path_processed);

avfbbs=dir('*/*/FBB*');
avfbbs=strcat({avfbbs.folder}','/',{avfbbs.name}');

for i=1:size(avfbbs,1)
    
    srcfbbsuvr=dir(strcat(avfbbs{i,1},'/LDS*suvr_cbl.nii'));
    
    if size(srcfbbsuvr,1)==0
        
        srcfbbd=dir(strcat(avfbbs{i,1},'/LDS*FBB*.nii')); 
        fbbscan=strcat(srcfbbd.folder,'/',srcfbbd.name);
        rfbbscan=strcat(srcfbbd.folder,'/r',srcfbbd.name);
        
        srcnu=dir(strcat(avfbbs{i,1},'/LDS*nu.nii')); nuscan=strcat(srcnu.folder,'/',srcnu.name);
        srcaparc=dir(strcat(avfbbs{i,1},'/LDS*raparc+aseg.nii')); aparcscan=strcat(srcaparc.folder,'/',srcaparc.name);

        pathfbb=avfbbs{i,1};
        
        run LEADS_PET_FBB_pipeline.m
        
        clear srcfbbd fbbscan rfbbscan srcnu srcaparc nuscan aparcscan pathfbb
        
    end % end if condition there is no SUVR file in the folder, this means the scan was not processed
    
    clear srcfbbsuvr
    
end % end for condition for each FBB folder

cd (path_processed);
%% Let's copy the all the available JPEG files to the shared petcore
system('cp -n   */*/FBB*/Multiaxial*.pdf   /shared/petcore/Projects/LEADS/data_f7p1/summary/fbb/suvr/.');
system('cp -n   */*/FBB*/Refreg*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/fbb/refonsuvr/.');



