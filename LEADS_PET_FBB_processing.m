%% All the new FBBs have been now moved and renamed, so it is time to process them
%% This script still scouts the whole LEADS database to try and process any scan in need.
%% Normally this would include only the scans that were just added
cd(path_processed);
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

    end % end if condition there is no SUVR file in the folder, this means the scan was not processed

end % end for condition for each FBB folder

cd(path_processed);
