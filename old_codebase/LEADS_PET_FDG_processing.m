%% All the new FDGs have been now moved and renamed, so it is time to process them

cd (path_processed);

avfdgs=dir('*/*/FDG*');
avfdgs=strcat({avfdgs.folder}','/',{avfdgs.name}');

for i=1:size(avfdgs,1)
    
    srcfdgsuvr=dir(strcat(avfdgs{i,1},'/LDS*suvr_pons.nii'));
    
    if size(srcfdgsuvr,1)==0
        
        srcfdgd=dir(strcat(avfdgs{i,1},'/LDS*FDG*.nii')); 
        fdgscan=strcat(srcfdgd.folder,'/',srcfdgd.name);
        rfdgscan=strcat(srcfdgd.folder,'/r',srcfdgd.name);
        
        srcnu=dir(strcat(avfdgs{i,1},'/LDS*nu.nii')); nuscan=strcat(srcnu.folder,'/',srcnu.name);
        srcbs=dir(strcat(avfdgs{i,1},'/LDS*rbrainstemSsLabels_v12_VoxSpace.nii')); bsscan=strcat(srcbs.folder,'/',srcbs.name);
        
        pathfdg=avfdgs{i,1};
        
        run LEADS_PET_FDG_pipeline.m
        
        clear srcfdgd fdgscan rfdgscan srcnu srcbs nuscan bsscan pathfdg
        
    end % end if condition there is no SUVR file in the folder, this means the scan was not processed
    
    clear srcfdgsuvr
    
end % end for condition for each FDG folder

cd (path_processed);
%% Let's copy the all the available JPEG files to the shared petcore
system('cp -n   */*/FDG*/Multiaxial*.pdf   /shared/petcore/Projects/LEADS/data_f7p1/summary/fdg/suvr/.');
system('cp -n   */*/FDG*/Refreg*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/fdg/refonsuvr/.');

