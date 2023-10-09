% Script to generate a path database for LEADS
% Leonardo.Iaccarino@ucsf.edu, Sep 2021: Code Cleanup and Timepoint
% independent processing
% Different approach to generate path database

% 1. Let's look at all the available MRIs. This also tells us where we are
% in terms of how many timepoints are available in LEADS

mrinu=dir(strcat(path_processed,'LDS*/*/MRI*/LDS*nu.nii'));
mrinu = strcat({mrinu.folder}','/',{mrinu.name}');
mritps=cellfun(@(x) x(61:70), mrinu(:,1),'uniformoutput',0);
timepoints=unique(mritps);
dbt=array2table(unique(cellfun(@(x) x(50:59), mrinu(:,1),'uniformoutput',0)));
dbt.Properties.VariableNames(1)=cellstr('ID');

for t=1:size(timepoints,1)
    
% 1. MRI section

% 1a. NU files, find and merge with outer join

tempnu=dir(strcat(path_processed,'LDS*/',timepoints{t,1},'/MRI*/LDS*nu.nii'));
tempnu=strcat({tempnu.folder}','/',{tempnu.name}');
tempnu=array2table(horzcat(cellfun(@(x) x(50:59), tempnu(:,1),'uniformoutput',0),tempnu));
tempnu.Properties.VariableNames(1)=cellstr('ID');
tempnu.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_MRI_nu'));
dbt=outerjoin(dbt,tempnu,'MergeKeys',true);

% 1b. APARC files, find and merge with outer join

tempaparc=dir(strcat(path_processed,'LDS*/',timepoints{t,1},'/MRI*/LDS*raparc+aseg.nii'));
tempaparc=strcat({tempaparc.folder}','/',{tempaparc.name}');
tempaparc=array2table(horzcat(cellfun(@(x) x(50:59), tempaparc(:,1),'uniformoutput',0),tempaparc));
tempaparc.Properties.VariableNames(1)=cellstr('ID');
tempaparc.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_MRI_aparc'));
dbt=outerjoin(dbt,tempaparc,'MergeKeys',true);

% 1c. s8iso files, find and merge with outer join

temps8iso=dir(strcat(path_processed,'LDS*/',timepoints{t,1},'/MRI*/s8iso*LDS*nu.nii'));
temps8iso=strcat({temps8iso.folder}','/',{temps8iso.name}');
temps8iso=array2table(horzcat(cellfun(@(x) x(50:59), temps8iso(:,1),'uniformoutput',0),temps8iso));
temps8iso.Properties.VariableNames(1)=cellstr('ID');
temps8iso.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_MRI_s8iso_mwc1'));
dbt=outerjoin(dbt,temps8iso,'MergeKeys',true);

% 1d. waffine files, find and merge with outer join

tempwaff=dir(strcat(path_processed,'LDS*/',timepoints{t,1},'/MRI*/w_affineLDS*nu.nii'));
tempwaff=strcat({tempwaff.folder}','/',{tempwaff.name}');
tempwaff=array2table(horzcat(cellfun(@(x) x(50:59), tempwaff(:,1),'uniformoutput',0),tempwaff));
tempwaff.Properties.VariableNames(1)=cellstr('ID');
tempwaff.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_MRI_waffnu'));
dbt=outerjoin(dbt,tempwaff,'MergeKeys',true);

% 1e. wmap files, find and merge with outer join

tempwmap=dir(strcat(path_processed,'LDS*/',timepoints{t,1},'/MRI*/Wmap_s8iso*LDS*nu.nii'));
tempwmap=strcat({tempwmap.folder}','/',{tempwmap.name}');
tempwmap=array2table(horzcat(cellfun(@(x) x(50:59), tempwmap(:,1),'uniformoutput',0),tempwmap));
tempwmap.Properties.VariableNames(1)=cellstr('ID');
tempwmap.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_MRI_wmap'));
dbt=outerjoin(dbt,tempwmap,'MergeKeys',true);

clear tempnu tempaparc temps8iso tempwaff tempwmap

% 2. PET images - create another loop to make things easier

mods={'FBB','FTP','FDG'}; refs={'cbl','infcblg','pons'};

     for m=1:size(mods,2)

    %2a. SUVR

    tempsuvr=dir(strcat(path_processed,'LDS*/',timepoints{t,1},'/',mods{m},'*/LDS*suvr_',refs{m},'.nii'));
    tempsuvr=strcat({tempsuvr.folder}','/',{tempsuvr.name}');
    
    if size(tempsuvr,1)==0
    tempsuvr=array2table(unique(cellfun(@(x) x(50:59), mrinu(:,1),'uniformoutput',0)));
    tempsuvr.Empty(:)=NaN;
    tempsuvr.Properties.VariableNames(1)=cellstr('ID');
    tempsuvr.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_',mods{m},'_suvr'));
    else
    tempsuvr=array2table(horzcat(cellfun(@(x) x(50:59), tempsuvr(:,1),'uniformoutput',0),tempsuvr));
    tempsuvr.Properties.VariableNames(1)=cellstr('ID');
    tempsuvr.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_',mods{m},'_suvr'));
    end % end if condition Images were found
    
    dbt=outerjoin(dbt,tempsuvr,'MergeKeys',true);

    %2b. wSUVR

    tempwsuvr=dir(strcat(path_processed,'LDS*/',timepoints{t,1},'/',mods{m},'*/wLDS*suvr_',refs{m},'.nii'));
    tempwsuvr=strcat({tempwsuvr.folder}','/',{tempwsuvr.name}');
    
    if size(tempwsuvr,1)==0
    tempwsuvr=array2table(unique(cellfun(@(x) x(50:59), mrinu(:,1),'uniformoutput',0)));
    tempwsuvr.Empty(:)=NaN;
    tempwsuvr.Properties.VariableNames(1)=cellstr('ID');
    tempwsuvr.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_',mods{m},'_wsuvr'));
    else
    tempwsuvr=array2table(horzcat(cellfun(@(x) x(50:59), tempwsuvr(:,1),'uniformoutput',0),tempwsuvr));
    tempwsuvr.Properties.VariableNames(1)=cellstr('ID');
    tempwsuvr.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_',mods{m},'_wsuvr'));
    end % end if condition Images were found
    
    dbt=outerjoin(dbt,tempwsuvr,'MergeKeys',true);

    %2c. waffSUVR

    tempwaffsuvr=dir(strcat(path_processed,'LDS*/',timepoints{t,1},'/',mods{m},'*/w_affineLDS*suvr_',refs{m},'.nii'));
    tempwaffsuvr=strcat({tempwaffsuvr.folder}','/',{tempwaffsuvr.name}');
    
    if size(tempwaffsuvr,1)==0
    tempwaffsuvr=array2table(unique(cellfun(@(x) x(50:59), mrinu(:,1),'uniformoutput',0)));
    tempwaffsuvr.Empty(:)=NaN;
    tempwaffsuvr.Properties.VariableNames(1)=cellstr('ID');
    tempwaffsuvr.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_',mods{m},'_waffsuvr'));
    else
    tempwaffsuvr=array2table(horzcat(cellfun(@(x) x(50:59), tempwaffsuvr(:,1),'uniformoutput',0),tempwaffsuvr));
    tempwaffsuvr.Properties.VariableNames(1)=cellstr('ID');
    tempwaffsuvr.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_',mods{m},'_waffsuvr'));
    end % end if condition Images were found
    
    dbt=outerjoin(dbt,tempwaffsuvr,'MergeKeys',true);

    %2d. Wmap

    tempwmap=dir(strcat(path_processed,'LDS*/',timepoints{t,1},'/',mods{m},'*/Wmap_wLDS*suvr_',refs{m},'.nii'));
    tempwmap=strcat({tempwmap.folder}','/',{tempwmap.name}');
    if size(tempwmap,1)==0
    tempwmap=array2table(unique(cellfun(@(x) x(50:59), mrinu(:,1),'uniformoutput',0)));
    tempwmap.Empty(:)=NaN;
    tempwmap.Properties.VariableNames(1)=cellstr('ID');
    tempwmap.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_',mods{m},'_wmap'));
    else
    tempwmap=array2table(horzcat(cellfun(@(x) x(50:59), tempwmap(:,1),'uniformoutput',0),tempwmap));
    tempwmap.Properties.VariableNames(1)=cellstr('ID');
    tempwmap.Properties.VariableNames(2)=cellstr(strcat(timepoints{t,1},'_Path_',mods{m},'_wmap'));
    end % end if condition Images were found

    dbt=outerjoin(dbt,tempwmap,'MergeKeys',true);

    clear tempsuvr tempwsuvr tempwaffsuvr tempwmap

     end % end for loop for each PET modality

end % end for loop for each timepoints

filename = sprintf('LEADS_DatabasePaths_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
writetable(dbt,char(strcat(path_dbs,filename)),'WriteRowNames',false);

%% create a symbolic link in the shared petcore for easier access

newfname=strcat('/shared/petcore/Projects/LEADS/data_f7p1/dbs/',filename);
copyfile(strcat(path_dbs,filename),newfname);

%%

fprintf(1,'Paths Database successfully created!\n');
clear;