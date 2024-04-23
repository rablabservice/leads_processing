%%% Script to update and process PET scans
%%% Let's find already available PET folders %%%
cd(path_processed);

% list of available scans %
listfbbs=dir('*/*/FBB*/LDS*suvr_cbl.nii');
listfbbs = {listfbbs.name}';
listfbbs=cellfun(@(fnosuvrcbl) fnosuvrcbl(1:end-13), listfbbs, 'uniformoutput',0);

listftps=dir('*/*/FTP*/LDS*suvr_infcblg.nii');
listftps = {listftps.name}';
listftps=cellfun(@(fnosuvrinfcblg) fnosuvrinfcblg(1:end-17), listftps, 'uniformoutput',0);

listfdgs=dir('*/*/FDG*/LDS*suvr_pons.nii');
listfdgs = {listfdgs.name}';
listfdgs=cellfun(@(fnosuvrpons) fnosuvrpons(1:end-14), listfdgs, 'uniformoutput',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% FBB Section %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find all the new available FBBs
cd(path_newfbbs);
newfbbs=dir('*/*/*/*/*.nii');
newfbbs=transpose(struct2cell(newfbbs));
tempnewids=regexp(newfbbs(:,2),'LDS\d{7}','match','once');
tempnewfbbdates=regexp(newfbbs(:,2),'\d{4}-\d{2}-\d{2}','match','once');
listnewfbbs=strcat(tempnewids,'_FBB_',tempnewfbbdates);

%% Check whether some of the new fbbs were actually already stored and processed
chkfbbs=intersect(listnewfbbs,listfbbs);
if size(chkfbbs,1)>0
    fprintf(2,'\n\nWarning! One or more FBBs have already been processed. See list below:\n\n');
    tchkfbbs=array2table(chkfbbs, ...
        'VariableNames', {'Subject'});
    disp(tchkfbbs);
    fprintf(2,'I will ignore and skip the corresponding newly downloaded FBBs. Press ENTER to acknowledge and continue:\n\n');
    pause;
    newfbbs=setdiff(listnewfbbs,listfbbs); % Removes the old processed subjects from the list of new fbbs
elseif size(chkfbbs,1)==0
    newfbbs=listnewfbbs;
end

%% Let's copy the new FBB scans we are interested in
% we have to be careful about dates here to match the scan correctly
% Sep 2021 change: code is now unlinked to number of timepoints
for i=1:size(newfbbs,1)
    temp_id=newfbbs{i}(1:10);
    temp_date=newfbbs{i}(16:25);

    % let's grab path to old file we will have to move first
    oldfname=dir(strcat(temp_id,'/*/',temp_date,'*/*/*',temp_id,'*.nii'));
    oldfname=strcat(oldfname.folder,'/',oldfname.name);

    % now let's find all available timepoints and respective dates
    temptps=dir(strcat(path_processed,temp_id,'/*/MRI*'));
    if size(temptps,1)>0
        temptps=strcat({temptps.folder}','/',{temptps.name}');
        temptps(:,2)=cellfun(@(grabtp) grabtp(61:70), temptps(:,1), 'uniformoutput',0);
        temptps(:,3)=cellfun(@(grabdate) grabdate(end-9:end), temptps(:,1), 'uniformoutput',0);
        temptps(:,4)=num2cell(abs(caldays(between(datetime(temp_date),datetime(temptps(:,3)),'days'))));
        temptps(:,5)=num2cell(cell2mat(temptps(:,4))<300);
        temptps(:,6)=cellfun(@(x) petexist(temp_id,x,'FBB'), temptps(:,2), 'uniformoutput',0);
        temptps=temptps(cell2mat(temptps(:,5)),:); % selecting only timepoints <300 days
        temptps=temptps(cell2mat(temptps(:,6))==0,:); % selecting only timepoints without FBB
        temptps=temptps(cell2mat(temptps(:,4))==min(cell2mat(temptps(:,4))),:); % if multiple still meet criteria, pick the closest to pet
        if size(temptps,1)==1
            mkdir(strcat(path_processed,temp_id,'/',temptps{1,2},'/FBB_',temp_date));
            newfname=strcat(path_processed,temp_id,'/',temptps{1,2},'/FBB_',temp_date,'/',temp_id,'_FBB_',temp_date,'.nii');
            movefile(oldfname,newfname); %% renamed and moved the file to the output folder

            %%% we also want to create the symbolic links for the MRI files (nu &
            %%% aparc+aseg) that will be used later for FBB processing.
            srcnu=dir(strcat(temptps{1,1},'/LDS*nu.nii'));
            nulnk=strcat(temptps{1,1},'/',srcnu.name);
            srcaparc=dir(strcat(temptps{1,1},'/LDS*raparc+aseg.nii'));
            aparclnk=strcat(temptps{1,1},'/',srcaparc.name);

            nucmd=strcat('cp -a',{' '},nulnk,{' '}, strcat(path_processed,temp_id,'/',temptps{1,2},'/FBB_',temp_date));
            system(char(nucmd));
            aparccmd=strcat('cp -a',{' '},aparclnk,{' '},strcat(path_processed,temp_id,'/',temptps{1,2},'/FBB_',temp_date));
            system(char(aparccmd));

            fprintf(1,'** Florbetaben scan %s was added to %s without warnings.\n', newfbbs{i}, temptps{1,2});
        else
            fprintf(2,'Warning! There was not a unique timepoint to which assign the new FBB scan done on %s for %s. I will skip the subject!\n', temp_date, temp_id);
        end
    else
        fprintf(2,'Warning! There were no MRIs available for %s. I will skip the subject!\n',temp_id);
    end
end

run LEADS_PET_FBB_processing.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% FTP Section %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find all the new available FTPs
cd(path_newftps);
newftps=dir('*/*/*/*/*.nii');
newftps=transpose(struct2cell(newftps));
tempnewids=regexp(newftps(:,2),'LDS\d{7}','match','once');
tempnewftpdates=regexp(newftps(:,2),'\d{4}-\d{2}-\d{2}','match','once');
listnewftps=strcat(tempnewids,'_FTP_',tempnewftpdates);

%% Check whether some of the new ftps were actually already stored and processed
chkftps=intersect(listnewftps,listftps);
if size(chkftps,1)>0
    fprintf(2,'\n\nWarning! One or more FTPs have already been processed. See list below:\n\n');
    tchkftps=array2table(chkftps, ...
        'VariableNames', {'Subject'});
    disp(tchkftps);
    fprintf(2,'I will ignore and skip the corresponding newly downloaded FTPs. Press ENTER to acknowledge and continue:\n\n');
    pause;
    newftps=setdiff(listnewftps,listftps); % Removes the old processed subjects from the list of new ftps
elseif size(chkftps,1)==0
    newftps=listnewftps;
end

%% Let's copy the new FTP scans we are interested in
% we have to be careful about dates here to match the scan correctly
% Sep 2021 change: code is now unlinked to number of timepoints
for i=1:size(newftps,1)
    temp_id=newftps{i}(1:10);
    temp_date=newftps{i}(16:25);

    % let's grab path to old file we will have to move first
    oldfname=dir(strcat(temp_id,'/*/',temp_date,'*/*/*',temp_id,'*.nii'));
    oldfname=strcat(oldfname.folder,'/',oldfname.name);

    % now let's find all available timepoints and respective dates
    temptps=dir(strcat(path_processed,temp_id,'/*/MRI*'));
    if size(temptps,1)>0
        temptps=strcat({temptps.folder}','/',{temptps.name}');
        temptps(:,2)=cellfun(@(grabtp) grabtp(61:70), temptps(:,1), 'uniformoutput',0);
        temptps(:,3)=cellfun(@(grabdate) grabdate(end-9:end), temptps(:,1), 'uniformoutput',0);
        temptps(:,4)=num2cell(abs(caldays(between(datetime(temp_date),datetime(temptps(:,3)),'days'))));
        temptps(:,5)=num2cell(cell2mat(temptps(:,4))<300);
        temptps(:,6)=cellfun(@(x) petexist(temp_id,x,'FTP'), temptps(:,2), 'uniformoutput',0);
        temptps=temptps(cell2mat(temptps(:,5)),:); % selecting only timepoints <300 days
        temptps=temptps(cell2mat(temptps(:,6))==0,:); % selecting only timepoints without FTP
        temptps=temptps(cell2mat(temptps(:,4))==min(cell2mat(temptps(:,4))),:); % if multiple still meet criteria, pick the closest to pet
        if size(temptps,1)==1
             mkdir(strcat(path_processed,temp_id,'/',temptps{1,2},'/FTP_',temp_date));
             newfname=strcat(path_processed,temp_id,'/',temptps{1,2},'/FTP_',temp_date,'/',temp_id,'_FTP_',temp_date,'.nii');
             movefile(oldfname,newfname); %% renamed and moved the file to the output folder

            %%% we also want to create the symbolic links for the MRI files (nu &
            %%% aparc+aseg) that will be used later for FTP processing.
            srcnu=dir(strcat(temptps{1,1},'/LDS*nu.nii'));
            nulnk=strcat(temptps{1,1},'/',srcnu.name);
            srcaparc=dir(strcat(temptps{1,1},'/LDS*raparc+aseg.nii'));
            aparclnk=strcat(temptps{1,1},'/',srcaparc.name);

            nucmd=strcat('cp -a',{' '},nulnk,{' '}, strcat(path_processed,temp_id,'/',temptps{1,2},'/FTP_',temp_date));
            system(char(nucmd));
            aparccmd=strcat('cp -a',{' '},aparclnk,{' '},strcat(path_processed,temp_id,'/',temptps{1,2},'/FTP_',temp_date));
            system(char(aparccmd));

            % Store the reverse normalized matrix
            srcrevnorm=dir(strcat(temptps{1,1},'/iy_*'));
            revnormlnk=strcat(temptps{1,1},'/',srcrevnorm.name);
            revnormcmd=strcat('ln -s',{' '},revnormlnk, {' '}, strcat(path_processed,temp_id,'/',temptps{1,2},'/FTP_',temp_date,'/y_revnorm.nii'));
            system(char(revnormcmd));

            fprintf(1,'** Flortaucipir scan %s was added to %s without warnings.\n', newftps{i}, temptps{1,2});
        else
            fprintf(2,'Warning! There was not a unique timepoint to which assign the new FTP scan done on %s for %s. I will skip the subject!\n', temp_date, temp_id);
        end
    else
        fprintf(2,'Warning! There were no MRIs available for %s. I will skip the subject!\n',temp_id);
    end
end

run LEADS_PET_FTP_processing.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% FDG Section %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find all the new available FDGs
cd(path_newfdgs);
newfdgs=dir('*/*/*/*/*.nii');
newfdgs=transpose(struct2cell(newfdgs));
tempnewids=regexp(newfdgs(:,2),'LDS\d{7}','match','once');
tempnewfdgdates=regexp(newfdgs(:,2),'\d{4}-\d{2}-\d{2}','match','once');
listnewfdgs=strcat(tempnewids,'_FDG_',tempnewfdgdates);

%% Check whether some of the new fdgs were actually already stored and processed
chkfdgs=intersect(listnewfdgs,listfdgs);
if size(chkfdgs,1)>0
    fprintf(2,'\n\nWarning! One or more FDGs have already been processed. See list below:\n\n');
    tchkfdgs=array2table(chkfdgs, ...
        'VariableNames', {'Subject'});
    disp(tchkfdgs);
    fprintf(2,'I will ignore and skip the corresponding newly downloaded FDGs. Press ENTER to acknowledge and continue:\n\n');
    pause;
    newfdgs=setdiff(listnewfdgs,listfdgs); % Removes the old processed subjects from the list of new fdgs
elseif size(chkfdgs,1)==0
    newfdgs=listnewfdgs;
end

%% Let's copy the new FDG scans we are interested in
% we have to be careful about dates here to match the scan correctly
% Sep 2021 change: code is now unlinked to number of timepoints
for i=1:size(newfdgs,1)
    temp_id=newfdgs{i}(1:10);
    temp_date=newfdgs{i}(16:25);

    % let's grab path to old file we will have to move first
    oldfname=dir(strcat(temp_id,'/*/',temp_date,'*/*/*',temp_id,'*.nii'));
    oldfname=strcat(oldfname.folder,'/',oldfname.name);

    % now let's find all available timepoints and respective dates
    temptps=dir(strcat(path_processed,temp_id,'/*/MRI*'));
    if size(temptps,1)>0
        temptps=strcat({temptps.folder}','/',{temptps.name}');
        temptps(:,2)=cellfun(@(grabtp) grabtp(61:70), temptps(:,1), 'uniformoutput',0);
        temptps(:,3)=cellfun(@(grabdate) grabdate(end-9:end), temptps(:,1), 'uniformoutput',0);
        temptps(:,4)=num2cell(abs(caldays(between(datetime(temp_date),datetime(temptps(:,3)),'days'))));
        temptps(:,5)=num2cell(cell2mat(temptps(:,4))<300);
        temptps(:,6)=cellfun(@(x) petexist(temp_id,x,'FDG'), temptps(:,2), 'uniformoutput',0);
        temptps=temptps(cell2mat(temptps(:,5)),:); % selecting only timepoints <300 days
        temptps=temptps(cell2mat(temptps(:,6))==0,:); % selecting only timepoints without FDG
        temptps=temptps(cell2mat(temptps(:,4))==min(cell2mat(temptps(:,4))),:); % if multiple still meet criteria, pick the closest to pet
        if size(temptps,1)==1
            mkdir(strcat(path_processed,temp_id,'/',temptps{1,2},'/FDG_',temp_date));
            newfname=strcat(path_processed,temp_id,'/',temptps{1,2},'/FDG_',temp_date,'/',temp_id,'_FDG_',temp_date,'.nii');
            movefile(oldfname,newfname); %% renamed and moved the file to the output folder

            %%% we also want to create the symbolic links for the MRI files (nu &
            %%% aparc+aseg) that will be used later for FDG processing.
            srcnu=dir(strcat(temptps{1,1},'/LDS*nu.nii'));
            nulnk=strcat(temptps{1,1},'/',srcnu.name);
            srcaparc=dir(strcat(temptps{1,1},'/LDS*raparc+aseg.nii'));
            aparclnk=strcat(temptps{1,1},'/',srcaparc.name);
            srcbs=dir(strcat(temptps{1,1},'/LDS*rbrainstemSsLabels_v12_VoxSpace.nii'));
            bslnk=strcat(temptps{1,1},'/',srcbs.name);

            nucmd=strcat('cp -a',{' '},nulnk,{' '}, strcat(path_processed,temp_id,'/',temptps{1,2},'/FDG_',temp_date));
            system(char(nucmd));
            aparccmd=strcat('cp -a',{' '},aparclnk,{' '},strcat(path_processed,temp_id,'/',temptps{1,2},'/FDG_',temp_date));
            system(char(aparccmd));
            bscmd=strcat('cp -a',{' '},bslnk,{' '},strcat(path_processed,temp_id,'/',temptps{1,2},'/FDG_',temp_date));
            system(char(bscmd));

            fprintf(1,'** FDG scan %s was added to %s without warnings.\n', newfdgs{i}, temptps{1,2});
        else
            fprintf(2,'Warning! There was not a unique timepoint to which assign the new FDG scan done on %s for %s. I will skip the subject!\n', temp_date, temp_id);
        end
    else
        fprintf(2,'Warning! There were no MRIs available for %s. I will skip the subject!\n',temp_id);
    end
end

run LEADS_PET_FDG_processing.m
