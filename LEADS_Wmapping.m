%% Wmapping Script - Sep 2021 Code Cleanup

%% store variables we will need later 

path_processed='/mnt/coredata/Projects/LEADS/data_f7p1/processed/';
path_wscore_meta='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/meta/';

itcpt_mri='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/mri_agetiv/beta_0001.nii';
betaage_mri='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/mri_agetiv/beta_0002.nii';
betativ_mri='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/mri_agetiv/beta_0003.nii';
sdmap_mri='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/mri_agetiv/SDmap.nii';

itcpt_fbb='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/fbb_age/beta_0001.nii';
betaage_fbb='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/fbb_age/beta_0002.nii';
sdmap_fbb='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/fbb_age/SDmap.nii';

itcpt_ftp='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/ftp_age/beta_0001.nii';
betaage_ftp='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/ftp_age/beta_0002.nii';
sdmap_ftp='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/ftp_age/SDmap.nii';

itcpt_fdg='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/fdg_age/beta_0001.nii';
betaage_fdg='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/fdg_age/beta_0002.nii';
sdmap_fdg='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/calibration/fdg_age/SDmap.nii';

emask='/mnt/coredata/Projects/LEADS/script_f7p1/templates/EM_final.nii';

%% Warn the user of the approach

fprintf(1,'The Wmapping script now looks for the LEADS_PTDEMOG.csv file in the petcore drive,\nin /shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/\nPress enter to acknowledge and continue:\n');
pause;

% new approach, look for the demographic database in the shared/petcore
% drive, keeping the approach to grab the latest if multiple are present

dbage=dir('/shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/LEADS_PTDEMOG*');
dbage=struct2cell(dbage)';
dbage=dbage(:,[1 3]);
dbage=array2table(dbage);
dbage.Properties.VariableNames(1) = cellstr(strcat('Filename'));
dbage.Properties.VariableNames(2) = cellstr(strcat('DateCreated'));
dbage.DateCreated=datetime(dbage.DateCreated);
filt=max(dbage.DateCreated); % Store most recent date
dbage=dbage(dbage.DateCreated==filt,:); % Subset to select the most recent file 
clear filt;

% some data handling to get to the same database format we used to have,
% since it worked

dbagef=readtable(strcat('/shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/',char(dbage{1,1}))); %% ready to check differences now
dbagef=dbagef(:,[1 7]);
dbagef.Properties.VariableNames(1) = cellstr(strcat('ID'));
dbagef.Properties.VariableNames(2) = cellstr(strcat('dob'));
dbagef.ID=cellfun(@upper,dbagef.ID,'UniformOutput',false);
dbagef.dob=datetime(dbagef.dob);

% Our database is now ready

% We have to make an extra step for the MRI Wmapping, we need the TIVs
% before computing

% Get the list of available seg8.mat files

srcseg8s=dir(strcat(path_processed,'LDS*/*/MRI*/*seg8.mat'));
allseg8s = strcat({srcseg8s.folder}','/',{srcseg8s.name}');

% Prepare for later

allseg8s=array2table(horzcat(allseg8s, regexp(allseg8s, 'LDS\d{7}_\w{3}_\w{2}_\d{4}-\d{2}-\d{2}','match','once')));
allseg8s.Properties.VariableNames(1) = cellstr(strcat('File'));
allseg8s.Properties.VariableNames(2) = cellstr(strcat('Identifier'));
allseg8s=sortrows(allseg8s,'Identifier');

% For most of these MRIs, we already have the TIV available since it was
% previously calculated, so now we want to just find who is new

%%% We now have the database with seg8 files ready

olddbtiv=dir(strcat(path_wscore_meta,'LEADS_TIV*')); % Time to check for existing databases, so we can run new TIV calculations only on new stuff

olddbtiv=struct2cell(olddbtiv)';
olddbtiv=olddbtiv(:,[1 3]);
olddbtiv=array2table(olddbtiv);
olddbtiv.Properties.VariableNames(1) = cellstr(strcat('Filename'));
olddbtiv.Properties.VariableNames(2) = cellstr(strcat('DateCreated'));
olddbtiv.DateCreated=datetime(olddbtiv.DateCreated);
filt=max(olddbtiv.DateCreated); % Store most recent date
olddbtiv=olddbtiv(olddbtiv.DateCreated==filt,:); % Subset to select the most recent file
oldinfotiv=readtable(strcat(path_wscore_meta,char(table2cell(olddbtiv(1,1)))),'Delimiter',','); %% ready to check differences now
oldinfotiv=sortrows(oldinfotiv,'Identifier');

checkinfo=isequal(allseg8s(:,{'Identifier'}),oldinfotiv(:,{'Identifier'})); % Check for overlap between identifiers
if checkinfo==0
    
    Indexa=not(ismember(allseg8s(:,{'Identifier'}),oldinfotiv(:,{'Identifier'}))); %% logical indexing for the new cases. 
    newcases=allseg8s(Indexa,:); %% Creates a table with the New MRI seg8 scans
    
    if size(newcases,1)>0

        fprintf(1,'**New MRI seg8 files were found:\n\n');
        disp(newcases.Identifier);
        newcases=table2array(newcases);
        newcases=newcases(:,1); % prepare for SPM

        clear matlabbatch
        spm('defaults','PET');
        matlabbatch{1}.spm.util.tvol.matfiles = newcases;
        matlabbatch{1}.spm.util.tvol.tmax = 3;
        matlabbatch{1}.spm.util.tvol.mask = {'/mnt/neuroimaging/SPM/spm12/tpm/mask_ICV.nii,1'};
        matlabbatch{1}.spm.util.tvol.outf = char(strcat(path_wscore_meta,'tempTIV.csv'));
        spm_jobman('run',matlabbatch); clear matlabbatch;
        
        % read the new values, compute the TIV
        dbtiv=readtable(char(strcat(path_wscore_meta,'tempTIV.csv')),'Delimiter',',');
        dbtiv.TIV=sum(dbtiv{:,2:end},2);
        dbtiv.Identifier=regexp(dbtiv.File, 'LDS\d{7}_\w{3}_\w{2}_\d{4}-\d{2}-\d{2}','match','once');
        
        %append to previous database, remove temporary file
        dbtiv=vertcat(oldinfotiv,dbtiv);
        dbtiv=sortrows(dbtiv,'Identifier');
        filenamedbtivs = sprintf('LEADS_TIVs_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
        writetable(dbtiv,char(strcat(path_wscore_meta,filenamedbtivs)),'WriteRowNames',false);
        delete(char(strcat(path_wscore_meta,'tempTIV.csv')));
        
    end % end if condition there are new MRIs for which we need TIVs
    
elseif checkinfo==1
    
    fprintf(2,'**There are no MRI scans to update\n\n');
    dbtiv=oldinfotiv;
    
end % end if condition there are no differences between old and new seg8 file lists

%%% Ready to Wmap!

%First thing: let's look in our "processed" subjects to see what is
%ready to be wmapped

listimgs=[dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/MRI*/s8iso*LDS*nu.nii');dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/FBB*/wLDS*suvr_cbl.nii');dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/FTP*/wLDS*suvr_infcblg.nii');dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/FDG*/wLDS*suvr_pons.nii')];
listimgs = strcat({listimgs.folder}','/',{listimgs.name}');

% Now let's work on the images

for v=1:size(listimgs,1)
    
  tempid=listimgs{v,1}(50:59);
  temptp=listimgs{v,1}(61:70);
  tempmod=listimgs{v,1}(72:74);

  if isequal(tempmod,'MRI')
  tempdate=datetime(listimgs{v,1}(79:88));
  else
  tempdate=datetime(listimgs{v,1}(76:85));
  end % end if condition what modality to grab date

    [p,f,e]=spm_fileparts(listimgs{v,1});
    
    % 1. First, check if the image already exists.
    
    if exist(strcat(p,'/Wmap_',f,e),'file')==0
        
        % 2. check that dob/age is available
        
        qcdob=datetime(dbagef(ismember(dbagef.ID, tempid),:).dob, 'Format','yyyy-MM-dd');
        
        if size(qcdob,1)==1 && qcdob<datetime('1978-01-01') && qcdob>datetime('1953-01-01')
            
        % we need to calculate the Wmap, create sets of files and
        % expression depending on the modality
            
        tempimg=listimgs{v,1};
        tempage=calyears(between(qcdob,tempdate,'years'));
        tempwmapname=strcat(p,'/Wmap_',f,e);
        
            if isequal(tempmod,'MRI')
            tempidentifier=cellstr(regexp(tempimg, 'LDS\d{7}_\w{3}_\w{2}_\d{4}-\d{2}-\d{2}','match','once'));
            temptiv=dbtiv(ismember(dbtiv.Identifier,tempidentifier),:).TIV;
            tempvols=vertcat(cellstr(tempimg),itcpt_mri,betaage_mri,betativ_mri,sdmap_mri);
            tempexpr=strcat('(i1-(i2+(i3*',num2str(tempage),')+(i4*',num2str(temptiv),')))./i5');
            elseif isequal(tempmod,'FBB')
            tempvols=vertcat(cellstr(tempimg),itcpt_fbb,betaage_fbb,sdmap_fbb);
            tempexpr=strcat('(i1-(i2+(i3*',num2str(tempage),')))./i4'); 
            elseif isequal(tempmod,'FTP')
            tempvols=vertcat(cellstr(tempimg),itcpt_ftp,betaage_ftp,sdmap_ftp);
            tempexpr=strcat('(i1-(i2+(i3*',num2str(tempage),')))./i4'); 
            elseif isequal(tempmod,'FDG')
            tempvols=vertcat(cellstr(tempimg),itcpt_fdg,betaage_fdg,sdmap_fdg);
            tempexpr=strcat('(i1-(i2+(i3*',num2str(tempage),')))./i4'); 
            end % end if condition with which modality we are working
            
        tempvols2=vertcat(cellstr(tempwmapname),cellstr(emask));
        tempwmapname_gm=char(strcat(p,'/Wmap_',f,'_GM',e));

        % Regardless of modality, we are ready to compute
           
        clear matlabbatch
        spm('defaults','PET'); 
        matlabbatch{1}.spm.util.imcalc.input = tempvols;
        matlabbatch{1}.spm.util.imcalc.output = tempwmapname;
        matlabbatch{1}.spm.util.imcalc.outdir = '';
        matlabbatch{1}.spm.util.imcalc.expression =tempexpr;
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        matlabbatch{2}.spm.util.imcalc.input = tempvols2;
        matlabbatch{2}.spm.util.imcalc.output = tempwmapname_gm;
        matlabbatch{2}.spm.util.imcalc.outdir = '';
        matlabbatch{2}.spm.util.imcalc.expression = 'i1.*i2';
        matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{2}.spm.util.imcalc.options.mask = 0;
        matlabbatch{2}.spm.util.imcalc.options.interp = 1;
        matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch); clear matlabbatch;
            
        % Now create link
        
        symlink1=strcat('ln -s',{' '},tempwmapname,{' '},strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/wmap',lower(tempmod),'/Wmap_',f,e)); system(char(symlink1));
        
        % Generate Multiaxial
        
        run LEADS_Wmap_MultiAxial_service.m
        
        % Generate 3D rendering
        
        run LEADS_Wmap_3DRend_service.m
            
        elseif size(qcdob,1)==1 && (qcdob>datetime('1978-01-01') || qcdob<datetime('1953-01-01'))
        
            fprintf(2,'Warning! For %s the dob was %s, it is not in the 40-65 LEADS range.\nI am skipping the subject. Press Enter to acknowledge and continue:\n',tempid,char(qcdob));
            pause
            
        end
    end % end if condition existence of the Wmap in the folder
    
    clear tempid temptp tempmod tempdate p f e qcdob tempimg tempage tempwmapname tempidentifier temptiv tempvols tempexpr tempvols2 tempwmapname_gm symlink1
    
end % end for loop for each image

cd (path_processed);
%% Let's copy the all the available JPEG files to the shared petcore
system('cp -n   **/**/FTP*/Wmap*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/ftp/wmap/.');
system('cp -n   **/**/FBB*/Wmap*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/fbb/wmap/.');
system('cp -n   **/**/FDG*/Wmap*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/fdg/wmap/.');
system('cp -n   **/**/MRI*/Wmap*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/mri/wmap/.');

system('cp -n   **/**/FTP*/3DRend_W*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/ftp/wmap_3d/.');
system('cp -n   **/**/FBB*/3DRend_W*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/fbb/wmap_3d/.');
system('cp -n   **/**/FDG*/3DRend_W*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/fdg/wmap_3d/.');
system('cp -n   **/**/MRI*/3DRend_W*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/mri/wmap_3d/.');

disp('Finished, enjoy your Wmaps!');
clear;
