%%% Master script to deal with EOnonAD M36 amyloid PET within LEADS project.
%%% Nidhi Mundada, UCSF; nidhi.mundada@ucsf.edu
%%% Derived from LEADS Screening Script
%%% April 2022

%%% Restore maltlab defaults %% Nidhi May 2022
restoredefaultpath;
addpath('/home/mac/nmundada/matlab');
addpath('/mnt/neuroimaging/FreeSurfer/V7.1.0/matlab');
addpath('/mnt/neuroimaging/FreeSurfer/V7.1.0/fsfast/toolbox');
addpath(genpath('/mnt/coredata/Projects/Resources/scripts/leotools/')) ;

%%% Addpath for something we will need later
addpath('/mnt/neuroimaging/SPM/spm12/');
addpath('/mnt/coredata/Projects/LEADS/screening_M36/service/');

%%% Create some variables we will use later, dirs template and ROIs
scrdir='/shared/petcore/Projects/LEADS/screening_M36/scr_folders/';
dwndir='/mnt/coredata/Projects/LEADS/screening_M36/new_downloads/';
db_file = '/mnt/coredata/Projects/LEADS/screening_M36/Database_LEADS.csv';
fbbtemplate = '/mnt/coredata/Projects/LEADS/screening_M36/service/TemplateFBB_ADNI72_s8mm.nii';
rois={'/mnt/coredata/Projects/LEADS/screening_M36/service/ROI/Cingulate.nii', ...
    '/mnt/coredata/Projects/LEADS/screening_M36/service/ROI/Frontal.nii', ...
    '/mnt/coredata/Projects/LEADS/screening_M36/service/ROI/Parietal.nii', ...
    '/mnt/coredata/Projects/LEADS/screening_M36/service/ROI/Temporal.nii', ...
    '/mnt/coredata/Projects/LEADS/screening_M36/service/ROI/WholeCerebellum.nii'}';


%%% Load the table
db=readtable(db_file); % Reads the database that should have been updated after the last analysis

%%% Get list of images already screened
allimgs=dir(strcat(scrdir,'LDS*')); allimgs={allimgs.name}';

%%% Check for new files to be screened %%%
newimgs=dir(strcat(dwndir,'LDS*')); newimgs={newimgs.name}';

%%% Get the subset of images that were uploaded and indeed have not been
%%% screened yet
newimgs = setdiff(newimgs, allimgs);

% Process new scans if any
if isempty(newimgs)==0

newdb={};

    for i=1:size(newimgs,1)

    % module to convert DICOM to NIFTI since LONI change
    % June 2023; nidhi.mundada@ucsf.edu
    dcmdir=dir(strcat(dwndir, newimgs{i,1},'/*/*/*/*.dcm'));
    dcmdir=dcmdir.folder;
    cd (dcmdir)
    dicm2nii(pwd,pwd,'3D.nii;')

    cd('/mnt/coredata/Projects/LEADS/screening_M36/')

    % Create new output directory
    newscrdir=char(strcat(scrdir,cellstr(newimgs{i}))); mkdir(newscrdir);

    % grab the image
    oldfname=dir(strcat(dwndir,newimgs{i},'/*/*/*/*.nii')); oldfname=strcat(oldfname.folder,'/',oldfname.name);

    % move and rename file
    temp_img=strcat(newscrdir,'/',newimgs{i},'.nii');
    movefile(oldfname,temp_img); %% renamed and moved the file to the output folder

    % Create dir in which to store the quantification result
    mkdir(strcat(newscrdir,'/quantification/results'))

    % attempt first coreg to template by first set origin to center of mass
    % nii_setOrigin(temp_img); % method failed
    % used method used in other scripts to reset origin % June 2023 change
    file = deblank(temp_img);
    st.vol = spm_vol(file);
    vs = st.vol.mat\eye(4);
    vs(1:3,4) = (st.vol.dim+1)/2;
    spm_get_space(st.vol.fname,inv(vs));

    % Warp the image with the FBBPET template
    clear matlabbatch;
    spm('defaults','PET');
    matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = cellstr(temp_img);
    matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
    matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = cellstr(temp_img);
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = cellstr(fbbtemplate);
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [-90 -130 -80
                                                             90 90 90];
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [2 2 2];
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
    spm_jobman('run',matlabbatch); clear matlabbatch;

    % Create the Warped Image object and import
    new_wimg=strcat(newscrdir,'/w',newimgs{i},'.nii');
    img=spm_vol(new_wimg);
    imgv=spm_read_vols(img);

    % extract values from ROIs
    for f = 1:size(rois,1) % Looping through the ROIs
      roi=spm_vol(rois{f,:}); % Reading img header from the loop
      roiv=spm_read_vols(roi);
      img_mask=imgv.*roiv; % Creating the masked image
      temp=nonzeros(img_mask);
      temp=temp(~isnan(temp));
      ext_val=mean(temp);
      vals(f)=ext_val; % Store values for the individual ROIs across the images
      clear roi roiv img_mask temp ext_val
    end % end for loop for each macrolobar ROI needed for screening

    %%%% Store results of extraction and compute mean and suvr, plus amyloid
    %%%% positivity based on quantification

    tvals=cell2table(num2cell(vals),'VariableNames',{'cingulate' 'frontal' 'parietal' 'temporal' 'wholecerebellum'});
    tvals.composite=mean(tvals{:,1:4},2);
    tvals.composite_suvr=tvals.composite./tvals.wholecerebellum;
    tvals.amypos_quant(tvals.composite_suvr >= 1.18) = 1; % Cut-off threshold based on the new ADNI72 Validation

    %%%%% Create a table to store substrings of filenames
    dems=table(newimgs(i),cellstr(newimgs{i}(4:6)),cellstr(newimgs{i}(7:end)),'VariableNames',{'filename' 'site' 'atri_id'});

    %%%%%% Append to new database %%%%%%%

    demstvals=[dems tvals]; %% combine filename,site and id with calculations
    newdb=vertcat(newdb,demstvals);

    %%% create quantification output %%%
    acq_date=strsplit(oldfname,'/'); acq_date=acq_date{10}(1:10);
    quant=table(newimgs(i),cellstr(acq_date),tvals.composite_suvr, 1.18, tvals.amypos_quant,'VariableNames',{'filename' 'Date_FBB' 'Composite_Score','Threshold','Amy_Positivity'});
    quantfile=strcat(newscrdir,'/quantification/results/',newimgs{i},'_score.csv');
    writetable(quant,quantfile) %% Write estimation for the reads

    %%%%%%%%%%% Create warped SUVR images %%%%%%%%%%%%%

    suvr_wfname=char(strcat(newscrdir,'/SUVR_w',cellstr(newimgs{i})));
    exp=char(strcat('i1./',num2str(tvals.wholecerebellum(1))));

    clear matlabbatch;
    spm('defaults','PET');
    matlabbatch{1}.spm.util.imcalc.input = cellstr(new_wimg);
    matlabbatch{1}.spm.util.imcalc.output = suvr_wfname;
    matlabbatch{1}.spm.util.imcalc.outdir = {''};
    matlabbatch{1}.spm.util.imcalc.expression = exp;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = spm_type('float32');
    spm_jobman('run',matlabbatch); clear matlabbatch;

    suvr_fname=char(strcat(newscrdir,'/SUVR_',cellstr(newimgs{i})));

    clear matlabbatch;
    spm('defaults','PET');
    matlabbatch{1}.spm.util.imcalc.input = cellstr(temp_img);
    matlabbatch{1}.spm.util.imcalc.output = suvr_fname;
    matlabbatch{1}.spm.util.imcalc.outdir = {''};
    matlabbatch{1}.spm.util.imcalc.expression = exp;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = spm_type('float32');
    spm_jobman('run',matlabbatch); clear matlabbatch;

    clear dicmdir newscrdir temp_img oldfname newfname new_wimg img imgv tvals dems demstvals acq_date quant quantfile suvr_wfname exp suvr_fname

    end % end for loop for each new image

% Write the updated database

db.site=cellstr(num2str(db.site));
db.atri_id=cellstr(num2str(db.atri_id));

finaltable=vertcat(db,newdb);%% merge with existing info in the database
writetable(finaltable,'/mnt/coredata/Projects/LEADS/screening_M36/Database_LEADS.csv') %% overwrite Database

% Cleanup the new downloads folder

tobedeleted=dir(strcat(dwndir,'LDS*'));
tobedeleted = strcat({tobedeleted.folder}','/',{tobedeleted.name}');

    for i=1:size(tobedeleted,1)
    rmdir(tobedeleted{i},'s')
    end % end for loop for each folder in the newly downloaded dir

% Feedback to the user

fprintf(1,'Done. Processed %d new image(s) for LEADS screening.\n.',size(newimgs,1));

clear;
close all;

else

    clear;
    close all
    fprintf(1,'No new images were found. Analyses are up-to-date!\n');

end % end if condition there are images that need to be screened

% Close matlab
exit;
