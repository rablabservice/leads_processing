%%%%%% Script MRI processing for LEADS 04/22/2019
%%%%%% Updated 3/6/2020 to add smoothing for the mwc1
%%%%%% Updated 6/24/2020 to handle processing of consecutive timepoints at
%%%%%% the same time and brainstem segmentation
%%%%%% Updated 9/13/2020 to handle Timepoint3 MRIs
%%%%%% Sep 2021: Code Cleanup
%%%%%% Feb 2024: Major update with the transition to 6mm PET data

% Get a list of all processed MRIs (subject ID and scan date)
processed_mris = dir(fullfile(PATHS('processed'), '**', 'MRI_T1_*', '*nu.nii*'));
processed_mris = {processed_mris.name}';
processed_mris = cellfun(@(x) x(1:end-7), processed_mris, 'uniformoutput', 0);

%% Find newly uploaded MRIs
new_mris = dir(fullfile(PATHS('newdata'), 'mri', '**', '*.nii*'));

new_mris = dir(fullfile(PATHS('newdata'), 'mri', '*'));
isSubdir = [new_mris.isdir];
notSpecialDirs = ~ismember({new_mris.name}, {'.', '..'});
new_mris = new_mris(isSubdir & notSpecialDirs);
new_mris = fullfile({new_mris.folder}', {new_mris.name}');
newmri_tab = getScanInfo()
% Get the directories one level deep and the files within them
new_mris = fullfile({new_mris.folder}', {new_mris.name}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newmris = dir('*/*/*/*/*.nii');
newmris = transpose(struct2cell(newmris));
tempnewids = regexp(newmris(:,2), 'LDS\d{7}', 'match', 'once');
tempnewmridates = regexp(newmris(:,2), '\d{4}-\d{2}-\d{2}', 'match', 'once');
listnewmris = strcat(tempnewids, '_MRI_T1_', tempnewmridates);
newmris = listnewmris;

%% We can use the newmris object for the processing
%% enter each folder in a loop to grab the MRI and rename it
%% Sep 2021: Adding the Freesurfer job filemaking here
cd(path_processed)

fname=sprintf('NewMRIProcessing_%s.txt', datestr(now,'mm-dd-yyyy_HH-MM'));
fid = fopen(fname,'wt');

fname2=sprintf('NewMRIProcessing_Brainstem_%s.txt', datestr(now,'mm-dd-yyyy_HH-MM'));
fid2 = fopen(fname2,'wt');

for i=1:size(newmris,1)
    fprintf(1,'** Now Moving and Renaming MRI %s\n', newmris{i});

    %% Grab path of the new file to be moved
    temp_id=newmris{i}(1:10);
    temp_date=newmris{i}(19:28);
    newfolder=strcat(path_newmris,temp_id);

    cd(newfolder);
    oldfname=dir(strcat('*/',temp_date,'*/*/*',temp_id,'*.nii'));
    oldfname=strcat(oldfname.folder,'/',oldfname.name);

    %% Quick module to create a new folder if needed
    profolder=char(strcat(path_processed, temp_id));
    if exist(profolder,'dir')==0
        mkdir(profolder)
    end

    %% Create required Timepoint and move the file over
    temptps=dir(strcat(profolder,'/Timepoint*'));
    temptps={temptps.name}';
    if size(temptps,1)==0 % There is no "Timepoint" folder in the processed one, this means it is a fresh subject
        mkdir(strcat(profolder,'/MRI_T1_',temp_date));
        newfname=strcat(profolder,'/MRI_T1_',temp_date,'/',temp_id,'_MRI_T1_',temp_date,'.nii');
        movefile(oldfname,newfname); %% renamed and moved the file to the output folder
    end

    %% Small module to populate Freesurfer jobs
    cmd_reconall=strcat('recon-all -all -i',{' '}, newfname,{' '}, '-sd',{' '},path_freesurfer,{' '},'-s',{' '},newmris{i});
    fprintf(fid, '%s\n', char(cmd_reconall));

    cmd_brainstem=strcat('segmentBS.sh',{' '},newmris{i},{' '},path_freesurfer);
    fprintf(fid2, '%s\n', char(cmd_brainstem));

    cd(path_newmris);
end

fclose(fid);
fclose(fid2);
cd(path_processed);

%% We just created two txt files with single lines showing the calls for freesurfer. This is done
%% so that they will be run in parallel
%% Let's run it
numc = input(['\n\n Choose how many cores (parallel jobs) (max 16): ', ...
    '\n     --> ']);
if numc>16
    fprintf(2,'Warning! You selected a number of cores >16 \nDefaulting to 16 (maximum)\n');
    numc=16;
end

disp('--> Starting Freesurfer! Come back later!');
cmd_fsf_exec=strcat('parallel -j',{' '}, num2str(numc), {' '},'-- <',{' '}, fname);
system(char(cmd_fsf_exec))

cmd_fsf_exec2=strcat('parallel -j',{' '}, num2str(numc), {' '},'-- <',{' '}, fname2);
system(char(cmd_fsf_exec2))

%% let's grab the results for each freesurfer computation
listfsf=strcat(path_freesurfer,newmris);
for ii=1:size(listfsf,1)
    if exist(listfsf{ii},'dir')==7 % The Freesurfer folder actually exists, thats good
        cd (listfsf{ii})
        cd mri

        %% Freesurfer output File conversion
        convaparc=strcat('mri_convert -it mgz -ot nii --out_orientation RAS nu.mgz nu.nii');
        convnu=strcat('mri_convert -it mgz -ot nii --out_orientation RAS aparc+aseg.mgz aparc+aseg.nii');
        convbs=strcat('mri_convert -it mgz -ot nii --out_orientation RAS brainstemSsLabels.v12.FSvoxelSpace.mgz brainstemSsLabels_v12_VoxSpace.nii');

        system(char(convaparc))
        system(char(convnu))
        system(char(convbs))

        %% recentering origin to center of mass
        % changed strcat to strvcat because script was crashing at this line
        % July 2022
        nii_setOriginMRI(strvcat('nu.nii','aparc+aseg.nii','brainstemSsLabels_v12_VoxSpace.nii'));

        % Edit August 2021 - previous steps were useless, only copying now to
        % keep filenaming consistent. All previous data is valid.
        copyfile('aparc+aseg.nii','raparc+aseg.nii');
        copyfile('brainstemSsLabels_v12_VoxSpace.nii','rbrainstemSsLabels_v12_VoxSpace.nii');

        % adding brainstem overlay and copy in petcore for the qc
        fslovrly2=strcat('overlay 1 1 nu.nii -a rbrainstemSsLabels_v12_VoxSpace.nii 170 176  qcbs_overlay.nii.gz');
        system(char(fslovrly2));
        fslovrly2_ot=strcat('slices qcbs_overlay.nii.gz -o qcbs_overlay.png');
        system(char(fslovrly2_ot))

        oldovrlyname=strcat(listfsf{ii},'/mri/qcbs_overlay.png');
        newlocovrlyname=strcat(listfsf{ii},'/mri/','qc_brainstem_',newmris{ii},'.png');
        rnmcmd=strcat('mv',{' '},oldovrlyname,{' '},newlocovrlyname);
        system(char(rnmcmd));
        newdistovrlyname=strcat('/shared/petcore/Projects/LEADS/data_f7p1/summary/mri/qc_brainstem/qc_brainstem_',newmris{ii},'.png');
        cpcmd=strcat('cp',{' '},newlocovrlyname,{' '},newdistovrlyname);
        system(char(cpcmd));

        % Small module in case this is longitudinal processing

        % new approach- if this is not the same date of the MRI at Timepoint 1,
        % coreg
        mridirs = sort({dir(strcat(path_processed,newmris{ii}(1:10),'/MRI_T1_*')).name});
        temptp1 = mridirs{1};
        if ~isequal(temptp1(8:17),newmris{ii}(19:28))
            newnupath=strcat(listfsf{ii},'/mri/','nu.nii');
            newaparcpath=strcat(listfsf{ii},'/mri/','raparc+aseg.nii');
            newbspath=strcat(listfsf{ii},'/mri/','rbrainstemSsLabels_v12_VoxSpace.nii');
            oldnupath=strcat(path_freesurfer,listfsf{ii}(size(listfsf{ii},2)-27:size(listfsf{ii},2)-18),'_',temptp1,'/mri/nu.nii');

            spm('defaults','PET');
            clear matlabbatch_coregmris
            matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.ref = cellstr(oldnupath);
            matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.source = cellstr(newnupath);
            matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.other = vertcat(cellstr(newaparcpath),cellstr(newbspath));
            matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
            matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
            matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            spm_jobman('run',matlabbatch_coregmris);
            clear matlabbatch_coregmris
        end

        % quick module to know in which Timepoint are we
        tempmytp=dir(strcat(path_processed,newmris{ii}(1:10),'/Timepoint*/MRI_T1_*'));
        tempmytp=struct2cell(tempmytp)';
        tempmytp=tempmytp(strcmp(tempmytp(:,1),newmris{ii}(12:28)),:);

        oldfname=strcat(listfsf{ii},'/mri/','nu.nii');
        newfname=strcat(tempmytp{1,2},'/',tempmytp{1,1},'/',newmris{ii},'_nu.nii');
        symlink1=strcat('ln -s',{' '},oldfname,{' '},newfname);

        oldfname2=strcat(listfsf{ii},'/mri/','raparc+aseg.nii');
        newfname2=strcat(tempmytp{1,2},'/',tempmytp{1,1},'/',newmris{ii},'_raparc+aseg.nii');
        symlink2=strcat('ln -s',{' '},oldfname2,{' '},newfname2);

        oldfname3=strcat(listfsf{ii},'/mri/','rbrainstemSsLabels_v12_VoxSpace.nii');
        newfname3=strcat(tempmytp{1,2},'/',tempmytp{1,1},'/',newmris{ii},'_rbrainstemSsLabels_v12_VoxSpace.nii');
        symlink3=strcat('ln -s',{' '},oldfname3,{' '},newfname3);

        system(char(symlink1))
        system(char(symlink2))
        system(char(symlink3))

        %%% Let's segment the 'nu' file in the MRI folder
        clear matlabbatch;
        spm('defaults','PET');
        matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(strcat(tempmytp{1,2},'/',tempmytp{1,1},'/',newmris{ii},'_nu.nii'));
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/mnt/neuroimaging/SPM/spm12/tpm/TPM.nii,1'};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/mnt/neuroimaging/SPM/spm12/tpm/TPM.nii,2'};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/mnt/neuroimaging/SPM/spm12/tpm/TPM.nii,3'};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/mnt/neuroimaging/SPM/spm12/tpm/TPM.nii,4'};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/mnt/neuroimaging/SPM/spm12/tpm/TPM.nii,5'};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/mnt/neuroimaging/SPM/spm12/tpm/TPM.nii,6'};
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.025 0.1];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
        matlabbatch{2}.spm.spatial.smooth.data = cellstr(strcat(tempmytp{1,2},'/',tempmytp{1,1},'/mwc1',newmris{ii},'_nu.nii'));
        matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{2}.spm.spatial.smooth.dtype = 0;
        matlabbatch{2}.spm.spatial.smooth.im = 0;
        matlabbatch{2}.spm.spatial.smooth.prefix = 's8iso_';
        spm_jobman('run',matlabbatch);
        clear matlabbatch
    else
        fprintf(2,'Warning! The Freesurfer folder for %s does not exist, maybe it failed?\n',newmris{ii});
    end
end

cd(path_processed);
