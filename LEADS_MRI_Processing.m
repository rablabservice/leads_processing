%%%%%% Script MRI processing for LEADS 04/22/2019
%%%%%% Updated 3/6/2020 to add smoothing for the mwc1
%%%%%% Updated 6/24/2020 to handle processing of consecutive timepoints at
%%%%%% the same time and brainstem segmentation
%%%%%% Updated 9/13/2020 to handle Timepoint3 MRIs
%%%%%% Sep 2021: Code Cleanup
%%%%%% Feb 2024: Major update with the transition to 6mm PET data

%% Checking what subjects (MRIs) have been already processed

cd (path_processed);

listmris=dir('*/*/MRI*/LDS*nu.nii');
listmris = {listmris.name}';
listmris=cellfun(@(fnonu) fnonu(1:end-7), listmris, 'uniformoutput',0);

%% Looking for newly uploaded MRI folders

cd (path_newmris);

newmris=dir('*/*/*/*/*.nii');
newmris=transpose(struct2cell(newmris));
tempnewids=regexp(newmris(:,2),'LDS\d{7}','match','once');
tempnewmridates=regexp(newmris(:,2),'\d{4}-\d{2}-\d{2}','match','once');
listnewmris=strcat(tempnewids,'_MRI_T1_',tempnewmridates);

%% Small QC module, for a given subject there is duplicate new MRI
%% This is possible if the user inadvertedly downloaded two different sequences or two repeat scans

if size(unique(listnewmris),1)<size(listnewmris,1)

   fprintf(2,'\nWarning! There is at least one duplicated MRI in the new batch.\nTake a look below:\n\n');

    [~, uniqueIdx] =unique(listnewmris); % Find the indices of the unique strings
    duplicates = listnewmris; % Copy the original into a duplicate array
    duplicates(uniqueIdx) = []; % remove the unique strings, anything left is a duplicate
    duplicates = unique(duplicates); % find the unique duplicates
    disp(duplicates);
   error('I will stop here. Clear the workspace, cleanup data in the new MRI folder and re-run.');

end % end if condition QC there are no duplicate MRIs

%% Check whether some of the new mris were actually already stored and processed

chkmris=intersect(listnewmris,listmris);

    if size(chkmris,1)>0

        fprintf(2,'\n\nWarning! One or more MRIs have already been processed. See list below:\n\n');
        tchkmris=array2table(chkmris, ...
            'VariableNames', {'Subject'});
        disp(tchkmris);

        fprintf(2,'I will ignore and skip the corresponding newly downloaded MRIs. Press ENTER to acknowledge and continue:\n\n');
        pause;

        newmris=setdiff(listnewmris,listmris); % Removes the old processed subjects from the list of new mris

    elseif size(chkmris,1)==0

        newmris=listnewmris;

    end % end if condition some MRIs have already been processed

%% We can use the newmris object for the processing
%% enter each folder in a loop to grab the MRI and rename it
%% Sep 2021: Adding the Freesurfer job filemaking here

cd (path_processed)
weirdtps={};

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

    cd (newfolder);

    oldfname=dir(strcat('*/',temp_date,'*/*/*',temp_id,'*.nii'));
    oldfname=strcat(oldfname.folder,'/',oldfname.name);

    %% Quick module to create a new folder if needed

    profolder=char(strcat(path_processed,temp_id));

        if exist(profolder,'dir')==0
            mkdir (profolder)
        end % end if condition this is entirely new subject, in that case create the respective folder in the main processed dir

    temptps=dir(strcat(profolder,'/Timepoint*'));
    temptps={temptps.name}';

    %% Create required Timepoint and move the file over

    if size(temptps,1)==0 % There is no "Timepoint" folder in the processed one, this means it is a fresh subject

        mkdir (strcat(profolder,'/Timepoint1'));
        mkdir (strcat(profolder,'/Timepoint1/MRI_T1_',temp_date));
        newfname=strcat(profolder,'/Timepoint1/MRI_T1_',temp_date,'/',temp_id,'_MRI_T1_',temp_date,'.nii');
        movefile(oldfname,newfname); %% renamed and moved the file to the output folder

    else % there is at least one Timepoint available. Since this is a new MRI we will add sequentially

        %% Quick module to be sure nothing is wrong in the numbering of timepoints, e.g. a middle timepoint missing

    tpseq=cellfun(@(mytp) mytp(end), temptps, 'uniformoutput',0);
    tpseq=sort(tpseq);

    template_seq=(str2num(tpseq{1,1}):str2num(tpseq{size(temptps,1),1}))';

     if ~isequal(tpseq,cellstr(num2str(template_seq)))

        fprintf(2,'Warning! The Timepoint combination for %s is not as expected. I am skipping this file out of caution. \nPress ENTER to acknowledge and continue:\n\n',temp_id);
        pause
        weirdtps=vertcat(weirdtps,newmris{i});
        continue

     end % end if condition QC Timepoint numbering is not normal

        temptps=sort(temptps);
        last_tp=temptps{size(temptps,1),1}(end);
        new_tp=num2str(str2num(last_tp)+1);

        mkdir (strcat(profolder,'/Timepoint',new_tp));
        mkdir (strcat(profolder,'/Timepoint',new_tp,'/MRI_T1_',temp_date));
        newfname=strcat(profolder,'/Timepoint',new_tp,'/MRI_T1_',temp_date,'/',temp_id,'_MRI_T1_',temp_date,'.nii');
        movefile(oldfname,newfname); %% renamed and moved the file to the output folder

    end % end if condition number of timepoints available

    %% Small module to populate Freesurfer jobs

        cmd_reconall=strcat('recon-all -all -i',{' '}, newfname,{' '}, '-sd',{' '},path_freesurfer,{' '},'-s',{' '},newmris{i});
        fprintf(fid, '%s\n', char(cmd_reconall));

        cmd_brainstem=strcat('segmentBS.sh',{' '},newmris{i},{' '},path_freesurfer);
        fprintf(fid2, '%s\n', char(cmd_brainstem));

    %%

cd (path_newmris);

clear temp_id temp_date newfolder oldfname profolder temptps tpseq template_seq newfname last_tp new_tp

end % end for loop for each new MRI

fclose(fid);
fclose(fid2);
cd (path_processed);

%% We just created two txt files with single lines showing the calls for freesurfer. This is done
%% so that they will be run in parallel (max 8 at a time)
%% Let's run it

  numc = input(['\n\n Choose how many cores (parallel jobs) (max 16): ', ...
        '\n     --> ']);

    if numc>16
        fprintf(2,'Warning! You selected a number of cores >16 \nDefaulting to 16 (maximum)\n');
        numc=16;
    end % end if condition number of parallel jobs selected by the user

    disp('--> Starting Freesurfer! Come back later!');
    cmd_fsf_exec=strcat('parallel -j',{' '}, num2str(numc), {' '},'-- <',{' '}, fname);
    system(char(cmd_fsf_exec))

    cmd_fsf_exec2=strcat('parallel -j',{' '}, num2str(numc), {' '},'-- <',{' '}, fname2);
    system(char(cmd_fsf_exec2))

%% Quick module, if MRIs were skipped let's remove them %%

if size(weirdtps,1)>0

newmris=setdiff(newmris,weirdtps);

end % end if condition MRIs were skipped because of weird Timepoints

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

    fslovrly2=strcat('overlay 1 1 nu.nii -a rbrainstemSsLabels_v12_VoxSpace.nii 170 176  qcbs_overlay.nii.gz'); system(char(fslovrly2));
    fslovrly2_ot=strcat('slices qcbs_overlay.nii.gz -o qcbs_overlay.png'); system(char(fslovrly2_ot))

    oldovrlyname=strcat(listfsf{ii},'/mri/qcbs_overlay.png');
    newlocovrlyname=strcat(listfsf{ii},'/mri/','qc_brainstem_',newmris{ii},'.png');
    rnmcmd=strcat('mv',{' '},oldovrlyname,{' '},newlocovrlyname); system(char(rnmcmd));
    newdistovrlyname=strcat('/shared/petcore/Projects/LEADS/data_f7p1/summary/mri/qc_brainstem/qc_brainstem_',newmris{ii},'.png');
    cpcmd=strcat('cp',{' '},newlocovrlyname,{' '},newdistovrlyname); system(char(cpcmd));

    % Small module in case this is longitudinal processing


    % new approach- if this is not the same date of the MRI at Timepoint 1,
    % coreg


    temptp1=dir(strcat(path_processed,newmris{ii}(1:10),'/Timepoint1/MRI_T1_*'));

    if ~isequal(temptp1.name(8:17),newmris{ii}(19:28))

        newnupath=strcat(listfsf{ii},'/mri/','nu.nii');
        newaparcpath=strcat(listfsf{ii},'/mri/','raparc+aseg.nii');
        newbspath=strcat(listfsf{ii},'/mri/','rbrainstemSsLabels_v12_VoxSpace.nii');
        oldnupath=strcat(path_freesurfer,listfsf{ii}(size(listfsf{ii},2)-27:size(listfsf{ii},2)-18),'_',temptp1.name,'/mri/nu.nii');

        spm('defaults','PET');
        clear matlabbatch_coregmris
        matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.ref = cellstr(oldnupath);
        matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.source = cellstr(newnupath);
        matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.other = vertcat(cellstr(newaparcpath),cellstr(newbspath));
        matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch_coregmris{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm_jobman('run',matlabbatch_coregmris); clear matlabbatch_coregmris

    end % end if condition we are not processing the baseline

    % quick module to know in which Timepoint are we %

    tempmytp=dir(strcat(path_processed,newmris{ii}(1:10),'/Timepoint*/MRI_T1_*'));
    tempmytp=struct2cell(tempmytp)';
    tempmytp=tempmytp(strcmp(tempmytp(:,1),newmris{ii}(12:28)),:);

    %

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
        spm_jobman('run',matlabbatch); clear matlabbatch


%%%%% code to create Multi-axial views for reporting %%%%

% 1. Multislice MRI standard

    [pp,ff,~]=spm_fileparts(char(cellstr(strcat(tempmytp{1,2},'/',tempmytp{1,1},'/',newmris{ii},'_nu.nii'))));
    slovname=char(strcat(pp,'/Multiaxial_',ff,'.jpg'));
    temp_rend=char(cellstr(strcat(tempmytp{1,2},'/',tempmytp{1,1},'/',newmris{ii},'_nu.nii')));
    o = slover;
    o.cbar = 1;
    o.img(1).vol=spm_vol(temp_rend);
    o.img(1).type='structural';
    o.img(1).prop=1;
    o.transform = 'axial';
    o.figure = spm_figure('GetWin','Graphics');
    o = fill_defaults (o);
    o.slices = -30:6:58;
    o = paint(o);
    crdate=char(pp(size(pp,2)-9:end)); crid=char(ff(1:10));
    jpeglab=strcat('ID:',{' '},crid,{' '},'***',{' '},'MRI T1',{' '},'***',{' '},'LEADS.PETCORE@ucsf.edu');
    hTitAx = axes('Parent',o.figure,'Position',[0 0.97 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
    jpeglab2=strcat('Scan Date:',{' '},datestr(crdate),{' '},'***',{' '},'Quantification Date:',{' '},date);
    hTitAx2 = axes('Parent',o.figure,'Position',[0 0.95 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab2,'Parent',hTitAx2,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
    jpeglab3=strcat('3D T1 MPRAGE',{' '},'***',{' '},'Left is Left');
    hTitAx3 = axes('Parent',o.figure,'Position',[0 0.93 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab3,'Parent',hTitAx3,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
    print(slovname,'-djpeg','-r300');

    % 2. Aparc+Aseg QC

    [pp,ff,~]=spm_fileparts(char(cellstr(strcat(tempmytp{1,2},'/',tempmytp{1,1},'/',newmris{ii},'_raparc+aseg.nii'))));
    slovname2=char(strcat(pp,'/Multiaxial_axial_',ff,'.jpg'));

    temp_back=char(cellstr(strcat(tempmytp{1,2},'/',tempmytp{1,1},'/',newmris{ii},'_nu.nii')));
    temp_ovrly=char(cellstr(strcat(tempmytp{1,2},'/',tempmytp{1,1},'/',newmris{ii},'_raparc+aseg.nii')));

    o = slover;
    o.cbar = 2;
    o.img(1).vol=spm_vol(temp_back);
    o.img(1).type='structural';
    o.img(1).prop=1;
    o.img(2).vol = spm_vol (temp_ovrly);
    o.img(2).type = 'truecolour';
    o.img(2).cmap = 'actc';
    o.img(2).range = [0 2030];
    o.img(2).prop=0.8;
    o.transform = 'axial';
    o.figure = spm_figure('GetWin','Graphics');
    o = fill_defaults (o);
    o.slices = -30:6:58;
    o = paint(o);
    crdate=char(pp(size(pp,2)-9:end)); crid=char(ff(1:10));
    jpeglab=strcat('ID:',{' '},crid,{' '},'***',{' '},'Aparc on MRI T1',{' '},'***',{' '},'LEADS.PETCORE@ucsf.edu');
    hTitAx = axes('Parent',o.figure,'Position',[0 0.97 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
    jpeglab2=strcat('Scan Date:',{' '},datestr(crdate),{' '},'***',{' '},'Quantification Date:',{' '},date);
    hTitAx2 = axes('Parent',o.figure,'Position',[0 0.95 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab2,'Parent',hTitAx2,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
    jpeglab3=strcat('3D T1 MPRAGE',{' '},'***',{' '},'Left is Left');
    hTitAx3 = axes('Parent',o.figure,'Position',[0 0.93 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab3,'Parent',hTitAx3,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
    print(slovname2,'-djpeg','-r300');

    % 3. SPM Segmentation QC

    [pp,ff,~]=spm_fileparts(char(cellstr(strcat(tempmytp{1,2},'/',tempmytp{1,1},'/',newmris{ii},'_nu.nii'))));
    slovname3=char(strcat(pp,'/Multiaxial_c1nu_',ff,'.jpg'));
    temp_back=char(cellstr(strcat(tempmytp{1,2},'/',tempmytp{1,1},'/',newmris{ii},'_nu.nii')));
    temp_ovrly=char(cellstr(strcat(tempmytp{1,2},'/',tempmytp{1,1},'/c1',newmris{ii},'_nu.nii')));
    o = slover;
    o.cbar = 2;
    o.img(1).vol=spm_vol(temp_back);
    o.img(1).type='structural';
    o.img(1).prop=1;
    o.img(2).vol = spm_vol (temp_ovrly);
    o.img(2).type = 'split';
    o.img(2).cmap = 'winter';
    o.img(2).range = [0.2 0.6];

    o.transform = 'axial';
    o.figure = spm_figure('GetWin','Graphics');
    o = fill_defaults (o);
    o.slices = -30:6:58;
    o = paint(o);
    crdate=char(pp(size(pp,2)-9:end)); crid=char(ff(1:10));
    jpeglab=strcat('ID:',{' '},crid,{' '},'***',{' '},'SPM c1 on MRI T1',{' '},'***',{' '},'LEADS.PETCORE@ucsf.edu');
    hTitAx = axes('Parent',o.figure,'Position',[0 0.97 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
    jpeglab2=strcat('Scan Date:',{' '},datestr(crdate),{' '},'***',{' '},'Quantification Date:',{' '},date);
    hTitAx2 = axes('Parent',o.figure,'Position',[0 0.95 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab2,'Parent',hTitAx2,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
    jpeglab3=strcat('3D T1 MPRAGE',{' '},'***',{' '},'Left is Left');
    hTitAx3 = axes('Parent',o.figure,'Position',[0 0.93 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab3,'Parent',hTitAx3,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
    print(slovname3,'-djpeg','-r300');

    else

       fprintf(2,'Warning! %s Freesurfer folder does not exist, maybe it failed?\n',newmris{ii});

    end % end if condition freesurfer folder exists

    clear  oldovrlyname newlocovrlyname rncmd newdistovrlyname cpcmd temptp1 tempmytp oldfname newfname symlink1 oldfname2 newfname2 symlink2 oldfname3 newfname3 symlink3 pp ff slovname temp_rend slovname2 temp_back temp_ovrly slovname3

end % end for loop for each new freesurfer folder

cd (path_processed);

%% Let's copy the all the available JPEG files to the shared petcore
system('cp -n  */*/MRI*/Multiaxial_L*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/mri/t1mprage/.');
system('cp -n  */*/MRI*/Multiaxial_axial_L*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/mri/aparcont1/.');
system('cp -n  */*/MRI*/Multiaxial_c1nu_L*.jpg   /shared/petcore/Projects/LEADS/data_f7p1/summary/mri/c1ont1/.');
