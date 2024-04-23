%% Sep 2021: Code Cleanup - Timepoint independent approach

%% Find entire list of images available that "could" be processed

% This now includes MRI, FBB SUVR, FTP SUVR, FDG SUVR
fprintf(1,'Starting the MNI Warping script...\n');

% MRI images
fprintf(1,'Looking for all available MRIs images...\n');
listmris=dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/MRI*/LDS*nu.nii');
listmris = strcat({listmris.folder}','/',{listmris.name}');

% Concatenate FBB, FTP and FDG SUVR images
fprintf(1,'Looking for all available PET SUVR images...\n');
listpets=[
    dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/FBB*/LDS*suvr_cbl.nii');
    dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/FTP*/LDS*suvr_infcblg.nii');
    dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/FDG*/LDS*suvr_pons.nii')
];
listpets = strcat({listpets.folder}','/',{listpets.name}');

% Feedback for the User
fprintf(1,'Found %d MRI and %d PET SUVR images!\n', size(listmris,1), size(listpets,1));

% Working with MRIs first
for v=1:size(listmris,1)
    tempimg=listmris{v,1};
    [p,f,e]=spm_fileparts(tempimg);

    % identify Timepoint
    temptp=p(61:70);
    if exist(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp),'dir')==0
       mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp)));
    end

    % For MRIs we have to 1) check s8iso_mwc1 is there, in that case create
    % soft link to MNI folder. 2) affine transform.

    % 1. Deal with s8iso
    s8isotempimg=strcat(p,'/s8iso_mwc1',f,e);
    if exist(s8isotempimg,'file')==0
        fprintf(2,'Warning! For %s the s8iso_mwc1* image was not found. This is very weird, I am skipping the scan!\n',f);
        continue
    elseif exist(s8isotempimg,'file')==2
        s8isolink=strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/mri/s8iso_mwc1',f,e);
        if exist(s8isolink,'file')==0
            symlink=strcat('ln -s',{' '},s8isotempimg,{' '},s8isolink); system(char(symlink)); clear symlink;
        end
    end

    % 2. Deal with affine transform
    afftempimg=strcat(p,'/a',f,e);
    if exist(afftempimg,'file')==0
        clear matlabbatch
        spm('defaults','PET');
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj(1).source = cellstr(tempimg);
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj(1).wtsrc = '';
        matlabbatch{1}.spm.tools.oldnorm.estwrite.subj(1).resample = cellstr(tempimg);
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {'/home/mac/ycobigo/neuroimaging/neuroimaging_CentOS7/SPM/spm12/toolbox/OldNorm/T1.nii,1'};
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 0;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [-100 -130 -80
                                                                 100 100 110];
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [1 1 1];
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'a';
        spm_jobman('run',matlabbatch); clear matlabbatch;

        %%% create the corresponding link
        afflink=strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/mri/affine/a',f,e);
        if exist(afflink,'file')==0
            symlink=strcat('ln -s',{' '},afftempimg,{' '},afflink); system(char(symlink)); clear symlink;
        end
    end
    clear tempimg p f e s8isotempimg s8isolink temptp afftempimg afflink
end

% Working with PET Scans next
for v=1:size(listpets,1)
    tempimg=listpets{v,1};
    [p,f,e]=spm_fileparts(tempimg);

    % identify Timepoint
    temptp=p(61:70);
    if exist(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp),'dir')==0
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp)));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/mri/')));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/fbb/')));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/ftp/')));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/fdg/')));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/mri/affine/')));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/fbb/affine/')));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/ftp/affine/')));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/fdg/affine/')));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/wmapmri/')));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/wmapftp/')));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/wmapfbb/')));
        mkdir(char(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/wmapfdg/')));
    end

    % identify Modality
    tempmod=p(72:74);

    % identify Nu file
    srcnu=dir(strcat(p,'/LDS*nu.nii'));
    tempnu=strcat(srcnu.folder,'/',srcnu.name);

    % identify y file
    srcy=dir(strcat(p(1:71),'/MRI*/y_LDS*nu.nii'));
    tempy=strcat(srcy.folder,'/',srcy.name);

    % Working on the non-linearly warped file, they are processed together
    % for PET scans so one is enough
    wtempimg=strcat(p,'/w',f,e);
    afftempimg=strcat(p,'/a',f,e);
    if exist(wtempimg,'file')==0 && exist(afftempimg,'file')==0
        clear matlabbatch
        spm('defaults','PET');

        % nonlinear warping module
        matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = cellstr(tempy);
        matlabbatch{1}.spm.util.defs.comp{1}.inv.space = cellstr(tempnu);
        matlabbatch{1}.spm.util.defs.out{1}.push.fnames = cellstr(tempimg);
        matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
        matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = cellstr(p);
        matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {'/mnt/coredata/Projects/LEADS/script_f7p1/templates/icbm152.nii'};
        matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
        matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.push.prefix = '';

        % affine transform module
        matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(1).source = cellstr(tempnu);
        matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(1).wtsrc = '';
        matlabbatch{2}.spm.tools.oldnorm.estwrite.subj(1).resample = cellstr(tempimg);
        matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.template = {'/home/mac/ycobigo/neuroimaging/neuroimaging_CentOS7/SPM/spm12/toolbox/OldNorm/T1.nii,1'};
        matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
        matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
        matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
        matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
        matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
        matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.nits = 0;
        matlabbatch{2}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
        matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
        matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.bb = [-100 -130 -80
                                                                    100 100 110];
        matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.vox = [1 1 1];
        matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
        matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{2}.spm.tools.oldnorm.estwrite.roptions.prefix = 'a';
        spm_jobman('run',matlabbatch); clear matlabbatch;

        % create links
        wlink=strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/',lower(tempmod),'/w',f,e);
        if exist(wlink,'file')==0
            symlink=strcat('ln -s',{' '},wtempimg,{' '},wlink); system(char(symlink)); clear symlink;
        end

        afflink=strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/',lower(tempmod),'/affine/a',f,e);
        if exist(afflink,'file')==0
            symlink=strcat('ln -s',{' '},afftempimg,{' '},afflink); system(char(symlink)); clear symlink;
        end
    end
    clear tempimg p f e wtempimg temptp tempmod srcnu tempnu srcy tempy afftempimg wlink
end
fprintf(1,'Finished! Enjoy your warped images!\n');
