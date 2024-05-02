function outfiles = affine_mri_to_mni(nuf, overwrite, verbose)
    % Calculate the full affine transform from native MRI to MNI space
    %
    % Parameters
    % ----------
    % nuf : char or str array
    %   Path to the T1 image that will be used for affine estimation
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % verbose : logical, optional
    %   If true, print diagnostic information
    %
    % Files created
    % -------------
    % - <mri_dir>/a<scan_tag>_nu.nii
    % - <mri_dir>/<scan_tag>_affine-to-mni.mat (or some file like this that we can use
    %                                           later on PET files coreg'd to native MRI)
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        overwrite logical = false
        verbose logical = true
    end

    % Format paths
    infiles = cellstr(infiles);
    mri_dir = fileparts(infiles{1});
    outfiles = cell(size(infiles));


    % Check existence of input and output files
    anuf = fullfile(mri_dir, append('a', basename(nuf)));
    if exist(anuf, 'file') && ~overwrite
        return
    end

    % Calculate and affine transform to MNI space
    % interpolation: trilinear
    if verbose
        fprintf('- Calculating and applying affine transformation to %s \n',scan_tag);
    end

    clear matlabbatch;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.subj(1).source = cellstr(nuf);
    matlabbatch{1}.spm.tools.oldnorm.estwrite.subj(1).wtsrc = '';
    matlabbatch{1}.spm.tools.oldnorm.estwrite.subj(1).resample = cellstr(nuf);
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {[fileparts(which('spm')),'/toolbox/OldNorm/T1.nii,1']};
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 0;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [Inf Inf Inf
                                                         Inf Inf Inf];
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [1 1 1];
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'a';
    spm_jobman('run',matlabbatch);


    movefile( ...
        fullfile(mri_dir, append(scan_tag,'_nu_sn.mat')), ...
        fullfile(mri_dir, append(scan_tag,'_affine_to_mni.mat')) ...
    )
end
