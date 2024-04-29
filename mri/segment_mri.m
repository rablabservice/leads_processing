function outfiles = segment_mri(mri_dir, overwrite, verbose)
    % Segment nu.nii and save forward and inverse deformation files
    %
    % Also smooths the mwc1 image by 8mm^3 FWHM to create the smwc1
    % image.
    %
    % Parameters
    % ----------
    % mri_dir : char or str array
    %   The directory that contains the processed MRI data
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % verbose : logical, optional
    %   If true, print diagnostic information
    %
    % Files created
    % -------------
    % <mri_dir>/c1<scan_tag>_nu.nii
    % <mri_dir>/c2<scan_tag>_nu.nii
    % <mri_dir>/c3<scan_tag>_nu.nii
    % <mri_dir>/iy_<scan_tag>_nu.nii
    % <mri_dir>/smwc1<scan_tag>_nu.nii
    % <mri_dir>/y_<scan_tag>_nu.nii
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        overwrite logical = false
        verbose logical = true
    end

    % Format paths
    mri_dir = abspath(mri_dir);
    scan_tag = get_scan_tag(mri_dir);
    [~, nuf] = get_freesurfer_files(mri_dir, 'nii');

    % Get the output filenames
    outfiles.c1 = fullfile(mri_dir, append('c1', scan_tag, '_nu.nii'));
    outfiles.c2 = fullfile(mri_dir, append('c2', scan_tag, '_nu.nii'));
    outfiles.c3 = fullfile(mri_dir, append('c3', scan_tag, '_nu.nii'));
    outfiles.iy = fullfile(mri_dir, append('iy_', scan_tag, '_nu.nii'));
    outfiles.smwc1 = fullfile(mri_dir, append('smwc1', scan_tag, '_nu.nii'));
    outfiles.wnu = fullfile(mri_dir, append('w', scan_tag, '_nu.nii'));
    outfiles.y = fullfile(mri_dir, append('y_', scan_tag, '_nu.nii'));

    % Check if output files already exist
    if all(structfun(@(x) exist(x, 'file'), outfiles)) && ~overwrite
        return
    end

    % Check that the template TPM.nii file exists
    tpmf = fullfile(fileparts(which('spm')), 'tpm/TPM.nii');
    if ~exist(tpmf, 'file')
        error('Could not find TPM.nii in SPM directory')
    end

    % Segment nu.nii and save the forward and inverse deformation files
    if verbose
        fprintf('- Segmenting %s\n', basename(nuf));
    end
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(nuf);
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];  % don't save bias-corrected images or field
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[tpmf ',1']};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];  % save native space, not DARTEL
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 1];  % save modulated, not unmodulated
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[tpmf ',2']};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[tpmf ',3']};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[tpmf ',4']};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[tpmf ',5']};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[tpmf ',6']};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;  % light clean
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';  % ICBM European
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];  % save fwd and inv defs
    spm('defaults','PET');
    spm_jobman('run',matlabbatch);

    % smooth the mwc1 image
    clear matlabbatch;
    mwc1f = fullfile(mri_dir, append('mwc1', scan_tag, '_nu.nii'))
    matlabbatch{1}.spm.spatial.smooth.data = cellstr(mwc1f);
    matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;  % same as input
    matlabbatch{1}.spm.spatial.smooth.im = 0;  % no implicit masking
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm('defaults','PET');
    spm_jobman('run',matlabbatch);

    % Delete intermediary files
    delete(mwc1f);
    delete(fullfile(mri_dir, '*.mat'));
end
