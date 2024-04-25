function segment_mri(mri_dir, overwrite, verbose)
    % Segment nu.nii and save deformation files.
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
    % <mri_dir>/y_<scan_tag>_nu.nii
    % <mri_dir>/iy_<scan_tag>_nu.nii
    % <mri_dir>/smwc1<scan_tag>_nu_mni.nii
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

    % Check if segmentation/warp has already been run
    c1f = fullfile(mri_dir, append('c1', basename(nuf)));
    wnuf = fullfile(mri_dir, append('w', basename(nuf)));
    if exist(c1f, 'file') && exist(wnuf, 'file') && ~overwrite
        return
    end

    % Segment nu.nii and warp to MNI space
    if verbose
        fprintf('- Segmenting %s \n', scan_tag);
    end

    % segment using SPM12 Segment toolbox and TPM saving both inverse and
    % forward transformations
    clear matlabbatch;spm('defaults','PET');
    matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(nuf);
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[fileparts(which('spm')),'/tpm/TPM.nii,1']};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[fileparts(which('spm')),'/tpm/TPM.nii,2']};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[fileparts(which('spm')),'/tpm/TPM.nii,3']};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[fileparts(which('spm')),'/tpm/TPM.nii,4']};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[fileparts(which('spm')),'/tpm/TPM.nii,5']};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[fileparts(which('spm')),'/tpm/TPM.nii,6']};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];          
    spm_jobman('run',matlabbatch);

    % smooth modulated c1 image for further analysis
    clear matlabbatch;spm('defaults','PET');
    matlabbatch{1}.spm.spatial.smooth.data = cellstr(strcat(mri_dir,'/mwc1',scan_tag,'_nu.nii'));
    matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run',matlabbatch); clear matlabbatch
    movefile(strcat(mri_dir,'/smwc1',scan_tag,'_nu.nii'),strcat(mri_dir,'/smwc1',scan_tag,'_nu_mni.nii'));
    delete(strcat(mri_dir,'/mwc1',scan_tag,'_nu.nii'));
    delete('*.mat');
end
