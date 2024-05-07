function outfiles = segment_mri(nuf, overwrite)
    % Segment nu.nii and save forward and inverse deformations to MNI space
    %
    % Also smooths the mwc1 image by 8mm^3 FWHM to create the smwc1
    % image.
    %
    % Parameters
    % ----------
    % nuf : char or str array
    %   Path to the T1 image to segment
    % overwrite : logical, optional
    %   If true, overwrite existing files
    %
    % Returns
    % -------
    % outfiles : struct
    %   Struct array with paths to each output file
    %
    % Files created
    % -------------
    % c1<scan_tag>_nu.nii
    % c2<scan_tag>_nu.nii
    % c3<scan_tag>_nu.nii
    % c4<scan_tag>_nu.nii
    % c5<scan_tag>_nu.nii
    % iy_<scan_tag>_nu.nii
    % smwc1<scan_tag>_nu.nii
    % w<scan_tag>_nu.nii
    % y_<scan_tag>_nu.nii
    % ------------------------------------------------------------------
    arguments
        nuf {mustBeFile}
        overwrite logical = false
    end

    % Format parameters
    nuf = abspath(nuf);
    mri_dir = fileparts(nuf);
    scan_tag = get_scan_tag(mri_dir);

    % Define the output filenames
    outfiles.c1 = fullfile(mri_dir, append('c1', scan_tag, '_nu.nii'));
    outfiles.c2 = fullfile(mri_dir, append('c2', scan_tag, '_nu.nii'));
    outfiles.c3 = fullfile(mri_dir, append('c3', scan_tag, '_nu.nii'));
    outfiles.c4 = fullfile(mri_dir, append('c4', scan_tag, '_nu.nii'));
    outfiles.c5 = fullfile(mri_dir, append('c5', scan_tag, '_nu.nii'));
    outfiles.iy = fullfile(mri_dir, append('iy_', scan_tag, '_nu.nii'));
    outfiles.smwc1 = fullfile(mri_dir, append('smwc1', scan_tag, '_nu.nii'));
    outfiles.wnu = fullfile(mri_dir, append('w', scan_tag, '_nu.nii'));
    outfiles.y = fullfile(mri_dir, append('y_', scan_tag, '_nu.nii'));

    % Check if output files already exist
    if all(structfun(@(x) exist(x, 'file'), outfiles)) && ~overwrite
        fprintf('- Segmentation already complete, will not rerun\n')
        return
    else
        fprintf('- Segmenting %s\n', basename(nuf));
    end

    % Run Segment
    tpmf = fullfile(fileparts(which('spm')), 'tpm/TPM.nii');
    mustBeFile(tpmf);
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(nuf);  % input volumes to segment
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;  % bias regularizatoin
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;  % FWHM of Gaussian smoothness of bias
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];  % don't save bias-corrected images or field
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[tpmf ',1']};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;  % number of Gaussians to model intensity distribution
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];  % save native space Pr(GM), not DARTEL
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 1];  % save modulated, not unmodulated
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[tpmf ',2']};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];  % save native space Pr(WM)
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[tpmf ',3']};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];  % save native space Pr(CSF)
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[tpmf ',4']};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];  % save native space Pr(Bone)
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[tpmf ',5']};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];  % save native space Pr(Extracranial soft tissue)
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[tpmf ',6']};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;  % MRF cleanup strength
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;  % light clean
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];  % warping regularization
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';  % affine regularization to ICBM European template
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;  % smoothness fudge factor
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;  % sampling distance
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];  % save inverse and forward deformation fields
    spm_jobman('run',matlabbatch);

    % Smooth the mwc1 image
    mwc1f = fullfile(mri_dir, append('mwc1', scan_tag, '_nu.nii'))
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.smooth.data = cellstr(mwc1f);
    matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;  % same as input
    matlabbatch{1}.spm.spatial.smooth.im = 0;  % no implicit masking
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run', matlabbatch);

    % Delete the intermediary files
    seg_matf = fullfile(mri_dir, append(scan_tag, '_nu_seg8.mat'));
    if isfile(seg_matf)
        delete(seg_matf);
    end
    delete(mwc1f);
end
