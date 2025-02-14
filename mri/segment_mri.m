function outfiles = segment_mri(nuf, fid, overwrite)
    % Segment nu.nii and save forward and inverse deformations to MNI space
    %
    % Also saves the mwc1 image for use in Dartel template creation
    % and seg8.mat file for TIV calculation
    %
    % TIV volume is calculated as the sum of the GM, WM, and CSF volumes
    %
    % Parameters
    % ----------
    % nuf : char or str array
    %   Path to the T1 image to segment
    % fid : int, optional
    %   File identifier for logging (default is 1 for stdout)
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
    % mwc1<scan_tag>_nu.nii
    % y_<scan_tag>_nu.nii
    % <scan_tag>_nu_seg8.mat
    % ------------------------------------------------------------------
    arguments
        nuf {mustBeFile}
        fid {mustBeNumeric} = 1
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
    outfiles.mwc1 = fullfile(mri_dir, append('mwc1', scan_tag, '_nu.nii'));
    outfiles.wnu = fullfile(mri_dir, append('w', scan_tag, '_nu.nii'));
    outfiles.y = fullfile(mri_dir, append('y_', scan_tag, '_nu.nii'));
    outfiles.seg_matf = fullfile(mri_dir, append(scan_tag, '_nu_seg8.mat'));

    % Check if output files already exist
    if all(structfun(@(x) exist(x, 'file'), outfiles)) && ~overwrite
        log_append(fid, '- Segmentation already complete, will not rerun');
        return
    else
        log_append(fid, sprintf('- Segmenting %s', basename(nuf)));
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

    % Extract volumes for c1, c2, c3; see PMID: 25255942 for more info on methods
    clear matlabbatch;
    matlabbatch{1}.spm.util.tvol.matfiles = cellstr(fullfile(mri_dir, append(scan_tag, '_nu_seg8.mat')));
    matlabbatch{1}.spm.util.tvol.tmax = 3;
    matlabbatch{1}.spm.util.tvol.mask = cellstr(fullfile(fileparts(which('spm')), 'tpm/mask_ICV.nii,1'));
    matlabbatch{1}.spm.util.tvol.outf = fullfile(mri_dir, append(scan_tag, '_nu_seg8_TIV.csv'));
    spm_jobman('run',matlabbatch);

    % Add total volume to spreadsheet
    T = readtable(fullfile(mri_dir, append(scan_tag, '_nu_seg8_TIV.csv')));
    T.TIV_in_mL = zeros(height(T), 1);
    T.TIV_in_mL = sum(T{1, 2:4})*1000; %convert volume in litres to mL
    writetable(T, fullfile(mri_dir, append(scan_tag, '_nu_seg8_TIV.csv')));
end
