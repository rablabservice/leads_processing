function atf = estimate_mri_affine_to_mni(nuf, templatef, overwrite)
    % Estimate affine transformation from native MRI to MNI space
    %
    % Parameters
    % ----------
    % nuf : char or str array
    %   Path to the T1 image to estimate the affine transformation for
    %
    % Returns
    % -------
    % atf : char
    %   Path to the affine transformation file
    % ------------------------------------------------------------------
    arguments
        nuf {mustBeFile}
        templatef = []
        overwrite logical = false
    end

    % Format parameters
    nuf = abspath(nuf);
    if isempty(templatef)
        templatef = fullfile(fileparts(which('spm')), '/toolbox/OldNorm/T1.nii');
    end

    % Check if the output file exists
    mri_dir = fileparts(nuf);
    scan_tag = get_scan_tag(mri_dir);
    atf = fullfile(mri_dir, append('atf_', scan_tag,'.mat'));
    if isfile(atf) && ~overwrite
        fprintf('- Affine transformation file already exists, will not re-estimate\n');
        return
    else
        fprintf('- Estimating affine transformation to MNI space\n');
    end

    % Run Old Normalise: Estimate
    clear matlabbatch;
    matlabbatch{1}.spm.tools.oldnorm.estwrite.subj(1).source = cellstr(nuf);
    matlabbatch{1}.spm.tools.oldnorm.estwrite.subj(1).wtsrc = '';
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.template = cellstr(templatef);
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.weight = '';
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.smosrc = 8;
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.smoref = 0;
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.regtype = 'mni';
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.nits = 0;
    spm_jobman('run', matlabbatch);

    % Rename the output file
    movefile(fullfile(mri_dir, append(scan_tag,'_nu_sn.mat')), atf);
end
