function atf = estimate_mri_affine_to_mni(nuf, overwrite, templatef)
    % Estimate affine transformation from native MRI to MNI space
    %
    % Parameters
    % ----------
    % nuf : char or str array
    %   Path to the T1 image to estimate the affine transformation for
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % templatef : char or str, optional
    %   Path to the template that nuf will be transformed to. Default
    %   uses the SPM12 OldNOrm T1 template
    %
    % Returns
    % -------
    % atf : char
    %   Path to the affine transformation file
    % ------------------------------------------------------------------
    arguments
        nuf {mustBeFile}
        overwrite logical = false
        templatef = []
    end

    % Format parameters
    nuf = abspath(nuf);
    if isempty(templatef)
        templatef = fullfile(fileparts(which('spm')), '/toolbox/OldNorm/T1.nii');
    end

    % Check if the output file exists
    mri_dir = fileparts(nuf);
    scan_tag = get_scan_tag(mri_dir);
    atf = fullfile(mri_dir, append('atf_', scan_tag,'_nu.mat'));
    if isfile(atf) && ~overwrite
        fprintf('- Affine transformation file exists, will not re-estimate\n');
        return
    else
        fprintf('- Estimating affine transformation to MNI space\n');
    end

    % Run Old Normalise: Estimate
    clear matlabbatch;
    matlabbatch{1}.spm.tools.oldnorm.est.subj(1).source = cellstr(nuf);
    matlabbatch{1}.spm.tools.oldnorm.est.subj(1).wtsrc = '';
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.template = cellstr(templatef);
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.weight = '';
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.smosrc = 8;
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.smoref = 0;
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.regtype = 'mni';
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.nits = 0;
    matlabbatch{1}.spm.tools.oldnorm.est.eoptions.reg = Inf;
    spm_jobman('run', matlabbatch);

    % Rename the output file
    tmpf = fullfile(mri_dir, append(scan_tag,'_nu_sn.mat'));
    movefile(tmpf, atf);
end
