function coivox = mri_reset_origin_com(vols)
    % Recenter T1 to center-of-mass then coregister to the SPM template
    %
    % Strictly rigid-body coregistration (no reslicing is done).
    %
    % vols : char array of image file names
    %     Images to reorient, assumed to be from same individual and
    %     session. The first image should be the T1 MRI, and the
    %     transform is estimated only for this image. Same transform is
    %     applied to any additional images.
    % ------------------------------------------------------------------
    if verbose
        fprintf('- Resetting origin to center-of-mass\n')
        for v = 1:size(vols, 1)
            filepath = abspath(vols(v, :));
            if v == 1
                fprintf('  - Source image: %s\n', basename(filepath));
            elseif v == 2
                fprintf('    Other images: %s\n', basename(filepath));
            else
                fprintf('                  %s\n', basename(filepath));
            end
        end
    end

    % Check that all input vols exist
    for v = 1:size(vols, 1)
        if ~exist(abspath(vols(v, :)), 'file')
            error('File not found: %s', abspath(vols(v, :)));
        end
    end

    % Initialize SPM
    spm_jobman('initcfg');

    % Load the first image
    coivox = ones(4,1);
    [pth, nam, ext] = spm_fileparts(deblank(vols(1,:)));
    fname = fullfile(pth, [nam ext]);
    hdr = spm_vol([fname,',1']);
    img = spm_read_vols(hdr);
    img = img - min(img(:));
    img(isnan(img)) = 0;

    % Find the center-of-mass along each dimension
    sumTotal = sum(img(:));
    coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal;
    coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal;
    coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal;
    XYZ_mm = hdr.mat * coivox; % convert from voxels to mm

    % Update the origin in the header of each image
    for v = 1:size(vols, 1)
        fname = deblank(vols(v, :));
        if ~isempty(fname)
            [pth, nam, ext] = spm_fileparts(fname);
            fname = fullfile(pth, [nam ext]);
            hdr = spm_vol([fname ',1']);
            fname = fullfile(pth, [nam '.mat']);
            if exist(fname, 'file')
                destname = fullfile(pth, [nam '_old.mat']);
                copyfile(fname, destname);
            end
            hdr.mat(1,4) =  hdr.mat(1,4) - XYZ_mm(1);
            hdr.mat(2,4) =  hdr.mat(2,4) - XYZ_mm(2);
            hdr.mat(3,4) =  hdr.mat(3,4) - XYZ_mm(3);
            spm_create_vol(hdr);
            if exist(fname, 'file')
                delete(fname);
            end
        end
    end

    % Coregister images to the SPM template
    if verbose
        fprintf('- Coregistering %s to the SPM OldNorm/T1.nii template\n', basename(vols(1, :)));
    end
    coregSub(vols);
    for v = 1:size(vols,1)
        [pth, nam, ~] = spm_fileparts(deblank(vols(v, :)));
        fname = fullfile(pth,[nam '.mat']);
        if exist(fname, 'file')
            delete(fname);
        end
    end
end

function coregSub(vols)
    % Coregister vols to the T1 template
    template = fullfile(spm('Dir'), 'toolbox', 'OldNorm', 'T1.nii');
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {template};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {[deblank(vols(1,:)),',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = cellstr(vols(2:end,:));
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm_jobman('run', matlabbatch);
end