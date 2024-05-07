function apply_warp_to_mni(yf, infiles, bb, voxmm, interp, prefix, overwrite, verbose)
    % Warp images from native MRI to MNI space using an existing y_ file
    %
    arguments
        bb = [Inf Inf Inf; Inf Inf Inf];
        voxmm = 1.5;
        interp = 4;
        prefix = 'w';
        overwrite = false;
        verbose = true;
    end

    clear matlabbatch;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = yf;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = infiles;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [voxmm voxmm voxmm];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = interp;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = prefix;
    spm_jobman('run', matlabbatch);
end
