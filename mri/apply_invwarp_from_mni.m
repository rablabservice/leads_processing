function outfiles = apply_invwarp_from_mni(infiles, iyf, overwrite, interp, vox, prefix, bb)
    % Warp images from MNI to native MRI space using an existing iy_ file
    %
    % Parameters
    % ----------
    % infiles : char|string|cellstr|struct
    %   Paths to one or more scans to be warped to MNI space
    % iyf : char|string|cellstr|struct
    %   Path to the iy_* deformation field file created during
    %   segmentation that will now be used for warping
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % interp : int, optional
    %   Interpolation method (0-7):
    %     0: Nearest-neighbor
    %     1: Trilinear
    %     2: 2nd Degree B-Spline
    %     3: 3rd Degree B-Spline
    %     4: 4th Degree B-Spline
    %     5: 5th Degree B-Spline
    %     6: 6th Degree B-Spline
    %     7: 7th Degree B-Spline
    % vox : 1x3 numeric array, optional
    %   Voxel size to reslice the output files to, in mm
    % prefix : char|string, optional
    %   Prefix to append to the output filenames
    %
    % Files created
    % -------------
    % w<infiles{1}>
    % w<infiles{2}>
    % ...
    % ------------------------------------------------------------------
    arguments
        infiles
        iyf
        overwrite logical = false
        interp {mustBeMember(interp,0:7)} = 0
        vox (1,3) {mustBePositive} = [1 1 1]
        prefix {mustBeText} = 'v'
        bb (2,3) {mustBePositive} = [Inf Inf Inf; Inf Inf Inf]
    end

    % Check that all input files exist, and format them correctly
    infiles_cp = infiles;
    infiles = abspath(cellvec(infiles));
    mustBeFile(infiles);
    iyf = abspath(cellvec(iyf));
    mustBeFile(iyf);

    % Check existence of output files
    outfiles = add_presuf(infiles, prefix);
    if all(isfile(outfiles)) && ~overwrite
        fprintf('- Inverse warp files exist, will not overwrite\n')
        outfiles = format_outfiles(infiles_cp, prefix);
        return
    end

    % Run Normalise: Write (MNI -> subject)
    fprintf('- Warping images from MNI to subject space:\n')
    for ii = 1:length(infiles)
        if ~isempty(prefix)
            fprintf('  * %s -> %s\n', basename(infiles{ii}), basename(outfiles{ii}));
        else
            fprintf('  * %s\n', basename(infiles{ii}));
        end
    end
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = iyf;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = infiles;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vox;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = interp;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = prefix;
    spm_jobman('run', matlabbatch);

    % Format the output filenames
    outfiles = format_outfiles(infiles_cp, prefix);
end
