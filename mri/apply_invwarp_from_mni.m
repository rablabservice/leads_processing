function outfiles = apply_invwarp_from_mni(infiles, iyf, interp, vox, prefix, bb, overwrite, verbose)
    % Warp images from MNI to subject MRI space using an existing iy_ file
    %
    % Parameters
    % ----------
    % infiles : char|str or cell array of char|str
    %   Path to scan(s) to be warped to MNI space
    % iyf : char|str or cell array of char|str
    %   Path to the iy_* deformation field file created during
    %   segmentation that will now be used for warping
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
    % prefix : char|str, optional
    %   Prefix to append to the output filenames
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % verbose : logical, optional
    %   If true, print runtime information
    %
    % Files created
    % -------------
    % w<infiles{1}>
    % w<infiles{2}>
    % ...
    % ------------------------------------------------------------------
    arguments
        infiles {mustBeFile}
        iyf {mustBeFile}
        interp {mustBeMember(interp,0:7)} = 0
        vox (1,3) {mustBePositive} = [1 1 1]
        prefix {mustBeText} = 'v'
        bb (2,3) {mustBePositive} = [Inf Inf Inf; Inf Inf Inf]
        overwrite logical = false
        verbose logical = true
    end

    % Format parameters
    infiles = cellfun(@(x) abspath(x), cellstr(infiles), 'UniformOutput', false);
    iyf = cellfun(@(x) abspath(x), cellstr(iyf), 'UniformOutput', false);
    outfiles = cellfun(@(x) fullfile(fileparts(x), append(prefix, basename(x))), infiles, 'UniformOutput', false);

    % Check existence of output files
    if all(isfile(outfiles)) && ~overwrite
        if verbose
            fprintf('- Warping to subject space already complete, will not overwrite existing files\n')
        end
        return
    end

    % Run Normalise: Write (MNI -> subject)
    if verbose
        fprintf('- Warping images from MNI to subject space:\n')
        for ii = 1:length(infiles)
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
end
