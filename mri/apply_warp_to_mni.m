function outfiles = apply_warp_to_mni(infiles, yf, fid, overwrite, interp, vox, prefix, bb)
    % Warp images from native MRI to MNI space using an existing y_ file
    %
    % Parameters
    % ----------
    % infiles : char|string|cellstr|struct
    %   Path to scan(s) to be warped to MNI space
    % yf : char|string|cellstr|struct
    %   Path to the y_* deformation field file created during
    %   segmentation that will now be used for warping
    % fid : int, optional
    %   File identifier for logging (default is 1 for stdout)
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
    % Returns
    % -------
    % outfiles : cellstr|struct
    %   Cell array with paths to the warped images, in the same order
    %   as the input images
    %
    % Files created
    % -------------
    % w<infiles{1}>
    % w<infiles{2}>
    % ...
    % ------------------------------------------------------------------
    arguments
        infiles
        yf
        fid {mustBeNumeric} = 1
        overwrite logical = false
        interp {mustBeMember(interp,0:7)} = 4
        vox (1,3) {mustBePositive} = [1.5 1.5 1.5]
        prefix {mustBeText} = 'w'
        bb (2,3) {mustBePositive} = [Inf Inf Inf; Inf Inf Inf]
    end

    % Check that all input files exist, and format them correctly
    infiles_cp = infiles;
    infiles = abspath(cellvec(infiles));
    mustBeFile(infiles);
    yf = abspath(cellvec(yf));
    mustBeFile(yf);

    % Check existence of output files
    outfiles = add_presuf(infiles, prefix);
    if all(isfile(outfiles)) && ~overwrite
        log_append(fid, '- Warped files exist, will not overwrite');
        outfiles = format_outfiles(infiles_cp, prefix);
        return
    else
        log_append(fid, '- Warping images to MNI space:');
        for ii = 1:length(infiles)
            if ~isempty(prefix)
                msg = sprintf( ...
                    '  * %s ->\n              %s', ...
                    basename(infiles{ii}), ...
                    basename(outfiles{ii}) ...
                );
                log_append(fid, msg);
            else
                log_append(fid, sprintf('  * %s', basename(infiles{ii})));
            end
        end
    end

    % Run Normalise: Write (subject -> MNI)
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = yf;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = infiles;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vox;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = interp;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = prefix;
    spm_jobman('run', matlabbatch);

    % Format the output filenames
    outfiles = format_outfiles(infiles_cp, prefix);
end
