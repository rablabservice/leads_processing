function outfiles = apply_affine_to_mni(infiles, atf, fid, overwrite, interp, vox, prefix, bb)
    % Transform images from native MRI to MNI space using an existing affine
    %
    % Parameters
    % ----------
    % infiles : char or str array
    %   Path to scan(s) to be affine transformed to MNI space
    % atf : char or str array
    %   Path to the affine transform .mat file
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
    % Files created
    % -------------
    % a<infiles{1}>
    % a<infiles{2}>
    % ...
    % ------------------------------------------------------------------
    arguments
        infiles
        atf
        fid {mustBeNumeric} = 1
        overwrite logical = false
        interp {mustBeMember(interp,0:7)} = 4
        vox (1,3) {mustBePositive} = [1.5 1.5 1.5]
        prefix {mustBeText} = 'a'
        bb (2,3) {mustBePositive} = [Inf Inf Inf; Inf Inf Inf]
    end

    % Check that all input files exist, and format them correctly
    infiles_cp = infiles;
    infiles = abspath(cellvec(infiles));
    mustBeFile(infiles);
    atf = abspath(cellvec(atf));
    mustBeFile(atf);

    % Check existence of output files
    outfiles = add_presuf(infiles, prefix);
    outfile_exists = isfile(outfiles);
    if any(outfile_exists) && ~overwrite
        if all(outfile_exists)
            log_append(fid, '- Affine transformed files exist, will not overwrite');
            outfiles = format_outfiles(infiles_cp, prefix);
            return
        else
            log_append(fid, '- Some affine transformed files exist, will not overwrite them');
            infiles = infiles(~outfile_exists);
            outfiles = outfiles(~outfile_exists);
        end
    end

    % Log the files that will be warped
    log_append(fid, '- Affine transforming images to MNI space:');
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

    % Run Old Normalise: Write (subject -> MNI)
    clear matlabbatch;
    matlabbatch{1}.spm.tools.oldnorm.write.subj(1).matname = atf;
    matlabbatch{1}.spm.tools.oldnorm.write.subj(1).resample = infiles;
    matlabbatch{1}.spm.tools.oldnorm.write.roptions.preserve = 0;
    matlabbatch{1}.spm.tools.oldnorm.write.roptions.bb = bb;
    matlabbatch{1}.spm.tools.oldnorm.write.roptions.vox = vox;
    matlabbatch{1}.spm.tools.oldnorm.write.roptions.interp = interp;
    matlabbatch{1}.spm.tools.oldnorm.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.tools.oldnorm.write.roptions.prefix = prefix;
    spm_jobman('run',matlabbatch);

    % Format the output filenames
    outfiles = format_outfiles(infiles_cp, prefix);
end
