function outfiles = apply_affine_to_mni(infiles, atf, overwrite, interp, vox, prefix, bb)
    % Transform images from native MRI to MNI space using an existing affine
    %
    % Parameters
    % ----------
    % infiles : char or str array
    %   Path to scan(s) to be affine transformed to MNI space
    % atf : char or str array
    %   Path to the affine transform .mat file
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
    if all(isfile(outfiles)) && ~overwrite
        fprintf('- Affine transformed files exist, will not overwrite\n');
        outfiles = format_outfiles(infiles_cp, prefix);
        return
    else
        fprintf('- Affine transforming images to MNI space:\n')
        for ii = 1:length(infiles)
            if ~isempty(prefix)
                fprintf('  * %s -> %s\n', basename(infiles{ii}), basename(outfiles{ii}));
            else
                fprintf('  * %s\n', basename(infiles{ii}));
            end
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
