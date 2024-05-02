function outfiles = apply_affine_to_mni(infiles, atf, prefix, overwrite, verbose)
    % Warp images from subject MRI to MNI space using an existing y_ file
    %
    % Parameters
    % ----------
    % infiles : char or str array
    %   Path to scan(s) to be affine transformed to MNI space
    % atf : char or str array
    %   Path to the atf_ file with inverse deformation parameters
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % verbose : logical, optional
    %   If true, print diagnostic information
    %
    % Files created
    % -------------
    % <mri_dir>/w<infiles{1}>
    % <mri_dir>/w<infiles{2}>
    % ...
    % ------------------------------------------------------------------
    arguments
        infiles {mustBeFile}
        atf {mustBeFile}
        interp {mustBeMember(interp,0:7)} = 0
        vox (1,3) {mustBePositive} = [1 1 1]
        prefix {mustBeText} = 'a'
        bb (2,3) {mustBePositive} = [Inf Inf Inf; Inf Inf Inf]
        overwrite logical = false
        verbose logical = true
    end

    % Format parameters
    infiles = cellfun(@(x) abspath(x), cellstr(infiles), 'UniformOutput', false);
    atf = cellfun(@(x) abspath(x), cellstr(atf), 'UniformOutput', false);
    outfiles = cellfun(@(x) fullfile(fileparts(x), append(prefix, basename(x))), infiles, 'UniformOutput', false);

    infiles = cellstr(infiles);
    mri_dir = fileparts(infiles{1});
    outfiles = cell(size(infiles));

    % Check existence of input and output files
    for ii = 1:length(infiles)
        infile = abspath(infiles{ii});
        if ~exist(infile, 'file')
            error('File not found: %s', infile)
        end
        outfiles{ii} = fullfile(fileparts(infile), append('w', basename(infile)));
        if exist(outfiles{ii}, 'file') && ~overwrite
            if verbose
                fprintf('File already exists, will not overwrite: %s\n', outfiles{ii})
            end
            return
        end
    end

    % Run Old Normalise: Write (subject -> MNI)
    if verbose
        for ii = 1:length(infiles)
            fprintf('- Warping %s to MNI space\n', basename(infiles{ii}));
        end
    end
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
end
