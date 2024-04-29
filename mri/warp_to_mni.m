function outfiles = warp_to_mni(infiles, iyf, overwrite, verbose)
    % Warp any file to MNI space using an existing iy_ file
    %
    % Parameters
    % ----------
    % infiles : char or str array
    %   Path to the scan that needs to be warped to MNI space
    % iyf : char or str array
    %   Path to the iy_ file with inverse deformation parameters
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % verbose : logical, optional
    %   If true, print diagnostic information
    %
    % Files created
    % -------------
    % <mri_dir>/w<infiles>
    % ------------------------------------------------------------------
    arguments
        infiles {mustBeFile}
        iyf {mustBeFile}
        overwrite logical = false
        verbose logical = true
    end

    % Format paths
    infiles = cellstr(infiles);
    outfiles = cell(size(infiles));
    for ii = 1:length(infiles)
        infile = abspath(infiles{ii});
        % Check that input file exists
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

    % Check that the iy_ file exists
    if ~exist(iyf, 'file')
        error('File not found, run segmentation first: %s', iyf)
    end

    % Apply pushforward transformation (subject -> MNI)
    if verbose
        for ii = 1:length(infiles)
            fprintf('- Warping %s to MNI space\n', basename(infiles{ii}));
        end
    end
    clear matlabbatch;
    matlabbatch{1}.spm.util.defs.comp{1}.def = cellstr(iyf);
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames = infiles;
    matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = cellstr(mri_dir);
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.bb = [Inf Inf Inf
                                                         Inf Inf Inf];
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.vox = [1 1 1];
    matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
    matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix = '';
    spm('defaults','PET');
    spm_jobman('run',matlabbatch);
end



