function [atlasf, out_dir] = format_mask_inputs(atlasf, in_dir, out_dir)
    % Given a known input path, return paths to atlasf and out_dir.
    %
    % Only one input is used to infer atlasf. Order of precedence is:
    % atlasf > in_dir > out_dir. Out_dir is set to atlasf directory if
    % not specified. Raises an error if atlasf or out_dir do not exist.
    %
    % Parameters
    % ----------
    % atlasf : char or str array
    %     Path to the atlas file
    % in_dir : char or str array
    %     The input directory. If atlasf is empty, this is where the
    %     function looks for the atlas file. This parameter is
    %     disregarded if atlasf is not empty
    % out_dir : char or str array
    %     The output directory. If out_dir is empty, the mask is saved
    %     in the same directory as the atlas file
    %
    % Returns
    % -------
    % atlasf : char
    %   Nicely formatted absolute path to the atlas file
    % out_dir : char
    %   Nicely formatted absolute path to the output directory
    % ------------------------------------------------------------------
    arguments
        atlasf {mustBeText} = ''
        in_dir {mustBeText} = ''
        out_dir {mustBeText} = ''
    end

    % Format inputs
    atlasf = char(atlasf);
    in_dir = char(in_dir);
    out_dir = char(out_dir);

    % Check that the input file and output dir exist
    if isempty(atlasf)
        if isempty(in_dir)
            if isempty(out_dir)
                error('At least one of atlasf, in_dir, or out_dir must be specified');
            else
                out_dir = abspath(out_dir);
                in_dir = out_dir;
            end
        else
            in_dir = abspath(in_dir);
        end
        mri_files = get_freesurfer_files(in_dir);
        atlasf = mri_files.aparc;
    else
        atlasf = abspath(atlasf);
    end
    if ~isfile(atlasf)
        error('File not found: %s', atlasf);
    end
    if isempty(out_dir)
        out_dir = fileparts(atlasf);
    end
    if ~isfolder(out_dir)
        error('Output directory not found: %s', out_dir);
    end
end
