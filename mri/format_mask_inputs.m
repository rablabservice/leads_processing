function [atlasf, out_dir] = format_mask_inputs(atlasf, in_dir, out_dir)
    % Return atlasf and out_dir and raise an error if they don't exist
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
        [~, ~, atlasf, ~] = get_freesurfer_files(in_dir);
    else
        atlasf = abspath(atlasf);
    end
    if ~exist(atlasf, 'file')
        error('File not found: %s', atlasf);
    end
    if isempty(out_dir)
        out_dir = fileparts(atlasf);
    end
    if ~exist(out_dir, 'dir')
        error('Output directory not found: %s', out_dir);
    end
end
