function outfiles = abspath(infiles)
    % Return the absolute paths for one or more filepaths.
    %
    % Parameters
    % ----------
    % infiles : char|string|cellstr|struct of filenames
    %    Paths to files or directories
    %
    % Returns
    % -------
    % outfiles : char|string|cellstr|struct of filenames
    %     Absoute paths for the input files, in the same order and
    %     object class as the input files
    % ------------------------------------------------------------------
    function outfile = fmt_path(infile)
        % Format a single path
        infile = deblank(infile);
        if strncmp(infile, '/', 1)
            % Use Java to normalize the path
            outfile = java.io.File(infile).getCanonicalPath();
        else
            % Construct full path assuming infile is relative
            fullPath = fullfile(pwd, infile);
            % Use Java to normalize the path
            outfile = java.io.File(fullPath).getCanonicalPath();
        end
        if ischar(infile)
            outfile = char(outfile);
        elseif isstring(infile)
            outfile = string(outfile);
        end
    end

    % Return outfiles in the correct format
    if ischar(infiles)
        outfiles = fmt_path(infiles);
    elseif isstring(infiles)
        if isscalar(infiles)
            outfiles = fmt_path(infiles);
        else
            outfiles = arrayfun(@(x) fmt_path(x), infiles, 'UniformOutput', false);
        end
    elseif iscell(infiles)
        outfiles = cellfun(@(x) fmt_path(x), infiles, 'UniformOutput', false);
    elseif isstruct(infiles)
        outfiles = structfun(@(x) fmt_path(x), infiles, 'UniformOutput', false);
    else
        error('Input must be a char, string, cellstr, or struct')
    end
end
