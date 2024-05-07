function outfiles = add_presuf(infiles, prefix, suffix)
    % Add a prefix and/or suffix to one or more filenames
    %
    % Parameters
    % ----------
    % infiles : char|string|cellstr|struct of filenames
    %     Paths to files
    % prefix : char
    %     Prefix to add to the filenames
    %
    % Returns
    % -------
    % outfiles : char|string|cellstr|struct of filenames
    %     Paths to the input files with prefix and/or suffix added, in
    %     the same order and object class as the input files
    % ------------------------------------------------------------------
    arguments
        infiles
        prefix {mustBeText} = ''
        suffix {mustBeText} = ''
    end

    function outfile = fmt_path(infile)
        % Format a single path
        [d, n, x] = fileparts(deblank(infile));
        if strcmp(x, '.gz')
            parts = strsplit(n, '.');
            n = strjoin(parts(1:end-1), '.');
            x = append('.', parts{end}, x);
        end
        outfile = fullfile(d, append(prefix, n, suffix, x));
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
