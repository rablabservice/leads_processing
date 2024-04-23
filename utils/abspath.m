function outfile = abspath(infile)
    % Return the absolute path to infile as a char array, regardless of
    % whether infile exists or is passed as an absolute or relative path
    % ------------------------------------------------------------------

    % Check if the path is already absolute
    infile = deblank(infile);
    if strncmp(infile, '/', 1)
        % Use Java to normalize the path
        outfile = char(java.io.File(infile).getCanonicalPath());
    else
        % Construct full path assuming infile is relative
        fullPath = fullfile(pwd, infile);
        % Use Java to normalize the path
        outfile = char(java.io.File(fullPath).getCanonicalPath());
    end
end
