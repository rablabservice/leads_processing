function outfiles = format_outfiles(infiles, prefix, suffix)
    % Format outfiles to match the infiles class
    arguments
        infiles
        prefix = ''
        suffix = ''
    end

    % Return outfiles in the correct class
    outfiles = add_presuf(abspath(cellvec(infiles)), prefix, suffix);
    if isstruct(infiles)
        outfile_keys = add_presuf(fieldnames(infiles), prefix, suffix);
        outfiles = cell2struct(outfiles, outfile_keys);
    elseif ischar(infiles)
        outfiles = char(outfiles{1});
    elseif isstring(infiles) && isscalar(outfiles)
        outfiles = string(outfiles{1});
    end
end
