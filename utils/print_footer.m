function print_footer(caption, fid, indent)
    % Print a footer caption with a border before and after
    % ------------------------------------------------------------------
    arguments
        caption {mustBeText}
        fid {mustBeNumeric} = 1
        indent {mustBeNumeric} = 0
    end

    % Define defaults
    add_timestamp = false;

    % Figure out the border lengths
    border_len = length(caption);

    % Print the starting border
    log_append(fid, '', 0, 0);
    log_append(fid, repmat('_', 1, border_len), add_timestamp, indent);
    log_append(fid, repmat('-', 1, border_len), add_timestamp, indent);

    % Print the caption
    log_append(fid, caption, add_timestamp, indent);

    % Print the ending border
    if border_len < 12
        log_append(fid, repmat('=', 1, border_len), add_timestamp, indent);
        log_append(fid, '', 0, 0);
    else
        log_append( ...
            fid, append('+..', repmat('=', 1, border_len - 6), '..+'), add_timestamp, indent ...
        );
        log_append(fid, '', 0, 0);
    end
end
