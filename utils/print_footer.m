function print_footer(caption, fid, indent)
    % Print a footer caption with a border before and after
    % ------------------------------------------------------------------
    arguments
        caption {mustBeText}
        fid {mustBeNumeric} = 1
        indent {mustBeNumeric} = 10
    end

    % Figure out the border lengths
    border_len = length(caption);

    % Print the starting border
    log_append(fid, '', false, 0);  % Add a blank line
    log_append(fid, repmat('_', 1, border_len), false, indent);
    log_append(fid, repmat('-', 1, border_len), false, indent);

    % Print the caption
    log_append(fid, caption, false, indent);

    % Print the ending border
    if border_len < 12
        log_append(fid, repmat('=', 1, border_len), false, indent);
        log_append(fid, '', false, 0);
        log_append(fid, '', false, 0);
    else
        log_append( ...
            fid, append('+..', repmat('=', 1, border_len - 6), '..+'), false, indent ...
        );
        log_append(fid, '', false, 0);
        log_append(fid, '', false, 0);
    end
end
