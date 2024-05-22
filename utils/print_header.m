function print_header(header, subtitle, fid, indent)
    % Print a header and optional subtitle with a border before and after
    % ------------------------------------------------------------------
    arguments
        header {mustBeText}
        subtitle {mustBeText} = ''
        fid {mustBeNumeric} = 1
        indent {mustBeNumeric} = 0
    end

    % Define defaults
    add_timestamp = false;

    % Figure out the border lengths
    border_len = max(length(header), length(subtitle));

    % Print the starting border
    if border_len < 12
        topbot = append(repmat('=', 1, border_len));
    else
        topbot = append('+..', repmat('=', 1, border_len - 6), '..+');
    end
    log_append(fid, '', add_timestamp, 0);
    log_append(fid, topbot, add_timestamp, indent);
    log_append(fid, repmat('-', 1, border_len), add_timestamp, indent);

    % Print the header and subtitle
    if isempty(subtitle)
        if border_len < 12
            log_append(fid, header, add_timestamp, indent);
        else
            log_append(fid, '', add_timestamp, 0);
            log_append(fid, header, add_timestamp, indent);
            log_append(fid, '', add_timestamp, 0);
        end
    else
        pad = floor((border_len - min(length(header), length(subtitle))) / 2);
        if length(header) > length(subtitle)
            subtitle = append(repmat(' ', 1, pad), subtitle, repmat(' ', 1, pad));
        elseif length(header) < length(subtitle)
            header = append(repmat(' ', 1, pad), header, repmat(' ', 1, pad));
        end
        log_append(fid, header, add_timestamp, indent);
        log_append(fid, '', add_timestamp, indent);
        log_append(fid, subtitle, add_timestamp, indent);
    end

    % Print the ending border
    log_append(fid, repmat('_', 1, border_len), add_timestamp, indent);
    log_append(fid, topbot, add_timestamp, indent);
    log_append(fid, '', add_timestamp, 0);
end
