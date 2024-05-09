function print_title(title, subtitle)
    % Print a title and optional subtitle with a border before and after
    % ------------------------------------------------------------------
    arguments
        title {mustBeText}
        subtitle {mustBeText} = ''
    end

    % Figure out the border lengths
    border_len = max(length(title), length(subtitle));

    % Print the starting border
    if border_len < 12
        topbot = append(repmat('=', 1, border_len), '\n');

    else
        topbot = append('+..', repmat('=', 1, border_len - 6), '..+');
    end
    fprintf('\n%s\n', topbot);
    fprintf('%s\n', repmat('-', 1, border_len));

    % Print the title and subtitle
    if isempty(subtitle)
        fprintf('\n%s\n\n', title);
    else
        pad = floor((border_len - min(length(title), length(subtitle))) / 2);
        if length(title) > length(subtitle)
            subtitle = append(repmat(' ', 1, pad), subtitle, repmat(' ', 1, pad));
        elseif length(title) < length(subtitle)
            title = append(repmat(' ', 1, pad), title, repmat(' ', 1, pad));
        end
        fprintf('%s\n\n', title);
        fprintf('%s\n', subtitle);
    end

    % Print the ending border
    fprintf('%s\n', repmat('_', 1, border_len));
    fprintf('%s\n\n', topbot);
end
