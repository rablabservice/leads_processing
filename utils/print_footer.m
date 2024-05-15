function print_footer(caption)
    % Print a footer caption with a border before and after
    % ------------------------------------------------------------------
    arguments
        caption {mustBeText}
    end

    % Figure out the border lengths
    border_len = length(caption);

    % Print the starting border
    fprintf('\n%s\n', repmat('_', 1, border_len));
    fprintf('%s\n', repmat('-', 1, border_len));

    % Print the caption
    fprintf('%s\n', caption);

    % Print the ending border
    if border_len < 12
        fprintf('%s\n\n\n', repmat('=', 1, border_len));
    else
        fprintf('%s\n\n\n', append('+..', repmat('=', 1, border_len - 6), '..+'));
    end
end
