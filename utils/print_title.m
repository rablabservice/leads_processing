function print_title(title, subtitle, border_char_start, border_char_end, border_step)
    % Print a title and optional subtitle with a border before and after
    % ------------------------------------------------------------------
    arguments
        title {mustBeText}
        subtitle {mustBeText} = ''
        border_char_start {mustBeText} = '~'
        border_char_end {mustBeText} = '~'
        border_step {mustBePositive} = 4
    end

    % Figure out the border lengths
    border_len = max(length(title), length(subtitle));
    border_lens_rev = [border_len:-border_step:0];
    border_lens_fwd = flip(border_lens_rev);

    % Print the starting border
    fprintf('\n');
    if ~isempty(border_char_start)
        for ii = 1:numel(border_lens_fwd)
            fprintf('%s\n', repmat(border_char_start, 1, border_lens_fwd(ii)));
        end
    end

    % Print the title and subtitle
    fprintf('%s\n\n', title);
    if ~isempty(subtitle)
        fprintf('%s\n', subtitle);
    end

    % Print the ending border
    if ~isempty(border_char_end)
        for ii = 1:numel(border_lens_rev)
            fprintf('%s\n', repmat(border_char_end, 1, border_lens_rev(ii)));
        end
    end
    fprintf('\n');
end
