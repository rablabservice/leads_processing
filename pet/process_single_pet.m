function outfiles = process_single_pet(pet_dir, overwrite)
% Process a single PET scan
% ----------------------------------------------------------------------
    % Print the module header
    title = 'PET PROCESSING MODULE';
    subcap = append('SCAN = ', scan_tag);
    border_len = max(length(title), length(subcap));
    border_char = '~';
    fprintf('%s\n', repmat(border_char, 1, border_len));
    fprintf('%s\n\n', repmat(border_char, 1, border_len));
    fprintf('%s\n\n', title);
    fprintf('%s\n\n', subcap);
    while border_len > 0
        fprintf('%s\n', repmat(border_char, 1, border_len));
        border_len = border_len - 4;
    end
end
