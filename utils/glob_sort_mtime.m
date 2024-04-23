function files = glob_sort_mtime(globstr)
    % Return files matching pattern in most recent modified order
    % ------------------------------------------------------------------
    d = dir(globstr);
    [~, idx] = sort([d.datenum], 'descend');
    files = {d(idx).name};
    files = fullfile({d(idx).folder}, files);
end
