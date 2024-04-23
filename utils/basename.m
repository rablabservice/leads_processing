function f = basename(p)
    % Return the basename of a file path.
    % ------------------------------------------------------------------
    [~, n, x] = fileparts(deblank(p));
    f = append(n, x);
end
