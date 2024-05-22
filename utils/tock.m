function elapsed = tock()
    % Print how long the global timer has been running
    % ------------------------------------------------------------------
    elapsed = toc;
    if elapsed < 60
        fprintf('Elapsed time: %fs\n', elapsed);
        return;
    end

    hh = floor(elapsed / 3600);
    mm = floor(mod(elapsed, 3600) / 60);
    ss = floor(mod(elapsed, 60));

    % Print the elapsed time
    if hh > 0
        fprintf('Elapsed time: %dh, %dm, %ds\n', hh, mm, ss);
    else
        fprintf('Elapsed time: %dm, %ds\n', mm, ss);
    end
end
