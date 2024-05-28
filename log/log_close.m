function log_close(fid)
    % Close log file
    % ------------------------------------------------------------------
    arguments
        fid {mustBeNumeric}
    end

    % Record time elapsed since the last call to tic
    elapsed = toc;
    hh = floor(elapsed / 3600);
    mm = floor(mod(elapsed, 3600) / 60);
    ss = floor(mod(elapsed, 60));

    % Write the last log message
    add_timestamp = false;
    indent = 0;
    if elapsed > 3600
        log_append(fid, sprintf('Elapsed time: %dh, %dm, %ds', hh, mm, ss), add_timestamp, indent);
    elseif elapsed > 60
        log_append(fid, sprintf('Elapsed time: %dm, %ds', mm, ss), add_timestamp, indent);
    else
        log_append(fid, sprintf('Elapsed time: %.3fs', elapsed), add_timestamp, indent);
    end
    log_append(fid, repmat('-', 1, 88), add_timestamp, indent);
    log_append(fid, '', 0, 0);

    % Close the log file
    logf = fopen(fid);
    fclose(fid);
    fprintf('\n%s:\n', basename(logf));
    cmd = sprintf('cat "%s"', logf);
    system(cmd);
end
