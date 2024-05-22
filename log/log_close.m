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
    log_append(fid, repmat('_', 1, 41), add_timestamp, indent);
    log_append(fid, repmat('-', 1, 41), add_timestamp, indent);
    log_append(fid, sprintf('Elapsed time: %dh, %dm, %ds', hh, mm, ss));
    log_append(fid, append('+..', repmat('=', 1, 35), '..+'), add_timestamp, indent);

    % Close the log file
    fclose(fid);
end
