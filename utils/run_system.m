function status = run_system(cmd, fid, raise_error_on_fail, suppress_output)
    % Run a system command and check for errors
    % ------------------------------------------------------------------
    arguments
        cmd {mustBeText}
        fid {mustBeNumeric} = 1
        raise_error_on_fail logical = true
        suppress_output logical = false
    end

    % Redirect stderr to stdout
    if ~contains(cmd, '2>&1')
        cmd = append(cmd, ' 2>&1');
    end

    % Redirect stdout to log file
    if fid ~= 1
        logf = fopen(fid);
        cmd = append(cmd, ' >> ', logf);
    % Suppress output if verbose is false
    elseif suppress_output && ~contains(cmd, '> /dev/null')
        cmd = append(cmd, ' > /dev/null');
    end

    % Run system command and check for errors
    status = system(cmd);
    if status ~= 0
        err_msg = sprintf('\nError executing %s (exit code = %d)\n', cmd, status);
        if raise_error_on_fail
            error(err_msg)
        else
            fprintf(err_msg)
        end
    end
end
