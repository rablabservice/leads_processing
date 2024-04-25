function [status, cmdout] = run_system_cmd(cmd, raise_error_on_fail, verbose)
    % Run a system command and check for errors
    % ------------------------------------------------------------------
    arguments
        cmd {mustBeText}
        raise_error_on_fail logical = true
        verbose logical = true
    end

    % Redirect stderr to stdout
    if ~contains(cmd, '2>&1')
        cmd = cmd + " 2>&1";
    end
    if ~verbose && ~contains(cmd, '> /dev/null')
        cmd = cmd + " > /dev/null";
    end

    % Run system command and check for errors
    [status, cmdout] = system(cmd);
    if status ~= 0
        err_msg = sprintf('\nError executing %s (exit code = %d)\n%s\n', cmd, status, cmdout);
        if raise_error_on_fail
            error(err_msg)
        else
            fprintf(err_msg)
        end
    else
        if verbose
            fprintf('%s\n', cmdout)
        end
    end
end
