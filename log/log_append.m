function log_append(fid, msg, add_timestamp, indent)
    % Append to an open log file
    %
    % Parameters
    % ----------
    % fid : int
    %     File identifier for logging (default is 1 for stdout)
    % msg : str
    %     Message to append to the log file
    % add_timestamp : logical, optional
    %     If true, add a timestamp to the log message (default is true)
    % indent : int, optional
    %     Number of spaces to indent the message (default is 2)
    % ------------------------------------------------------------------
    arguments
        fid {mustBeNumeric} = 1
        msg {mustBeText} = ''
        add_timestamp logical = true
        indent {mustBeNumeric} = 2
    end

    % Format column spacing
    sep = repmat(' ', 1, indent);

    % Format the message
    msg = sprintf('%s\n', msg);

    % Get the current time in HH:mm:ss format
    if add_timestamp
        time_now = datetime('now', 'TimeZone', 'America/Los_Angeles');
        time_str = char(time_now, 'HH:mm:ss');
    else
        time_str = '';
    end

    % Write the timestamp and msg to the log file
    fprintf(fid, append(time_str, sep, msg));
end
