function [fid, logf] = log_start(scan_dir, logf)
    % Start a new log file
    % ------------------------------------------------------------------
    arguments
        scan_dir {mustBeText} = ''
        logf {mustBeText} = ''
    end

    % Start the program timer
    tic;

    % Get the full path to the log file
    time_now = datetime('now', 'TimeZone', 'America/Los_Angeles');
    if isempty(logf)
        mustBeFolder(scan_dir);
        scan_tag = get_scan_tag(scan_dir);
        log_dir = fullfile(scan_dir, 'log');
        if ~isfolder(log_dir)
            mkdir(log_dir);
        end
        datetime_str = char(time_now, 'yyyy-MM-dd-HH-mm-ss');
        logf = fullfile(log_dir, append(scan_tag, '_processed_', datetime_str, '.log'));
    else
        log_dir = fileparts(logf);
        scan_dir = fileparts(log_dir);
    end

    % Create the log directory if it doesn't exist
    mustBeFolder(scan_dir);
    if ~isfolder(log_dir)
        mkdir(log_dir);
    end

    % Open the log file for writing
    fprintf('Starting log file: %s\n', logf);
    fid = fopen(logf, 'w');

    % Write the log file header
    date_str = char(time_now, 'yyyy-MM-dd');
    time_str = char(time_now, 'HH:mm:ss');
    user = getenv('USER');
    [~, hostname] = system('hostname');
    hostname = deblank(hostname);
    add_timestamp = false;
    indent = 0;
    if exist('scan_tag')
        log_append(fid, sprintf('SCAN: %s', scan_tag), add_timestamp, indent);
    end
    log_append(fid, sprintf('DATE: %s', date_str), add_timestamp, indent);
    log_append(fid, sprintf('TIME: %s', time_str), add_timestamp, indent);
    log_append(fid, sprintf('USER: %s', user), add_timestamp, indent);
    log_append(fid, sprintf('HOST: %s', hostname), add_timestamp, indent);
    log_append(fid, repmat('_', 1, 41), add_timestamp, indent);
    log_append(fid, append('+..', repmat('=', 1, 35), '..+'), add_timestamp, indent);
end
