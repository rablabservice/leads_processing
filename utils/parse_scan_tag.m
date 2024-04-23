function [subj, scan_type, scan_date] = parse_scan_tag(scan_tag)
    % Parse the scan tag and return subject ID, scan type, and scan date
    % ------------------------------------------------------------------
    parts = strsplit(char(scan_tag), '_');
    subj = parts{1};
    scan_type = parts{2};
    scan_date = parts{3};
end
