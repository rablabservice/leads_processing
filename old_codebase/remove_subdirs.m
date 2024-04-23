function remove_subdirs(paths, prefix, verbose)
% Remove subdirectories with a specified prefix in given paths.
%
% Parameters
% ----------
%   paths : cell array of strings
%     Each string is a path in which directories with the specified
%     prefix will be searched for and deleted. Required.
%   prefix : string, optional
%     String specifying the prefix of the directories to delete.
%     Default is '' (all subdirectories are deleted).
%   verbose : logical, optional
%     If true, prints the path of each directory that is deleted.
%     Default is true.
%

% Set default values for optional parameters.
switch nargin
    case 1
        prefix = '';
        verbose = true;
    case 2
        verbose = true;
end

% Ensure paths is a cell array.
if ischar(paths) || isstring(paths)
    paths = {paths};
end

% Remove subdirectories with the specified prefix.
for ii = 1:length(paths)
    parentDir = paths{ii};
    if ~isfolder(parentDir)
        continue;
    end
    subDirs = dir(fullfile(parentDir, [prefix '*']));
    for jj = 1:length(subDirs)
        if strcmp(subDirs(jj).name, '.') || strcmp(subDirs(jj).name, '..')
            continue;
        end
        if subDirs(jj).isdir
            subDir = fullfile(parentDir, subDirs(jj).name);
            % rmdir(subDir, 's');
            if verbose
                fprintf('Deleted %s\n', subDir);
            end
        end
    end
end
end
