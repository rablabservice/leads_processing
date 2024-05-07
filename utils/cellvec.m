function outdat = cellvec(indat)
    % Recast indat as an n x 1 cell array, for various indat classes
    %
    % If indat is a char, string, or cell array of only text elements,
    % return a cell array of character vectors
    % ------------------------------------------------------------------
    if iscellstr(indat)
        outdat = indat;
    elseif ischar(indat) || isstring(indat)
        outdat = cellstr(indat);
    elseif isstruct(indat)
        outdat = struct2cell(indat);
        % If all elements are char or string, convert cell to cellstr
        if all(cellfun(@(x) ischar(x) || isstring(x), outdat))
            outdat = cellstr(outdat);
        end
    elseif iscell(indat)
        if all(cellfun(@(x) ischar(x) || isstring(x), indat))
            outdat = cellstr(indat);
        else
            outdat = indat;
        end
    elseif isnumeric(indat)
        outdat = num2cell(indat);
    elseif islogical(indat)
        outdat = num2cell(double(indat));
    elseif istabular(indat)
        outdat = table2cell(indat);
    else
        error('Unsupported input type');
    end
    outdat = outdat(:);
end
