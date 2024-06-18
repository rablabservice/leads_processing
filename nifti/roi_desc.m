function output = roi_desc(scanfs, maskfs, subrois, aggf)
    % Calculate descriptive statistics for scan within one or more ROIs
    %
    % Parameters
    % ----------
    % scanfs : char or cell array of char
    %   Path to the image file(s) to extract data from
    % maskfs : char or cell array of char
    %   Path to the mask file(s) to extract data from
    % subrois : struct
    %   Struct with sub-ROI values to extract data from
    % aggf : function handle or cell array of function handles
    %   Function(s) to apply to the data within each ROI
    % ------------------------------------------------------------------
    arguments
        scanfs = []
        maskfs = []
        subrois struct = struct()
        aggf function_handle = @mean
    end

    % Format parameters
    scanfs = cellvec(abspath(scanfs));
    maskfs = cellvec(abspath(maskfs));
    aggfs = struct();
    if isa(aggf, 'function_handle')
        aggfs.(functions(aggf).function) = aggf;
    elseif iscell(aggf)
        for ii = 1:numel(aggf)
            aggfs.(functions(aggf{ii}).function) = aggf{ii};
        end
    elseif isstruct(aggf)
        aggfs = aggf;
    end

    % Check that all files exist
    mustBeFile(scanfs);
    mustBeFile(maskfs);

    % Initialize the output array
    n_scanfs = numel(scanfs);
    n_maskfs = numel(maskfs);
    subroi_names = fieldnames(subrois);
    n_subrois = numel(subroi_names);
    aggf_names = fieldnames(aggfs);
    n_aggfs = numel(aggf_names);
    if n_subrois > 0
        cols = cat(1, {'image_file'; 'aparc_file'; 'roi'}, aggf_names, {'voxel_count'});
        n_rows = n_scanfs * n_subrois * n_aggfs;
    else
        cols = cat(1, {'image_file'; 'mask_file'}, aggf_names, {'voxel_count'});
        n_rows = n_scanfs * n_maskfs * n_aggfs;
    end
    n_cols = numel(cols);
    output = cell(n_rows, n_cols);

    % Loop over each scan, region, and function
    iRow = 1;
    if n_subrois > 0
        disp('Sub-ROIs not yet implemented');
        return
    else
        for iScan = 1:n_scanfs
            % Load and flatten the data array
            scanf = scanfs{iScan};
            dat = spm_read_vols(spm_vol(scanf));
            dat_shp = size(dat);
            dat = dat(:);
            dat_is_finite = isfinite(dat);
            for iMask = 1:n_maskfs
                % Load and flatten the mask array
                maskf = maskfs{iMask};
                mask = spm_read_vols(spm_vol(maskf));
                mask_shp = size(mask);
                mask = mask(:);

                % Check that dat and mask have equal shape
                assert(isequal(dat_shp, mask_shp), 'Data and mask shape must be equal');

                % Select voxel_count that are finite and within the mask
                mask_is_gt0 = mask > 0;

                % Calculate the number of voxel_count in the mask
                mask_vol = sum(mask_is_gt0);

                % Select values to aggregate over
                dat_mask = dat(dat_is_finite & mask > 0);

                % Loop over each aggregation function
                for iAggf = 1:n_aggfs
                    % Apply the aggregation function
                    aggf_name = aggf_names{iAggf};
                    aggf = aggfs.(aggf_name);
                    aggf_val = aggf(dat_mask);

                    % Append value to the output array
                    outrow = {scanf, maskf, aggf_val, mask_vol};
                    output(iRow, :) = outrow;
                    iRow = iRow + 1;
                end
            end
        end
    end

    % Convert output to table
    output = cell2table(output, 'VariableNames', cols);
end
