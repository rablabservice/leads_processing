function multislice_file = save_fdg_multislice_pdf( ...
    pet_dir, ...
    overwrite, ...
    axial_slices, ...
    coronal_slices, ...
    clim, ...
    cmap, ...
    alpha, ...
    tmp_dir ...
)
    % Save multislice PDF of the affine-transformed FDG SUVR (pons)
    % overlaid on the affine-transformed nu MRI.
    %
    % Parameters
    % ----------
    % pet_dir : str
    %     Path to the processed PET directory
    % overwrite : bool, optional
    %     Overwrite output PDF if it already exists. Default is false.
    % axial_slices : array, optional
    %     Axial slices to display. Default is -30:6:58
    % coronal_slices : array, optional
    %     Coronal slices to display. Default is -65:8:50
    % clim : array, optional
    %     Colorscale limits for the SUVR image. Default is [0.1 2.2]
    % cmap : str, optional
    %     Colormap for the SUVR image. Default is 'nih.lut'
    % alpha : float, optional
    %     Overlay transparency of the SUVR image on MRI. Default is 0.7
    % tmp_dir : str, optional
    %     Path to the temporary directory. Default is '/mnt/tmp-scratch'
    %
    % Output
    % ------
    % multislice_file : str
    %     Path to the output multislice PDF file
    % ------------------------------------------------------------------
    arguments
        pet_dir {mustBeFolder}
        overwrite logical = false
        axial_slices = -30:6:58
        coronal_slices = -65:8:50
        clim = [0.1 2.2]
        cmap = 'nih.lut'
        alpha = 0.7
        tmp_dir = '/mnt/tmp-scratch'
    end

    % Get scan info
    pet_dir = abspath(pet_dir);
    scan_tag = get_scan_tag(pet_dir);
    [subj, scan_type, scan_date] = parse_scan_tag(scan_tag);

    % Get the input files
    pet_files = get_processed_pet_files(pet_dir);
    mri_files = get_processed_mri_files(fullfile(pet_dir, 'mri'));
    asuvr_file = pet_files.asuvr_pons;
    anu_file = mri_files.anu;
    mustBeFile(asuvr_file);
    mustBeFile(anu_file);

    % Check if the output file exists
    multislice_file = strrep(add_presuf(asuvr_file, '', '_multislice'), '.nii', '.pdf');
    if isfile(multislice_file) && ~overwrite
        return
    end

    % Set parameters for the multislice figure
    date_saved = char(datetime('today', 'Format', 'yyyy-MM-dd'));
    labels = { ...
        sprintf('Subject ID: %s | Scan Date: %s | Tracer: %s', subj, scan_date, scan_type), ...
        'SUVR Image (ref-region: pons) in neurological view (left is left)', ...
        sprintf('File saved: %s | Email: LEADS.PETCORE@ucsf.edu', date_saved), ...
        'PET now shown at 6mm FWHM. Comparison with older, 8mm data not advised' ...
    };
    fs = 10;
    lab_colors = {'black', 'black', 'black', 'red'};
    y_start = 0.985;
    y_step = 0.02;

    % ------------------------------------------------------------------
    % Create the axial multislice figure
    axial_slices_file = fullfile( ...
        tmp_dir, add_presuf(basename(multislice_file), '', '_axial') ...
    );

    o = slover;
    o.cbar = 2;
    o.img(1).vol = spm_vol(anu_file);
    o.img(1).type = 'structural';
    o.img(1).prop = 1;
    o.img(2).vol = spm_vol(asuvr_file);
    o.img(2).type = 'truecolour';
    o.img(2).cmap = cmap;
    o.img(2).range = clim;
    o.img(2).prop = alpha;
    o.transform = 'axial';
    o.figure = spm_figure('GetWin', 'Graphics');
    o = fill_defaults(o);
    o.slices = axial_slices;
    o = paint(o);
    o = add_figure_labels(o, labels, lab_colors, y_start, y_step, fs);

    % Save the figure
    print(axial_slices_file, '-dpdf', '-r300');

    % ------------------------------------------------------------------
    % Create the coronal multislice figure
    coronal_slices_file = fullfile( ...
        tmp_dir, add_presuf(basename(multislice_file), '', '_coronal') ...
    );

    o = slover;
    o.cbar = 2;
    o.img(1).vol = spm_vol(anu_file);
    o.img(1).type = 'structural';
    o.img(1).prop = 1;
    o.img(2).vol = spm_vol(asuvr_file);
    o.img(2).type = 'truecolour';
    o.img(2).cmap = cmap;
    o.img(2).range = clim;
    o.img(2).prop = alpha;
    o.transform = 'coronal';
    o.figure = spm_figure('GetWin', 'Graphics');
    o = fill_defaults(o);
    o.slices = coronal_slices;
    o = paint(o);
    o = add_figure_labels(o, labels, lab_colors, y_start, y_step, fs);

    % Save the figure
    print(coronal_slices_file, '-dpdf', '-r300');

    % ------------------------------------------------------------------
    % Merge the axial and coronal multislice figures
    cmd = sprintf( ...
        'pdfunite %s %s %s', axial_slices_file, coronal_slices_file, multislice_file ...
    );
    system(cmd);
    fprintf('Saved %s\n', multislice_file);

    % Delete the temporary files
    delete axial_slices_file coronal_slices_file;
end


function o = add_figure_labels(o, labels, lab_colors, y_start, y_step, fs)
    % Add labels to the multislice figure.
    %
    % Parameters
    % ----------
    % o : slover
    %     slover object
    % labels : cell
    %     Cell array of strings to add as labels
    % lab_colors : cell
    %     Cell array of colors for each label
    % y_start : float
    %     Starting y position for the labels
    % y_step : float
    %     Step size for each label
    % fs : int
    %     Font size for the labels
    % ------------------------------------------------------------------
    for ii = 1:length(labels)
        y_pos = y_start - ((ii - 1) * y_step);
        ax = axes( ...
            'Parent', o.figure, 'Position', [0 y_pos 0.06 0.02], 'Visible', 'off' ...
        );
        text( ...
            0.5, 0, labels{ii}, 'Parent', ax, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'baseline', 'Color', lab_colors{ii}, 'FontSize', fs ...
        );
    end
end
