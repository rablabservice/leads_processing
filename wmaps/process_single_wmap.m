function process_single_wmap(...
    run_dir, ...
    overwrite ...
)

 % Create W-map for a single entered scan
    %
    % Parameters
    % ----------
    % run_dir: cell
    %     Path to a single scan scheduled to create W-maps
    % overwrite : logical
    %     If true, overwrite existing W-map
    % ------------------------------------------------------------------
    arguments
        run_dir = {}
        overwrite logical = false
    end

    % Format paths
    run_dir = abspath(run_dir);
    code_dir = fileparts(fileparts(mfilename('fullpath')));
    ref_region_file = fullfile(code_dir, 'config', 'ref_regions.csv');
    ref_regions = readtable(ref_region_file);

    % Initialize SPM jobman and PET parameter defaults
    spm_jobman('initcfg');
    spm('defaults','PET');

    % ------------------------------------------------------------------
    % If processing is already complete and overwrite is false, then return
    if processed_wmap_exist(run_dir) && ~overwrite
        fprintf('%s processing already complete, returning output files\n', scan_tag)
        return
    end

    % Get the metadata for entered scan
    scan_tag = get_scan_tag(run_dir);
    [subj, scan_type, scan_date] = parse_scan_tag(scan_tag);

    % Get age, sex, and TIV information
    metadata_file = dir(fullfile(fileparts(code_dir), 'metadata/loni/', strcat('LEADS_PTDEMOG_*.csv')));
    metadata = readtable(fullfile(metadata_file(1).folder,metadata_file(1).name),'VariableNamingRule','preserve');

    if any(strcmpi(metadata.subject_code, subj))
        dob = metadata.ptdob(strcmpi(metadata.subject_code,subj));
        age=round(years(datetime(scan_date)-dob),1); % age in years

        sex = metadata.ptgender(strcmpi(metadata.subject_code,subj));
        sex(sex==2)=0; %convert to the common convention (0 for female, 1 for male)
    else
        fprintf('Subject %s is not found in the LEADS_PTDEMOG file\n', subj)
        return
    end
    
    if scan_type == "MRI-T1"
        tiv_file = fullfile(run_dir, strcat(scan_tag,'_nu_seg8_TIV.csv'));
        if isfile(tiv_file)
            tiv = readtable(tiv_file,'VariableNamingRule','preserve','delimiter','comma').TIV_in_mL;
        else
            fprintf('%sRequired TIV value is not found for: \n', scan_tag)
            return
        end
    end


    % ------------------------------------------------------------------
    % Start the log file
    fid = log_start(run_dir);

    try
        % --------------------------------------------------------------
        % Print the module header
        log_append(fid, '', 0, 0);
        log_append(fid, 'START W-MAP PROCESSING MODULE', 0, 0);
        log_append(fid, '---------------------------', 0, 0);
        log_append(fid, '', 0, 0);

        % Print path to the current code file
        log_append(fid, 'Code file:', 0, 0);
        log_append(fid, sprintf('%s.m\n', mfilename('fullpath')), 0, 0);

        % Print the input parameters
        log_append(fid, 'Input parameters:', 0, 0);
        log_append(fid, sprintf('run_dir = %s', run_dir), 0, 0);
        log_append(fid, sprintf('overwrite = %d', overwrite), 0, 0);
        log_append(fid, sprintf('age = %d', age), 0, 0);
        log_append(fid, sprintf('sex = %d', sex), 0, 0);
        if scan_type == "MRI-T1"
            log_append(fid, sprintf('TIV = %d', tiv), 0, 0);
        end
        log_append(fid, '', 0, 0);

        % --------------------------------------------------------------
        % Define the equation for the W-map
        log_append(fid, '- Defining equation for W-map.');
        f='(i1-(i2+i3.*';
	    f=strcat(f,num2str(age));
	    f=strcat(f,'+i4.*');
        f=strcat(f, num2str(sex));
        if scan_type == "MRI-T1"
            f=strcat(f,'+i5.*');
            f=strcat(f,num2str(tiv));
            f = strcat(f,'))./i6');
        else
            f = strcat(f,'))./i5');
        end
        log_append(fid, sprintf('-- %s', f));

        % --------------------------------------------------------------
        % Define inpiut and output files
        
        if scan_type == "MRI-T1"
            w_files = cell(6,1);
            mwc1_file = abspath(fullfile(run_dir, strcat('mwc1',scan_tag, '_nu.nii')));
            if isfile(mwc1_file)
                log_append(fid, sprintf('- Smoothing mwc1 file to 8mm: %s', mwc1_file));
                % Smooth mwc1 file to 8mm
                clear matlabbatch;
                matlabbatch{1}.spm.spatial.smooth.data = {mwc1_file};
                matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
                matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                matlabbatch{1}.spm.spatial.smooth.im = 0;
                matlabbatch{1}.spm.spatial.smooth.prefix = 's8';
                spm_jobman('run',matlabbatch);

                w_files{1} = add_presuf(mwc1_file,'s8');
                w_files{5} = fullfile(code_dir, 'wmaps', scan_type, strcat(scan_type,'_beta_0004.nii'));
                w_files{6} = fullfile(code_dir, 'wmaps', scan_type, strcat(scan_type,'_sd.nii'));
                
                wmap_filename=cellstr(strcat('W-map_s8mwc1',scan_tag,'_nu.nii'));
            else
                log_append(fid, '- Cannot complete W-map processing until mwc1 file exists');
                return
            end
        else
            w_files = cell(5,1);
            pet_outfiles = get_processed_pet_files(run_dir);
            ref_region = ref_regions(strcmp(ref_regions.tracer, scan_type), :);
            ref_region = ref_region.ref_region(1); % cross-sectional ref region selection

            if isempty(pet_outfiles.("wsuvr_" + ref_region)) 
                log_append(fid, '- Cannot complete W-map processing until warped PET files exist');
                return
            end

            wmap_filename=cellstr(strcat('W-map_wr',scan_tag,'_suvr-',ref_region,'.nii'));
            w_files{1} = pet_outfiles.("wsuvr_" + ref_region);
            log_append(fid, sprintf(' - Warped PET used for W-map: %s', w_files{1}));

            w_files{5} = fullfile(code_dir, 'wmaps', scan_type, strcat(scan_type,'_sd.nii'));
        end

        w_files{2} = fullfile(code_dir, 'wmaps', scan_type, strcat(scan_type,'_beta_0001.nii'));
        w_files{3} = fullfile(code_dir, 'wmaps', scan_type, strcat(scan_type,'_beta_0002.nii'));
        w_files{4} = fullfile(code_dir, 'wmaps', scan_type, strcat(scan_type,'_beta_0003.nii'));
        log_append(fid, '- Using this files for W-map processing:');
        cellfun(@(x) log_append(fid, sprintf('    - %s', x)), w_files);

        % --------------------------------------------------------------
        % Create the W-map
        log_append(fid, '- Creating W-map');
        clear matlabbatch;
        matlabbatch{1}.spm.util.imcalc.input = w_files;
        matlabbatch{1}.spm.util.imcalc.output = char(wmap_filename);
        matlabbatch{1}.spm.util.imcalc.outdir = cellstr(run_dir);
        matlabbatch{1}.spm.util.imcalc.expression = f;
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = spm_type('float32');
        spm_jobman('run',matlabbatch); 
       
        log_append(fid, sprintf('- W-map saved as: %s', char(wmap_filename)));

        % Print the module footer
        log_append(fid, '', 0, 0);
        log_append(fid, '-------------------------', 0, 0);
        log_append(fid, 'END W-MAP PROCESSING MODULE', 0, 0);
        log_append(fid, '', 0, 0);

        % --------------------------------------------------------------
        % Close the log file
        log_close(fid);
    catch ME
        % Print the error message
        log_append(fid, '!! ERROR !!');
        log_append(fid, getReport(ME, 'extended', 'hyperlinks', 'off'), 0, 0);
 
        % Close the log file
        log_append(fid, '\nClosing log file early due to error', 0, 0);
        log_close(fid);
 
        % Rethrow the error
        rethrow(ME);
    end
end