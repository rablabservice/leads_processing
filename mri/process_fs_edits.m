function process_fs_edits( ...
    mri_dirs, ...
    scans_to_edit_dir, ...
    segment_brainstem, ...
    process_post_freesurfer, ...
    run_qc)
    % Process all MRIs that have been edited in FreeSurfer
    %
    % Parameters
    % ----------
    % mri_dirs : cell
    %     List of MRI directories to process. This variable overrides
    %     the default behavior of selecting scans to process based on
    %     the latest raw_MRI-T1_index file in scans_to_process_dir
    % scans_to_edit_dir : char or str
    %     The directory that stores edited MRI FS files
    % segment_brainstem : logical
    %     If true, segment the brainstem using segmentBS.sh
    % process_post_freesurfer : logical
    %     If true, run post-FreeSurfer processing
    % run_qc : logical, optional
    %     If true, create QC image and add new QC eval file. Default is
    %     true
    % ------------------------------------------------------------------
    arguments
        mri_dirs = {}
        scans_to_edit_dir {mustBeText} = '/mnt/coredata/processing/leads/data/freesurfer_edits'
        segment_brainstem logical = true
        process_post_freesurfer logical = true
        run_qc logical = true
    end
   
    fs_edited = true;
    overwrite=false;
    process_freesurfer=true;
    
    % Format paths
    scans_to_edit_dir = abspath(scans_to_edit_dir);

    % Schedule scans to process
    fprintf('Creating list of scans to process...\n');
    if isempty(mri_dirs)
        % Find LDS folders in the scans_to_edit_dir
        mris = dir(fullfile(scans_to_edit_dir, 'LDS*/MRI-T1*'));
        mri_dirs = fullfile({mris.folder}, {mris.name})';
    else
        mri_dirs = abspath(cellvec(mri_dirs));
        mri_dirs = mri_dirs(cellfun(@isfolder, mri_dirs));
    end
   
    raw_mrif = '';

    % Are there are any scans to process?
    if isempty(mri_dirs)
        fprintf('No MRIs to process\n');
        return;
    elseif length(mri_dirs) == 1
        fprintf('Preparing to process 1 MRI\n');
        cellfun(@(x) fprintf('  %s\n', x), mri_dirs);
    else
        fprintf('Preparing to process %d MRIs\n', length(mri_dirs));
        cellfun(@(x) fprintf('  %s\n', x), mri_dirs);
    end

   % Run processing for each MRI
    for i=1:length(mri_dirs)
        mri_dir = strrep(mri_dirs{i},'freesurfer_edits','processed');
        
        % Remove the processed MRI files
        outfiles = get_processed_mri_files(mri_dir);
        fprintf('Removing processed MRI files in %s\n', mri_dir);
        fields = fieldnames(outfiles); % Get field names
        for ii = 1:numel(fields)
            filename = outfiles.(fields{ii}); % Get filename
            if isfile(filename) % Check if the file exists
                delete(filename); % Delete the file
            end
        end

        % Find any linked PET directories
        linked_pet_dirs = pet_linked_to_mri(mri_dir);
        
        % Remove the processed PET files
        for ii=1:length(linked_pet_dirs)
            outfiles_pet = get_processed_pet_files(linked_pet_dirs{ii});
            fprintf('Removing processed PET files in %s\n', linked_pet_dirs{ii});
            fields = fieldnames(outfiles_pet); % Get field names
            
            fields(strcmp(fields, 'pet')) = []; %keep raw recentered PET
            for j = 1:numel(fields)
                filename = outfiles_pet.(fields{j}); % Get filename
                if isfile(filename) % Check if the file exists
                    delete(filename); % Delete the file
                end
            end

            if contains(linked_pet_dirs{ii},"FTP")
                if ~isempty(dir(fullfile(linked_pet_dirs{ii}, "arLDS*FTP*.pdf")))
                    filename_old=dir(fullfile(linked_pet_dirs{ii}, "arLDS*FTP*multislice.pdf")); 
                    filename_new=fullfile(filename_old.folder,...
                        [filename_old.name(1:end-4),'_',filename_old.date(1:11),'.pdf']);
                    movefile(fullfile(filename_old.folder,filename_old.name),filename_new);
                end
            end
        end

        % Move the FS edits to the processed directory
        fprintf('Move FS edits in %s\n', mri_dir);
        movefile(fullfile(mri_dirs{i},'freesurfer_7p1'), fullfile(mri_dirs{i},'freesurfer_7p1_edited'));
        movefile(fullfile(mri_dirs{i},'freesurfer_7p1_edited'), mri_dir);
        rmdir(mri_dirs{i}, 's');

        % Process MRI with new FS segmentation
        process_single_mri( ...
            mri_dir, ...
            overwrite, ...
            raw_mrif, ...
            segment_brainstem, ...
            process_freesurfer, ...
            process_post_freesurfer, ...
            run_qc, ...
            fs_edited...
        );
    end  