function myres = petexist(id,tp,mod)

path_processed='/mnt/coredata/Projects/LEADS/data_f7p1/processed/';

    % Check if PET data exists for a given subject at a given timepoint

    x = dir(strcat(path_processed,id,'/',tp,'/',mod,'*'));

    myres = size(x,1)>0;
end

