% List of open inputs
% Smooth: Images to smooth - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/mnt/coredata/processing/leads/code/batch/spm_batch_smooth_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Smooth: Images to smooth - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
