% List of open inputs
% Segment: Volumes - cfg_files
% Deformations: Deformation Field - cfg_files
% Deformations: Apply to - cfg_files
% Deformations: Field of View - cfg_choice
nrun = X; % enter the number of runs here
jobfile = {'/mnt/coredata/processing/leads/code/batch/spm_batch_segment_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(4, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Volumes - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Deformations: Deformation Field - cfg_files
    inputs{3, crun} = MATLAB_CODE_TO_FILL_INPUT; % Deformations: Apply to - cfg_files
    inputs{4, crun} = MATLAB_CODE_TO_FILL_INPUT; % Deformations: Field of View - cfg_choice
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
