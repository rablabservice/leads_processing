clear matlabbatch;
matlabbatch{1}.spm.spatial.normalise.write.subj = struct('def', {}, 'resample', {});
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [Inf Inf Inf
                                                          Inf Inf Inf];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
spm_jobman('run', matlabbatch);
