clear matlabbatch;
matlabbatch{1}.spm.tools.oldnorm.write.subj = struct('matname', {}, 'resample', {});
matlabbatch{1}.spm.tools.oldnorm.write.roptions.preserve = 0;
matlabbatch{1}.spm.tools.oldnorm.write.roptions.bb = [Inf Inf Inf
                                                      Inf Inf Inf];
matlabbatch{1}.spm.tools.oldnorm.write.roptions.vox = [1 1 1];
matlabbatch{1}.spm.tools.oldnorm.write.roptions.interp = 4;
matlabbatch{1}.spm.tools.oldnorm.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.tools.oldnorm.write.roptions.prefix = 'a';
spm_jobman('run', matlabbatch);
