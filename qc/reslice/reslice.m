matlabbatch{1}.spm.spatial.coreg.write.ref = {'<RTPM_PATH>,1'};
matlabbatch{1}.spm.spatial.coreg.write.source = {
                                                 '<DATA_PATH>,1'                                                 
                                                 };
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'qc';

spm_jobman('run', matlabbatch);
