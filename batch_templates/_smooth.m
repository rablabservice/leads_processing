clear matlabbatch;
matlabbatch{7}.spm.spatial.smooth.data = '<UNDEFINED>';
matlabbatch{7}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{7}.spm.spatial.smooth.dtype = 0;
matlabbatch{7}.spm.spatial.smooth.im = 0;
matlabbatch{7}.spm.spatial.smooth.prefix = 's';
spm_jobman('run', matlabbatch);
