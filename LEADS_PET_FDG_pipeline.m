%%% script to run when in the FDG folder where an FDG needs to be processed
%%% outputs: Mask ref region pons, SUVR image
%%% inputs FDG scan, nu MRI, bs seg

%nii_setOrigin(fdgscan);
% Different module to set origin
file = deblank(fdgscan);
st.vol = spm_vol(file);
vs = st.vol.mat\eye(4);
vs(1:3,4) = (st.vol.dim+1)/2;
spm_get_space(st.vol.fname,inv(vs));

spm('defaults','PET');
clear matlabbatch;
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(nuscan);
matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(fdgscan);
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch); clear matlabbatch;

%%% Let's take into account the bs

r1n=spm_vol(bsscan); % reading the brainstem parcellation
r1=spm_read_vols(r1n);

img=spm_vol(rfdgscan); % Reading img header 
img1=spm_read_vols(img); % Reading img values 

mask = (r1 == 174); % Creating the logical mask to separate the pons ROI
img_mask=img1.*mask; % Creating the masked image
temp=nonzeros(img_mask);
temp=temp(~isnan(temp));                              
meanpons=mean(temp); %% stored the ref region value to create the SUVR image

exp=char(strcat('i1/',num2str(meanpons)));
newfname=char(strcat(fdgscan(1:end-4),'_suvr_pons.nii'));

spm('defaults','PET');
clear matlabbatch;
matlabbatch{1}.spm.util.imcalc.input = cellstr(bsscan);
matlabbatch{1}.spm.util.imcalc.output = 'pons_ref_mask';
matlabbatch{1}.spm.util.imcalc.outdir = cellstr(pathfdg);
matlabbatch{1}.spm.util.imcalc.expression = 'i1==174';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 0;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
matlabbatch{2}.spm.util.imcalc.input = cellstr(rfdgscan);
matlabbatch{2}.spm.util.imcalc.output = newfname;
matlabbatch{2}.spm.util.imcalc.outdir = {''};
matlabbatch{2}.spm.util.imcalc.expression = exp;
matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{2}.spm.util.imcalc.options.mask = 0;
matlabbatch{2}.spm.util.imcalc.options.interp = 1;
matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch); clear matlabbatch;

%%%%% code to create Multi-axial views for reporting %%%%

% 1. SUVR

[pp,ff,~]=spm_fileparts(newfname);
slovname=char(strcat(pp,'/Multiaxial_',ff,'.pdf'));
o = slover;
o.cbar = 2;
o.img(1).vol=spm_vol(nuscan);
o.img(1).type='structural';
o.img(1).prop=1;
o.img(2).vol = spm_vol(newfname);
o.img(2).type = 'truecolour';
o.img(2).cmap = 'nih.lut';
o.img(2).range = [0 2.2];
o.img(2).prop=0.7;
o.transform = 'axial';
o.figure = spm_figure('GetWin','Graphics');
o = fill_defaults (o);
o.slices = -30:6:58;
o = paint(o);
crdate=char(pp(size(pp,2)-9:end)); crid=char(ff(1:10));
jpeglab=strcat('ID:',{' '},crid,{' '},'***',{' '},'FDG-PET',{' '},'***',{' '},'LEADS.PETCORE@ucsf.edu');            
hTitAx = axes('Parent',o.figure,'Position',[0 0.97 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
jpeglab2=strcat('Scan Date:',{' '},datestr(crdate),{' '},'***',{' '},'Quantification Date:',{' '},date);
hTitAx2 = axes('Parent',o.figure,'Position',[0 0.95 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab2,'Parent',hTitAx2,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
jpeglab3=strcat('SUVR Map - Ref region: Pons',{' '},'***',{' '},'Left is Left');
hTitAx3 = axes('Parent',o.figure,'Position',[0 0.93 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab3,'Parent',hTitAx3,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
print(slovname,'-dpdf','-r300');

%%%%% code to create Multi-axial views for reporting QC refreg %%%%
slovnameref=char(strcat(pp,'/Refreg_',ff,'_ponsref.jpg'));
o = slover;
o.cbar = 2;
o.img(1).vol=spm_vol(newfname);
o.img(1).type='truecolour';
o.img(1).cmap = 'gray';
o.img(1).prop=1;
o.img(2).vol = spm_vol(char(strcat(pp,'/pons_ref_mask.nii')));
o.img(2).type = 'split';
o.img(2).cmap = char('actc');
o.img(2).range = [0.5 1.2];
o.transform = 'sagittal';
o.figure = spm_figure('GetWin','Graphics');
o = fill_defaults (o);
o.slices = -32:8:30;
o = paint(o);
jpeglab=strcat('ID:',{' '},crid,{' '},'***',{' '},'FDG-PET',{' '},'***',{' '},'LEADS.PETCORE@ucsf.edu');            
hTitAx = axes('Parent',o.figure,'Position',[0 0.97 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
jpeglab2=strcat('Scan Date:',{' '},datestr(crdate),{' '},'***',{' '},'Quantification Date:',{' '},date);
hTitAx2 = axes('Parent',o.figure,'Position',[0 0.95 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab2,'Parent',hTitAx2,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
jpeglab3=strcat('Ref region QC: Pons on FDG-PET',{' '},'***',{' '},'Left to Right');
hTitAx3 = axes('Parent',o.figure,'Position',[0 0.93 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab3,'Parent',hTitAx3,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
print(slovnameref,'-djpeg','-r300');

clear r1n r1 img img1 mask img_mask temp meanpons exp newfname slovname slovnameref


