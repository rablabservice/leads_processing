%%%%% code to create Multi-axial views for reporting %%%%
%%%%% Sep 2021: slight edit

vol1=tempwmapname;
vol2='/mnt/coredata/Projects/Resources/templates/mean_112CN_MRIT1.nii';

[pp,ff,ee]=spm_fileparts(vol1);
slovname=char(strcat(pp,'/Wmap_Multiaxial_',ff,'.jpg'));
    myclass = 'slover';

% Default object structure
defstruct = struct('img', [], ...
    'transform', 'axial', ...
    'slicedef', [], ...
    'slices', [], ...
    'figure', [], ...
    'figure_struct', [], ...
    'refreshf', 1, ...
    'clf', 1, ...
    'resurrectf', 1, ...
    'userdata', 1, ...
    'area', [], ...
    'xslices', [], ...
    'cbar', [], ...
    'labels', [], ...
    'callback', ';', ...
    'printstr', 'print -dpsc -painters -noui', ...
    'printfile', 'slices.ps');

others = [];
temp_rend=vol1;
temp_nu=vol2;
o = slover;

o.cbar = 2;

o.img(1).vol=spm_vol (temp_nu);
o.img(1).type='structural';
o.img(1).prop=1;
o.img(2).vol = spm_vol (temp_rend);
o.img(2).type = 'truecolour';
o.img(2).cmap = 'nih.lut';

if strcmp(tempmod,'FTP')==1
o.img(2).range = [1.65 30];
elseif strcmp(tempmod,'FBB')==1
o.img(2).range = [1.65 25];
elseif strcmp(tempmod,'FDG')==1
o.img(2).range = [-1.65 -5];
elseif strcmp(tempmod,'MRI')==1
o.img(2).range = [-1.65 -5];   
end

o.img(2).prop=0.7;

o.transform = 'axial';
o.figure = spm_figure('GetWin','Graphics');
o = fill_defaults (o);
o.slices = -30:6:58;
o = paint(o);
crdate=char(pp(size(pp,2)-9:end)); crid=tempid;
jpeglab=strcat('ID:',{' '},crid,{' '},'***',{' '},tempmod,{' '},'***',{' '},'LEADS.PETCORE@ucsf.edu');            
hTitAx = axes('Parent',o.figure,'Position',[0 0.97 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
jpeglab2=strcat('Scan Date:',{' '},datestr(crdate),{' '},'***',{' '},'Quantification Date:',{' '},date);
hTitAx2 = axes('Parent',o.figure,'Position',[0 0.95 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab2,'Parent',hTitAx2,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);

if strcmp(tempmod,'FTP')==1
jpeglab3=strcat('WMap - Age-corrected - 61CN as control group',{' '},'***',{' '},'Left is Left');
elseif strcmp(tempmod,'FBB')==1
jpeglab3=strcat('WMap - Age-corrected - 61CN as control group',{' '},'***',{' '},'Left is Left');
elseif strcmp(tempmod,'FDG')==1
jpeglab3=strcat('WMap - Age-corrected - 24CN as control group',{' '},'***',{' '},'Left is Left');
elseif strcmp(tempmod,'MRI')==1
jpeglab3=strcat('WMap - Age/TIV-corrected - 61CN as control group',{' '},'***',{' '},'Left is Left');
end

hTitAx3 = axes('Parent',o.figure,'Position',[0 0.93 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab3,'Parent',hTitAx3,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',12);
print(slovname,'-djpeg','-r300');
            