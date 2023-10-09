% Function to quickly generate Multislice type 2: Overlay: Image on custom background (e.g. PET on MRI in native space)
% Leonardo.Iaccarino@ucsf.edu July 9th 2021 %
% Updated December 2021 to allow automatic definition of range in images 

% Functions expects: 
% 1) Background image
% 2) Overlay Image
% 3) Plane
% 4) Alpha (Transparency on background)
% 5) Colorscale (possible choices: 'nih.lut','flow.lut','gray','inv_gray','hot','winter','bone','actc'
% 6) Range colorscale
% 7) Slices
% 8) Custom Info
% 9) Font size: Header
% 10) Font size: Custom Info

% Usage Example quickmultislice2('/my/path/to/the/file/mymri.nii','/my/path/to/the/file/mypet.nii','axial','0.7','nih.lut','0.5 3.5','-30 6 58','Centiloids=56','14','10') 

function [slovname] = quickmultislice2(bgimage,ovimage,plane,alpha,colorscale,rangecolorscale,slices,custlabel,fntsizeheader,fntsizecustom)

if nargin<2
  error('You must provide at least the background and the overlay images');
elseif nargin<3
  plane='axial';
  alpha='0.7';
  colorscale='nih.lut';
  rangecolorscale='0.5 2.5';
  slices='-30 6 58';
  custlabel='';
  fntsizeheader='14';
  fntsizecustom='14';
  fprintf(1,'****\nOnly %s argument found in the function call.\nDefaulting for remaining options\n****\n',num2str(nargin));   
elseif nargin<4
  alpha='0.7';
  colorscale='nih.lut';
  rangecolorscale='0.5 2.5';
  slices='-30 6 58';
  custlabel='';
  fntsizeheader='14';
  fntsizecustom='14';
  fprintf(1,'****\nOnly %s argument found in the function call.\nDefaulting for remaining options\n****\n',num2str(nargin));   
elseif nargin<5
  colorscale='nih.lut';
  rangecolorscale='0.5 2.5';
  slices='-30 6 58';
  custlabel='';
  fntsizeheader='14';
  fntsizecustom='14';
  fprintf(1,'****\nOnly %s argument found in the function call.\nDefaulting for remaining options\n****\n',num2str(nargin));   
elseif nargin<6
  rangecolorscale='0.5 2.5';
  slices='-30 6 58';
  custlabel='';
  fntsizeheader='14';
  fntsizecustom='14';
  fprintf(1,'****\nOnly %s argument found in the function call.\nDefaulting for remaining options\n****\n',num2str(nargin));   
elseif nargin<7
  slices='-30 6 58';
  custlabel='';
  fntsizeheader='14';
  fntsizecustom='14';
  fprintf(1,'****\nOnly %s argument found in the function call.\nDefaulting for remaining options\n****\n',num2str(nargin));   
elseif nargin<8
  custlabel='';
  fntsizeheader='14';
  fntsizecustom='14';
  fprintf(1,'****\nOnly %s argument found in the function call.\nDefaulting for remaining options\n****\n',num2str(nargin));   
elseif nargin<9
  fntsizeheader='14';
  fntsizecustom='14';
  fprintf(1,'****\nOnly %s argument found in the function call.\nDefaulting for remaining options\n****\n',num2str(nargin));   
end

if isequal(colorscale,'inv_gray')

colorscale='flipud(gray)';

end

if isequal(rangecolorscale,'auto')
tempimgvals=spm_read_vols(spm_vol(ovimage));
tempimgvals=nonzeros(tempimgvals);
tempimgvals=tempimgvals(~isnan(tempimgvals));
rangecolorscale=[min(tempimgvals) max(tempimgvals)];
else
rangecolorscale=str2num(rangecolorscale);
end

slices=str2num(slices);
fntsizeheader=str2num(fntsizeheader);
fntsizecustom=str2num(fntsizecustom);
alpha=str2num(alpha);

[pp,ff,~]=spm_fileparts(ovimage);
[~,ff2,~]=spm_fileparts(bgimage);
slovname=char(strcat(pp,'/Multislice_',cellstr(plane),'_',ff,'_on_',ff2,'.jpg'));

o = slover;
o.cbar = 2;

o.img(1).vol=spm_vol(bgimage);
o.img(1).type='structural';
o.img(1).prop=1;
o.img(2).vol = spm_vol(ovimage);
o.img(2).type = 'truecolour';
o.img(2).cmap = char(colorscale);
o.img(2).range = [rangecolorscale(1) rangecolorscale(2)];
o.img(2).prop=alpha;

o.transform = plane;
o.figure = spm_figure('GetWin','Graphics');
o = fill_defaults (o);
o.slices = slices(1):slices(2):slices(3);
o = paint(o);
jpeglab=strcat('Image: ',ff);            
hTitAx = axes('Parent',o.figure,'Position',[0 0.98 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',fntsizeheader);

jpeglab2=strcat('Background: ',ff2);            
hTitAx = axes('Parent',o.figure,'Position',[0 0.96 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab2,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',fntsizeheader);

if isempty(custlabel)==0

jpeglab3=strcat('Custom Info -->',{' '},custlabel);            
hTitAx = axes('Parent',o.figure,'Position',[0 0.94 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab3,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',fntsizecustom);

end

print(slovname,'-djpeg','-r300');

end

