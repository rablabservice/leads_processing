%%%% Quality Check for LEADS reference region Data %%%%%%%
%%%%%% Checking Overlap aparcaseg and PET image %%%%%%%%%

mkdir qc_aparc
output_dir=strcat(pwd,'/qc_aparc');
output_dir_spm=cellstr(output_dir);

%%% Load images %%%

 vols = spm_select(Inf,'image', 'Select images of interest'); %% Input images 
                   vols_spm=cellstr(vols);
                   numv1= size(vols_spm,1);
                   
                   vols_names=regexp(vols_spm,'B\d{2}-\d{3}\S*(?=.nii)','match');
                   vols_names_cstr=vertcat(vols_names{:});
                   vols_names_cstr=regexprep(vols_names_cstr, '/', '_')
                   
                    
                        
 aparcs = spm_select([1 inf],'image','Select co-registered aparc.aseg files:');
                    aparc_spm=cellstr(aparcs);


%%%%%%%%%%% Then create slover for each case and save the jpeg %%%%%%%%%%%

if size(vols_spm,1)==size(aparc_spm,1)

for i =1:size(vols_spm,1) 
    
    slovname=char(strcat(output_dir,'/QC_Overlay_APARC_',vols_names_cstr{i}));
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
temp_img=vols_spm{i};
temp_aparc=aparc_spm{i};
o = slover;

o.cbar = [];

o.img(1).vol = spm_vol (temp_img);
o.img(1).type = 'structural';
o.img(1).prop=1;

o.img(2).vol = spm_vol(temp_aparc);
o.img(2).type = 'contour';
o.img(2).cmap = 'red';
o.img(2).range = [999 2035];
o.img(2).prop=1;

o.transform = 'axial';
o.cbar=[1];
o.figure = spm_figure('GetWin','Graphics');
o = fill_defaults (o);
o.slices = -42:6:72;
o = paint(o);
jpeglab=strcat(vols_names_cstr{i});
hTitAx = axes('Parent',o.figure,'Position',[0.20 0.06 0.06 0.02],'Visible','off');
text(0.5,0,jpeglab,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','red','FontSize',14);
print(slovname,'-djpeg','-r150');
    
end

msg=strcat('Success!', {' '}, num2str(size(vols_spm,1)),{' '},'JPEG files were generated. Have fun!');
disp(char(msg));

else
    
    disp('Warning! You entered a different number of images and aparcaseg files. Try again!');
    
end

clear;