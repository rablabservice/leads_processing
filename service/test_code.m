

template = spm_select(1,'any', 'Select template file (*.nv)'); %% images 
images = spm_select(Inf,'image', 'Select images to render'); %% images 
images=cellstr(images);
% images = cellfun(@(x) x(1:end-2),images,'UniformOutput',false);
opt_file=spm_select(1,'any', 'Select BrainNet Option file'); %% images

     for i=1:size(images,1) 

         
         temp_img=images{i};
         [p,f,e]=spm_fileparts(char(temp_img));
         newfname=strcat(p,'/','3DRend_',f,'.jpg');


BrainNet_MapCfg(template,temp_img,opt_file,newfname);
close all
clear temp_img newfname p f e 

     end
     
     
     
     
     
   