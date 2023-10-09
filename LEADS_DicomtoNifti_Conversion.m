% Convert dicoms downloaded from LONI to nifti's 
% available in newdata folder and delete dicoms
% June 2023; nidhi.mundada@ucsf.edu

% find MRIs to convert 
cd (path_newmris)

dcmfiles=rdir('**/*.dcm');
dcmfiles=struct2cell(dcmfiles)';
dcmfiles=dcmfiles(:,1);
dcmfiles(:,2)=strcat(pwd,'/',cellfun(@fileparts,dcmfiles(:,1),'un',0));

mridirs=unique(dcmfiles(:,2));

for i=1:size(mridirs,1)

    cd (mridirs{i,1});
    
    dicm2nii(pwd,pwd,'3D.nii'); % this command will convert a dicom dataset into all the composing frames
    convdicoms=dir('*.nii');
    convdicoms=transpose(struct2cell(convdicoms));
    convdicoms=convdicoms(:,1); %% list of processed nii

    %rename nifti file to unique name including LEADS ID, series
    %description, date, and image ID for safety
    newname=strsplit((mridirs{i,1}),'/');
    newname=newname(end-4+1:end); % selecting last 4 folders from the filepath 
    newname=strcat(strjoin(newname,'_'),'.nii');
    movefile(char(convdicoms), newname)
    
    delete *.mat *.dcm 
    
    clear convdicoms dcmfiles newname i 

end 

fprintf(2,'%d MRI scans have been converted from dicom to nifti format. \n', size(mridirs,1));


% find FBBs to convert 
cd (path_newfbbs)

dcmfiles=rdir('**/*.dcm');
dcmfiles=struct2cell(dcmfiles)';
dcmfiles=dcmfiles(:,1);
dcmfiles(:,2)=strcat(pwd,'/',cellfun(@fileparts,dcmfiles(:,1),'un',0));

fbbdirs=unique(dcmfiles(:,2));

for i=1:size(fbbdirs,1)

    cd (fbbdirs{i,1});
    
    dicm2nii(pwd,pwd,'3D.nii'); % this command will convert a dicom dataset into all the composing frames
    convdicoms=dir('*.nii');
    convdicoms=transpose(struct2cell(convdicoms));
    convdicoms=convdicoms(:,1); %% list of processed nii

    %rename nifti file to unique name including LEADS ID, series
    %description, date, and image ID for safety
    newname=strsplit((fbbdirs{i,1}),'/');
    newname=newname(end-4+1:end); % selecting last 4 folders from the filepath 
    newname=strcat(strjoin(newname,'_'),'.nii');
    movefile(char(convdicoms), newname)
    
    delete *.mat *.dcm 
    
    clear convdicoms dcmfiles newname i 

end 

fprintf(2,'%d FBB-PET scans have been converted from dicom to nifti format. \n', size(fbbdirs,1));

% find FTPs to convert 
cd (path_newftps)

dcmfiles=rdir('**/*.dcm');
dcmfiles=struct2cell(dcmfiles)';
dcmfiles=dcmfiles(:,1);
dcmfiles(:,2)=strcat(pwd,'/',cellfun(@fileparts,dcmfiles(:,1),'un',0));

ftpdirs=unique(dcmfiles(:,2));

for i=1:size(ftpdirs,1)

    cd (ftpdirs{i,1});
    
    dicm2nii(pwd,pwd,'3D.nii'); % this command will convert a dicom dataset into all the composing frames
    convdicoms=dir('*.nii');
    convdicoms=transpose(struct2cell(convdicoms));
    convdicoms=convdicoms(:,1); %% list of processed nii
    
    %rename nifti file to unique name including LEADS ID, series
    %description, date, and image ID for safety
    newname=strsplit((ftpdirs{i,1}),'/');
    newname=newname(end-4+1:end); % selecting last 4 folders from the filepath 
    newname=strcat(strjoin(newname,'_'),'.nii');
    movefile(char(convdicoms), newname)

    delete *.mat *.dcm 
    
    clear convdicoms dcmfiles newname i 

end 

fprintf(2,'%d FTP-PET scans have been converted from dicom to nifti format. \n', size(ftpdirs,1));

% find FTPs to convert 
cd (path_newfdgs)

dcmfiles=rdir('**/*.dcm');
dcmfiles=struct2cell(dcmfiles)';
dcmfiles=dcmfiles(:,1);
dcmfiles(:,2)=strcat(pwd,'/',cellfun(@fileparts,dcmfiles(:,1),'un',0));

fdgdirs=unique(dcmfiles(:,2));

for i=1:size(fdgdirs,1)

    cd (fdgdirs{i,1});
    
    dicm2nii(pwd,pwd,'3D.nii'); % this command will convert a dicom dataset into all the composing frames
    convdicoms=dir('*.nii');
    convdicoms=transpose(struct2cell(convdicoms));
    convdicoms=convdicoms(:,1); %% list of processed nii
    
    %rename nifti file to unique name including LEADS ID, series
    %description, date, and image ID for safety
    newname=strsplit((fdgdirs{i,1}),'/');
    newname=newname(end-4+1:end); % selecting last 4 folders from the filepath 
    newname=strcat(strjoin(newname,'_'),'.nii');
    movefile(char(convdicoms), newname)

    delete *.mat *.dcm 
    
    clear convdicoms dcmfiles newname i 

end 

fprintf(2,'%d FDG-PET scans have been converted from dicom to nifti format. \n', size(fdgdirs,1));
