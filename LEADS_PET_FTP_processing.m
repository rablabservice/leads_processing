%% All the new FTPs have been now moved and renamed, so it is time to process them

cd (path_processed);

avftps=dir('*/*/FTP*');
avftps=strcat({avftps.folder}','/',{avftps.name}');

for i=1:size(avftps,1)

    srcftpsuvr=dir(strcat(avftps{i,1},'/LDS*suvr_infcblg.nii'));

    if size(srcftpsuvr,1)==0

        srcftpd=dir(strcat(avftps{i,1},'/LDS*FTP*.nii'));
        ftpscan=strcat(srcftpd.folder,'/',srcftpd.name);
        rftpscan=strcat(srcftpd.folder,'/r',srcftpd.name);

        srcnu=dir(strcat(avftps{i,1},'/LDS*nu.nii')); nuscan=strcat(srcnu.folder,'/',srcnu.name);
        srcaparc=dir(strcat(avftps{i,1},'/LDS*raparc+aseg.nii')); aparcscan=strcat(srcaparc.folder,'/',srcaparc.name);

        pathftp=avftps{i,1};

        run LEADS_PET_FTP_pipeline.m

        clear srcftpd ftpscan rftpscan srcnu srcaparc nuscan aparcscan pathftp

    end % end if condition there is no SUVR file in the folder, this means the scan was not processed

    clear srcftpsuvr

end

cd (path_processed);
