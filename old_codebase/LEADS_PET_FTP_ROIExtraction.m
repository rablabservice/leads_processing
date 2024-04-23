% Create a temp table with the info from the "updated" list of scans
% available in the processed folder
% Sep 2021: Code cleanup

% Feedback for the user

fprintf(1,'Checking what FTP scans are new...\n');

% Store IDs for regions that will be used later

ind_brk12_nohcp=[1006;2006];
ind_brk12=[1006;2006;17;53];
ind_brk34=[1016;1007;1013;18;2016;2007;2013;54;1015;1002;1026;1023;1010;1035;1009;1033;2015;2002;2026;2023;2010;2035;2009;2033];
ind_brk56=[1028;1012;1014;1032;1003;1027;1018;1019;1020;1011;1031;1008;1030;1029;1025;1001;1034;2028;2012;2014;2032;2003;2027;2018;2019;2020;2011;2031;2008;2030;2029;2025;2001;2034;1021;1022;1005;1024;1017;2021;2022;2005;2024;2017];
ind_metaroi=[1006;2006;18;54;1016;2016;1007;2007;1009;2009;1015;2015];

% Generate the entire list of images now available at all timepoints in the
% processed folder

newldsids=regexp(tdb_ftp(:,1),'LDS\d{7}','match','once');
newldsdates=regexp(tdb_ftp(:,1),'\d{4}-\d{2}-\d{2}','match','once');
newinfo=horzcat(newldsids,newldsdates);
newinfo=array2table(newinfo);
newinfo.Properties.VariableNames={'ID','FTPPET_Date'};

%% Let's save the tdb_ftp list also for the longitudinal ref processing
% we will use this to replicate the approach of reading last database and
% add on top also for the long database

tdb_ftp_longFTP=tdb_ftp;

%% Read the last extraction database, this will be used to know what's new %%
%% This command reads all the files in the extraction folder and then reads the latest one %%
%% The final aim is to cleanup our list to only extract from new scans. If no extraction files are found then it will extract from all the scans

olddb=dir('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/FTP_ROI_Extraction*');

    if size(olddb,1)>0

    olddb=struct2cell(olddb)';
    olddb=olddb(:,[1 3]);
    olddb=array2table(olddb);
    olddb.Properties.VariableNames(1) = cellstr(strcat('Filename'));
    olddb.Properties.VariableNames(2) = cellstr(strcat('DateCreated'));
    olddb.DateCreated=datetime(olddb.DateCreated);
    filt=max(olddb.DateCreated); % Store most recent date
    olddb=olddb(olddb.DateCreated==filt,:); % Subset to select the most recent file
    oldinfo=readtable(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/',char(table2cell(olddb(1,1))))); %% ready to check differences now

    newinfo=sortrows(newinfo);
    oldinfo=sortrows(oldinfo);

    checkinfo=isequal(newinfo(:,{'ID','FTPPET_Date'}),oldinfo(:,{'ID','FTPPET_Date'}));% Check if the b-ids columns in the new and old tables is unchanged

        if checkinfo==0

        Indexa=not(ismember(newinfo(:,{'ID','FTPPET_Date'}),oldinfo(:,{'ID','FTPPET_Date'}))); %% logical indexing for the new cases. Both are based on BID+date to include longitudinal scans
        Indexb=not(ismember(oldinfo(:,{'ID','FTPPET_Date'}),newinfo(:,{'ID','FTPPET_Date'}))); %% logical indexing for the removed cases

        newcases=newinfo(Indexa,:); %% Creates a table with the New FTP scans
        remcases=oldinfo(Indexb,:); %% Creates a table with eventually removed FTP scans

            if size(newcases,1)>0

                fprintf(2,'**New FTP-PET scans were found for Cross-Sectional Processing:\n\n');
                disp(newcases);
                tdb_ftp=tdb_ftp(Indexa,:);

            end % end if condition there are new scans available

            if size(remcases,1)>0

                fprintf(2,'Warning! %s cases were in the last FTP database and now are not found:\n\n', num2str(size(remcases,1)));

                choice = input(['\n\n Choose an action:' ,...
                    '\n     [1] Delete removed cases from present database',...
                    '\n     [2] Keep removed cases, I will check later'
                    '\n     --> ']);

                if choice==1

                    filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/Removed_Cases_FTP_ROI_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
                    writetable(remcases,filename,'WriteRowNames',true)
                    oldinfo(Indexb,:)=[];

                elseif choice==2

                        fprintf(1,'Fine! I will still save a csv file to keep track of that.\n\n');
                        filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/RemovedbutKept_Cases_FTP_ROI_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
                        writetable(remcases,filename,'WriteRowNames',true)
                end % end if condition choice on removed cases

                clear choice filename

            end % end if condition existence of removed cases

            elseif checkinfo==1

                fprintf(2,'No new FTP scans needed extraction!\n');
                tdb_ftp={};

         end % end if condition there is a difference between old extraction and new extraction list

    end % end if condition there is a previous FTP extraction database available

% To be safe, let's separately scan also the longitudinal database

olddb_longFTP=dir('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/FTP_ErodedWM_Extraction*');

    if size(olddb_longFTP,1)>0

    olddb_longFTP=struct2cell(olddb_longFTP)';
    olddb_longFTP=olddb_longFTP(:,[1 3]);
    olddb_longFTP=array2table(olddb_longFTP);
    olddb_longFTP.Properties.VariableNames(1) = cellstr(strcat('Filename'));
    olddb_longFTP.Properties.VariableNames(2) = cellstr(strcat('DateCreated'));
    olddb_longFTP.DateCreated=datetime(olddb_longFTP.DateCreated);
    filt=max(olddb_longFTP.DateCreated); % Store most recent date
    olddb_longFTP=olddb_longFTP(olddb_longFTP.DateCreated==filt,:); % Subset to select the most recent file
    oldinfo_longFTP=readtable(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/',char(table2cell(olddb_longFTP(1,1))))); %% ready to check differences now

    newinfo=sortrows(newinfo);
    oldinfo_longFTP=sortrows(oldinfo_longFTP);

    checkinfo=isequal(newinfo(:,{'ID','FTPPET_Date'}),oldinfo_longFTP(:,{'ID','FTPPET_Date'}));% Check if the b-ids columns in the new and old tables is unchanged

        if checkinfo==0

        Indexa=not(ismember(newinfo(:,{'ID','FTPPET_Date'}),oldinfo_longFTP(:,{'ID','FTPPET_Date'}))); %% logical indexing for the new cases. Both are based on BID+date to include longitudinal scans
        Indexb=not(ismember(oldinfo_longFTP(:,{'ID','FTPPET_Date'}),newinfo(:,{'ID','FTPPET_Date'}))); %% logical indexing for the removed cases

        newcases_longFTP=newinfo(Indexa,:); %% Creates a table with the New FTP scans
        remcases_longFTP=oldinfo_longFTP(Indexb,:); %% Creates a table with eventually removed FTP scans

            if size(newcases_longFTP,1)>0

                fprintf(2,'**New FTP-PET scans were found for Longitudinal Processing:\n\n');
                disp(newcases_longFTP);
                tdb_ftp_longFTP=tdb_ftp_longFTP(Indexa,:);

            end % end if condition there are new scans available

            if size(remcases_longFTP,1)>0

                fprintf(2,'Warning! %s cases were in the last Longitudinal FTP database and now are not found:\n\n', num2str(size(remcases,1)));

                choice = input(['\n\n Choose an action:' ,...
                    '\n     [1] Delete removed cases from present database',...
                    '\n     [2] Keep removed cases, I will check later'
                    '\n     --> ']);

                if choice==1

                    filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/Removed_Cases_FTP_CompWM_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
                    writetable(remcases_longFTP,filename,'WriteRowNames',true)
                    oldinfo_longFTP(Indexb,:)=[];

                elseif choice==2

                        fprintf(1,'Fine! I will still save a csv file to keep track of that.\n\n');
                        filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/RemovedbutKept_Cases_FTP_CompWM_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
                        writetable(remcases_longFTP,filename,'WriteRowNames',true)

                end % end if condition choice on removed cases

                clear choice filename

            end % end if condition existence of removed cases

            elseif checkinfo==1

                fprintf(2,'No new FTP scans needed Longitudinal extraction\n');
                tdb_ftp_longFTP={};

         end % end if condition there is a difference between old extraction and new extraction list

    end % end if condition there is a previous FTP Longitudinal extraction database available


    if size(tdb_ftp,1)>0

% 1. Cross-sectional Module

%%% Create the list of regions we are interested in, plus creating labels

fsids=[2;  4;	5;	7;	8;	10;	11;	12;	13;	14;	15;	16;	17;	18;	24;	26;	28;	30;	31;	41;	43;	44;	46;	47;	49;	50;	51;	52;	53;	54;	58;	60;	62;	63;	72;	77;	80;	85;	251;	252;	253;	254;	255;	1000;	1001;	1002;	1003;	1005;	1006;	1007;	1008;	1009;	1010;	1011;	1012;	1013;	1014;	1015;	1016;	1017;	1018;	1019;	1020;	1021;	1022;	1023;	1024;	1025;	1026;	1027;	1028;	1029;	1030;	1031;	1032;	1033;	1034;	1035;	2000;	2001;	2002;	2003;	2005;	2006;	2007;	2008;	2009;	2010;	2011;	2012;	2013;	2014;	2015;	2016;	2017;	2018;	2019;	2020;	2021;	2022;	2023;	2024;	2025;	2026;	2027;	2028;	2029;	2030;	2031;	2032;	2033;	2034;	2035];
fslabs={'Left_Cerebral_White_Matter'; 'Left_Lateral_Ventricle';	'Left_Inf_Lat_Vent';	'Left_Cerebellum_White_Matter';	'Left_Cerebellum_Cortex';	'Left_Thalamus_Proper';	'Left_Caudate';	'Left_Putamen';	'Left_Pallidum';	'Third_Ventricle';	'Fourth_Ventricle';	'Brain_Stem';	'Left_Hippocampus';	'Left_Amygdala';	'CSF';	'Left_Accumbens_area';	'Left_VentralDC';	'Left_vessel';	'Left_choroid_plexus';	'Right_Cerebral_White_Matter'; 'Right_Lateral_Ventricle';	'Right_Inf_Lat_Vent';	'Right_Cerebellum_White_Matter';	'Right_Cerebellum_Cortex';	'Right_Thalamus_Proper';	'Right_Caudate';	'Right_Putamen';	'Right_Pallidum';	'Right_Hippocampus';	'Right_Amygdala';	'Right_Accumbens_area';	'Right_VentralDC';	'Right_vessel';	'Right_choroid_plexus';	'Fifth_Ventricle';	'WM_hypointensities';	'non_WM_hypointensities';	'Optic_Chiasm';	'CC_Posterior';	'CC_Mid_Posterior';	'CC_Central';	'CC_Mid_Anterior';	'CC_Anterior';	'ctx_lh_unknown';	'ctx_lh_bankssts';	'ctx_lh_caudalanteriorcingulate';	'ctx_lh_caudalmiddlefrontal';	'ctx_lh_cuneus';	'ctx_lh_entorhinal';	'ctx_lh_fusiform';	'ctx_lh_inferiorparietal';	'ctx_lh_inferiortemporal';	'ctx_lh_isthmuscingulate';	'ctx_lh_lateraloccipital';	'ctx_lh_lateralorbitofrontal';	'ctx_lh_lingual';	'ctx_lh_medialorbitofrontal';	'ctx_lh_middletemporal';	'ctx_lh_parahippocampal';	'ctx_lh_paracentral';	'ctx_lh_parsopercularis';	'ctx_lh_parsorbitalis';	'ctx_lh_parstriangularis';	'ctx_lh_pericalcarine';	'ctx_lh_postcentral';	'ctx_lh_posteriorcingulate';	'ctx_lh_precentral';	'ctx_lh_precuneus';	'ctx_lh_rostralanteriorcingulate';	'ctx_lh_rostralmiddlefrontal';	'ctx_lh_superiorfrontal';	'ctx_lh_superiorparietal';	'ctx_lh_superiortemporal';	'ctx_lh_supramarginal';	'ctx_lh_frontalpole';	'ctx_lh_temporalpole';	'ctx_lh_transversetemporal';	'ctx_lh_insula';	'ctx_rh_unknown';	'ctx_rh_bankssts';	'ctx_rh_caudalanteriorcingulate';	'ctx_rh_caudalmiddlefrontal';	'ctx_rh_cuneus';	'ctx_rh_entorhinal';	'ctx_rh_fusiform';	'ctx_rh_inferiorparietal';	'ctx_rh_inferiortemporal';	'ctx_rh_isthmuscingulate';	'ctx_rh_lateraloccipital';	'ctx_rh_lateralorbitofrontal';	'ctx_rh_lingual';	'ctx_rh_medialorbitofrontal';	'ctx_rh_middletemporal';	'ctx_rh_parahippocampal';	'ctx_rh_paracentral';	'ctx_rh_parsopercularis';	'ctx_rh_parsorbitalis';	'ctx_rh_parstriangularis';	'ctx_rh_pericalcarine';	'ctx_rh_postcentral';	'ctx_rh_posteriorcingulate';	'ctx_rh_precentral';	'ctx_rh_precuneus';	'ctx_rh_rostralanteriorcingulate';	'ctx_rh_rostralmiddlefrontal';	'ctx_rh_superiorfrontal';	'ctx_rh_superiorparietal';	'ctx_rh_superiortemporal';	'ctx_rh_supramarginal';	'ctx_rh_frontalpole';	'ctx_rh_temporalpole';	'ctx_rh_transversetemporal';	'ctx_rh_insula'};

brklabsext={'Assigned_Braak_ADNIcutoffs','Braak_12_nohcp','Braak_12','Braak_34','Braak_56'};
brklabs={'Braak_12_nohcp','Braak_12','Braak_34','Braak_56'};

metaroilabext={'Assigned_MetaROI_ADNIcutoff_1p2','MetaROI'};
metaroilab={'MetaROI'};

% Generate list of images I will work with

vols=char(tdb_ftp(:,1));
vols_spm=tdb_ftp(:,1);
numv1= size(vols_spm,1);

aparcs = char(tdb_ftp(:,2));
aparc_spm = tdb_ftp(:,2);

%% Extraction section

regs = fsids; % We are extracting from all the 113 regions included in the object fsids
numf1=size(regs,1);

M = zeros(numv1,numf1); % create empty matrix in which to store values from all the ROIs
M_sz = zeros(numv1,numf1);

M_refreg=zeros(numv1,1); % creating an empty matrix storing ref region extraction
M_refreg_sz=zeros(numv1,1);

for v = 1:numv1 % Reading each image

r1n=spm_vol(aparcs(v,:)); % reading the aparc-aseg
r1=spm_read_vols(r1n);

img=spm_vol(vols(v,:)); % Reading img header
img1=spm_read_vols(img); % Reading img values

[~,f,e]=spm_fileparts(vols(v,:));
tempfname=char(strcat(f,e));
fprintf(1,'**Now Extracting from %s\n',tempfname);

%%%%%%%%%%%%%%% Braak Section %%%%%%%%%%%%%%%%

%%% Braak 12 no Hippocampus %%%

Mbrk12_nohcp = zeros(1,size(ind_brk12_nohcp,1)); % create empty matrix in which to store values from all the ROIs
Mbrk12_nohcp_sz = zeros(1,size(ind_brk12_nohcp,1));

for f = 1:size(ind_brk12_nohcp,1) % Looping through the ROIs
      reg = ind_brk12_nohcp(f);  % Selecting the ROI from the loop
      mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
      img_mask=img1.*mask; % Creating the masked image
      temp=nonzeros(img_mask);
      temp=temp(~isnan(temp));
      ext_val=mean(temp);
      vec_braak12_nohcp(f)=ext_val; % Store values for the individual ROIs across the images
      vec_braak12_nohcp_sz(f)=size(temp,1);
end % end for loop for Braak 1
      Mbrk12_nohcp(1,:)=vec_braak12_nohcp; % save the values ROI-wise (column-wise)
      Mbrk12_nohcp_sz(1,:)=vec_braak12_nohcp_sz;  % save the rois size

%%% Braak 12 with Hippocampus %%%

Mbrk12 = zeros(1,size(ind_brk12,1)); % create empty matrix in which to store values from all the ROIs
Mbrk12_sz = zeros(1,size(ind_brk12,1));

for f = 1:size(ind_brk12,1) % Looping through the ROIs
      reg = ind_brk12(f);  % Selecting the ROI from the loop
      mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
      img_mask=img1.*mask; % Creating the masked image
      temp=nonzeros(img_mask);
      temp=temp(~isnan(temp));
      ext_val=mean(temp);
      vec_braak12(f)=ext_val; % Store values for the individual ROIs across the images
      vec_braak12_sz(f)=size(temp,1);
end % end for loop Braak 12
      Mbrk12(1,:)=vec_braak12; % save the values ROI-wise (column-wise)
      Mbrk12_sz(1,:)=vec_braak12_sz;  % save the rois size

%%% Braak 34  %%%

Mbrk34 = zeros(1,size(ind_brk34,1)); % create empty matrix in which to store values from all the ROIs
Mbrk34_sz = zeros(1,size(ind_brk34,1));

for f = 1:size(ind_brk34,1) % Looping through the ROIs
      reg = ind_brk34(f);  % Selecting the ROI from the loop
      mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
      img_mask=img1.*mask; % Creating the masked image
      temp=nonzeros(img_mask);
      temp=temp(~isnan(temp));
      ext_val=mean(temp);
      vec_braak34(f)=ext_val; % Store values for the individual ROIs across the images
      vec_braak34_sz(f)=size(temp,1);
end % end for loop Braak 34
      Mbrk34(1,:)=vec_braak34; % save the values ROI-wise (column-wise)
      Mbrk34_sz(1,:)=vec_braak34_sz;  % save the rois size

%%% Braak 56  %%%

Mbrk56 = zeros(1,size(ind_brk56,1)); % create empty matrix in which to store values from all the ROIs
Mbrk56_sz = zeros(1,size(ind_brk56,1));

for f = 1:size(ind_brk56,1) % Looping through the ROIs
      reg = ind_brk56(f);  % Selecting the ROI from the loop
      mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
      img_mask=img1.*mask; % Creating the masked image
      temp=nonzeros(img_mask);
      temp=temp(~isnan(temp));
      ext_val=mean(temp);
      vec_braak56(f)=ext_val; % Store values for the individual ROIs across the images
      vec_braak56_sz(f)=size(temp,1);
end % end for loop Braak 56
      Mbrk56(1,:)=vec_braak56; % save the values ROI-wise (column-wise)
      Mbrk56_sz(1,:)=vec_braak56_sz;  % save the rois size

%%%%%%% Computing weighted means for the Braak
%%%%%%% stages at this point and saving sizes for
%%%%%%% each of them

wmean_brk12_nohcp = sum(Mbrk12_nohcp_sz.*Mbrk12_nohcp,2)./sum(Mbrk12_nohcp_sz,2);
wmean_brk12 = sum(Mbrk12_sz.*Mbrk12,2)./sum(Mbrk12_sz,2);
wmean_brk34 = sum(Mbrk34_sz.*Mbrk34,2)./sum(Mbrk34_sz,2);
wmean_brk56 = sum(Mbrk56_sz.*Mbrk56,2)./sum(Mbrk56_sz,2);
size_brk12_nohcp=sum(Mbrk12_nohcp_sz,2);
size_brk12=sum(Mbrk12_sz,2);
size_brk34=sum(Mbrk34_sz,2);
size_brk56=sum(Mbrk56_sz,2);

    if wmean_brk56(1)>1.15

    brkstg=56;

    elseif wmean_brk56(1)<1.15 && wmean_brk34(1)>1.25

    brkstg=34;

    elseif wmean_brk56(1)<1.15 && wmean_brk34(1)<1.25 && wmean_brk12(1)>1.31

    brkstg=12;

    elseif wmean_brk56(1)<1.15 && wmean_brk34(1)<1.25 && wmean_brk12(1)<1.31

    brkstg=0;

    else

    brkstg=999;

    end % end if condition to assign braak stage according to Maass et al. Neuroimage. Not provided anymore in our reports

    allbrks(v,:)=[brkstg wmean_brk12_nohcp wmean_brk12 wmean_brk34 wmean_brk56];
    allbrks_sz(v,:)=[size_brk12_nohcp size_brk12 size_brk34 size_brk56];

    %%%%%%%%%%%%%%% MetaROI Section %%%%%%%%%%%%%%
    %%% Jack et al., 2017 https://doi.org/10.1016/j.jalz.2016.08.005
    %%% Cut-off of 1.20 for ADNI data provided in Maass et al.,
    %%% 2017 http://dx.doi.org/10.1016/j.neuroimage.2017.05.058
    %%% Cut-off refers to Abeta+ MCI/AD dementia vs. Abeta-
    %%% old controls

    Mmetaroi = zeros(1,size(ind_metaroi,1)); % create empty matrix in which to store values from all the ROIs
    Mmetaroi_sz = zeros(1,size(ind_metaroi,1));

for f = 1:size(ind_metaroi,1) % Looping through the ROIs
      reg = ind_metaroi(f);  % Selecting the ROI from the loop
      mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
      img_mask=img1.*mask; % Creating the masked image
      temp=nonzeros(img_mask);
      temp=temp(~isnan(temp));
      ext_val=mean(temp);
      vec_metaroi(f)=ext_val; % Store values for the individual ROIs across the images
      vec_metaroi_sz(f)=size(temp,1);
end % end for loop metaROI
      Mmetaroi(1,:)=vec_metaroi; % save the values ROI-wise (column-wise)
      Mmetaroi_sz(1,:)=vec_metaroi_sz;  % save the rois size

    wmean_metaroi = sum(Mmetaroi_sz.*Mmetaroi,2)./sum(Mmetaroi_sz,2);
    size_metaroi=sum(Mmetaroi_sz,2);
    metaroi_pos=NaN;
%         if wmean_metaroi(1)>1.20
%          metaroi_pos=1;
%         else
%          metaroi_pos=0;
%         end % end if condition assign Status based on metaROI value of 1.2 based on Maass et al.

    defmetaroi(v,:)=[metaroi_pos wmean_metaroi];
    defmetaroi_sz(v,:)=size_metaroi;

for f = 1:numf1 % Looping through the ROIs
            reg = regs(f);  % Selecting the ROI from the loop
            mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
                  img_mask=img1.*mask; % Creating the masked image
                  temp=nonzeros(img_mask);
                  temp=temp(~isnan(temp));
                  ext_val=mean(temp);
                  vec(f)=ext_val; % Store values for the individual ROIs across the images
                  vec_sz(f)=size(temp,1);
end % end for loop for each Freesurfer APARC+ASEG region

  M(v,:)=vec; % save the values ROI-wise (column-wise)
  M_sz(v,:)=vec_sz;  % save the roi size

 %%%% Small module to get reference region
 %%%% values

 templdsid=regexp(vols_spm(v),'LDS\d{7}','match','once');
 tempftpdate=regexp(vols_spm(v),'\d{4}-\d{2}-\d{2}','match','once');
 tempftpfold=regexp(vols_spm(v),'\FTP_\d{4}-\d{2}-\d{2}','match','once');
 temptimepoint=regexp(vols_spm(v),'\Timepoint\d{1}','match','once');

 refreg=strcat('/mnt/coredata/Projects/LEADS/data_f7p1/processed/',templdsid,'/',temptimepoint,'/',tempftpfold,'/infcblg_ref_mask.nii');
 tempsuv=strcat('/mnt/coredata/Projects/LEADS/data_f7p1/processed/',templdsid,'/',temptimepoint,'/',tempftpfold,'/r',templdsid,'_FTP_',tempftpdate,'.nii');

 refregn=spm_vol(char(refreg));
 refreg1=spm_read_vols(refregn);
 refreg1=round(refreg1);

 suv=spm_vol(char(tempsuv));
 suv1=spm_read_vols(suv);

 mask = (refreg1 ==1); % creating the mask
 suv_mask=suv1.*mask; % Creating the masked image
 temp=nonzeros(suv_mask);
 temp=temp(~isnan(temp));
 refregval=mean(temp);
 refreg_sz=size(temp,1);

 M_refreg(v,1)=refregval;
 M_refreg_sz(v,1)=refreg_sz;

 clear r1n r1 img img1 p f e Mbrk12_nohcp Mbrk12_nohcp_sz Mbrk12 Mbrk12_sz Mbrk34 Mbrk34_sz Mbrk56 Mbrk56_sz Mmetaroi Mmetaroi_sz templdsid tempftpdate tempftpfold temptimepoint refreg tempsuv refregn regreg1 suv suv1 mask suv_mask temp refregval refreg_sz metaroi_pos wmean_metaroi size_metaroi brkstg wmean_brk12_nohcp wmean_brk12 wmean_brk34 wmean_brk56 size_brk12_nohcp size_brk12 size_brk34 size_brk56

end % end for loop for each image

T_refreg=array2table(M_refreg);
T_refreg_sz=array2table(M_refreg_sz);
T_refreg.Properties.VariableNames={'ScalingFactor_InfCerebGray'};
T_refreg_sz.Properties.VariableNames={'ScalingFactor_InfCerebGray_ClustSize'};

T_braaks=array2table(allbrks);
T_braaks_sz=array2table(allbrks_sz);
T_braaks.Properties.VariableNames = brklabsext;
T_braaks_sz.Properties.VariableNames = cellstr(strcat(brklabs,'_ClustSize'));

T_metarois=array2table(defmetaroi);
T_metarois_sz=array2table(defmetaroi_sz);
T_metarois.Properties.VariableNames=metaroilabext;
T_metarois_sz.Properties.VariableNames = cellstr(strcat(metaroilab,'_ClustSize'));

T = array2table(M); % create Tables for value reports
T_sz=array2table(M_sz);
T.Properties.VariableNames = fslabs;
T_sz.Properties.VariableNames = cellstr(strcat(fslabs,'_ClustSize'));


templdsids=regexp(vols_spm,'LDS\d{7}','match','once');
tempdates=regexp(vols_spm,'\d{4}-\d{2}-\d{2}','match','once');
tempmridates=regexp(aparc_spm,'(?<=MRI-T1_)\d{4}-\d{2}-\d{2}','match','once');


meta=horzcat(templdsids, tempdates, tempmridates);
T_meta=array2table(meta);
T_meta.Properties.VariableNames={'ID','FTPPET_Date','MRI_Date'};

qc=horzcat(vols_spm, aparc_spm);
T_qc=array2table(qc);
T_qc.Properties.VariableNames={'FTPPET_path','APARC_path'};


T=[T_meta T_qc T_refreg T_metarois T_braaks T T_refreg_sz T_metarois_sz T_braaks_sz T_sz];
if size(newcases,1)>0
T=vertcat(oldinfo,T);
end

filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/FTP_ROI_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
writetable(T,filename,'WriteRowNames',true)
copyfile(filename, '/shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/');

    end % end if condition new FTP scans needed extractions

    if size(tdb_ftp_longFTP,1)>0

% 2. Longitudinal Module

fprintf(1,'Now Starting longitudinal module for FTP...\n');

ftpscans_longFTP=char(tdb_ftp_longFTP(:,1)); ftpscans_longFTP_spm=cellstr(ftpscans_longFTP);
aparcscans_longFTP=char(tdb_ftp_longFTP(:,2)); aparcscans_longFTP_spm=cellstr(aparcscans_longFTP);

M_longFTP = zeros(size(ftpscans_longFTP,1),1); % create empty matrix in which to store values from all the ROIs
M_suvr_longFTP=zeros(size(ftpscans_longFTP,1),1);
M_sz_longFTP = zeros(size(ftpscans_longFTP,1),1);

for i=1:size(ftpscans_longFTP,1)

    fprintf(1,'Now starting Longitudinal processing for %s\n',ftpscans_longFTP_spm{i});

    tempimg=ftpscans_longFTP_spm{i};
    tempaparc=aparcscans_longFTP_spm{i};

    % working on the pieces needed
    % 1. WM parcellation from the aparc, smooth and threshold

   [appath,apfname,~]=spm_fileparts(char(tempaparc));
    wmaparcfname=strcat(appath,'/',apfname,'_wm.nii');
    swmaparcfname=strcat(appath,'/s',apfname,'_wm.nii');
    swmaparcfnamethres=strcat(appath,'/s',apfname,'_wm_thr0p7.nii');
    tempnu=strcat(appath,'/',apfname(1:29),'nu.nii');

    spm('defaults','PET');
    matlabbatch{1}.spm.util.imcalc.input = cellstr(tempaparc);
    matlabbatch{1}.spm.util.imcalc.output = char(wmaparcfname);
    matlabbatch{1}.spm.util.imcalc.outdir = {''};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1==2 | i1==41';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 0;
    matlabbatch{1}.spm.util.imcalc.options.dtype = spm_type('float32');
    matlabbatch{2}.spm.spatial.smooth.data = cellstr(wmaparcfname);
    matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{2}.spm.spatial.smooth.dtype = 0;
    matlabbatch{2}.spm.spatial.smooth.im = 0;
    matlabbatch{2}.spm.spatial.smooth.prefix = 's';
    matlabbatch{3}.spm.util.imcalc.input = cellstr(swmaparcfname);
    matlabbatch{3}.spm.util.imcalc.output = char(swmaparcfnamethres);
    matlabbatch{3}.spm.util.imcalc.outdir = {''};
    matlabbatch{3}.spm.util.imcalc.expression = 'i1>0.7';
    matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{3}.spm.util.imcalc.options.mask = 0;
    matlabbatch{3}.spm.util.imcalc.options.interp = 0;
    matlabbatch{3}.spm.util.imcalc.options.dtype = spm_type('uint8');
    spm_jobman('run',matlabbatch); clear matlabbatch;

    copyfile(char(swmaparcfnamethres),char(strcat(appath,'/erodedwm_ref_mask.nii')));

    % Send QC multislice to the shared petcore if it does not exist already

    srcmulti=dir(strcat(path_qccompwm,'/*',apfname,'*'));

    if size(srcmulti,1)==0

    tempbs_wm_slovname=quickmultislice3(tempnu, char(swmaparcfnamethres),'axial','nih.lut','0.5 1.5','-30 6 58');
    copyfile(tempbs_wm_slovname,path_qccompwm);

    end


    % Grab the rSCAN to get the scaling factor to be saved in the database

    [imgpath,imgfname,~]=spm_fileparts(char(tempimg));

    rimgfname=strcat(imgpath,'/r',imgfname(1:end-13),'.nii');

    % Extract value from r image

    roi=spm_vol(char(strcat(appath,'/erodedwm_ref_mask.nii')));
    roi1=spm_read_vols(roi);

      img=spm_vol(rimgfname); % Reading img header from the loop
      img1=spm_read_vols(img); % Reading img values
      img_mask=img1.*roi1; % Creating the masked image
      img_mask=nonzeros(img_mask);
      img_mask=img_mask(~isnan(img_mask));
      ext_val=mean(img_mask);
      M_longFTP(i,1)=ext_val;

  % Extract value from suvr image

  imgsuvr=spm_vol(char(tempimg)); % Reading imgsuvr header from the loop
  imgsuvr1=spm_read_vols(imgsuvr); % Reading imgsuvr values
  imgsuvr_mask=imgsuvr1.*roi1; % Creating the masked image
  imgsuvr_mask=nonzeros(imgsuvr_mask);
  imgsuvr_mask=imgsuvr_mask(~isnan(imgsuvr_mask));
  ext_val_suvr=mean(imgsuvr_mask);
  M_suvr_longFTP(i,1)=ext_val_suvr;

  M_sz_longFTP(i,1)=size(img_mask,1);

  % Create the new suvr image, last thing I want

    exp=char(strcat('i1/',num2str(ext_val)));
    newfname=char(strcat(imgpath,'/',imgfname(1:end-13),'_suvr_erodedWM.nii'));

    spm('defaults','PET');
    matlabbatch{1}.spm.util.imcalc.input = cellstr(rimgfname);
    matlabbatch{1}.spm.util.imcalc.output = newfname;
    matlabbatch{1}.spm.util.imcalc.outdir = {''};
    matlabbatch{1}.spm.util.imcalc.expression = exp;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = spm_type('float32');
    spm_jobman('run',matlabbatch); clear matlabbatch;

    clear tempnu srcmulti apfname appath bsaparcfname compwmfname compwmfs exp ext_val ext_val_suvr img img1 imgsuvr imgsuvr1 imgsuvr_mask img_mask imgfname imgpath newfname rimgfname roi roi1 swmaparcfname swmaparcfnamethres tempaparc tempimg wmaparcfname

end

T_longFTP=array2table(M_longFTP);
Tsuvr_longFTP=array2table(M_suvr_longFTP);
T_sz_longFTP=array2table(M_sz_longFTP);
T_longFTP.Properties.VariableNames={'ScalingFactor_ErodedWM'};
Tsuvr_longFTP.Properties.VariableNames={'SUVR_ErodedWM'};
T_sz_longFTP.Properties.VariableNames={'ScalingFactor_ErodedWM_ClustSize'};

templdsids_longFTP=regexp(ftpscans_longFTP_spm,'LDS\d{7}','match','once');
tempdates_longFTP=regexp(ftpscans_longFTP_spm,'\d{4}-\d{2}-\d{2}','match','once');
tempmridates_longFTP=regexp(aparcscans_longFTP_spm,'(?<=MRI-T1_)\d{4}-\d{2}-\d{2}','match','once');

meta_longFTP=horzcat(templdsids_longFTP, tempdates_longFTP, tempmridates_longFTP);
T_meta_longFTP=array2table(meta_longFTP);
T_meta_longFTP.Properties.VariableNames={'ID','FTPPET_Date','MRI_Date'};

qc_longFTP=horzcat(ftpscans_longFTP_spm, aparcscans_longFTP_spm);
T_qc_longFTP=array2table(qc_longFTP);
T_qc_longFTP.Properties.VariableNames={'FTPPET_path','APARC_path'};

T_longFTP=[T_meta_longFTP T_qc_longFTP T_longFTP Tsuvr_longFTP T_sz_longFTP];

if size(newcases,1)>0
T_longFTP=vertcat(oldinfo_longFTP,T_longFTP);
end

filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/FTP_ErodedWM_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
writetable(T_longFTP,filename,'WriteRowNames',true)
copyfile(filename,'/shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/');

    end % end if condition new FTP scans needed Longitudinal extractions

% Done, wrapping up - clear everything before proceeding to the FDG ROI
% Extraction

clearvars -except tdb_fdg path_processed path_extraction path_qccompwm
fprintf(1,'**Completed ROI and Longitudinal extraction for all the FTP SUVR images available!\n');
