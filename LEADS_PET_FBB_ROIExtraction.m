% Create a temp table with the info from the "updated" list of scans
% available in the processed folder
% Sep 2021: Code cleanup

% Feedback for the user

fprintf(1,'Checking what FBB scans are new...\n');

% Store IDs for regions that will be used later

ind_accpcc=[1002 1010 1023 1026 2002 2010 2023 2026]';
ind_front=[1003 1032 1012 1014 1018 1019 1020 1027 1028 2003 2032 2012 2014 2018 2019 2020 2027 2028]';
ind_temp=[1015 1030 2015 2030]';
ind_pariet=[1008 1025 1029 1031 2008 2025 2029 2031]';
ind_cbl=[7 8 46 47]';

% Generate the entire list of images now available at all timepoints in the
% processed folder

newldsids=regexp(tdb_fbb(:,1),'LDS\d{7}','match','once');
newldsdates=regexp(tdb_fbb(:,1),'\d{4}-\d{2}-\d{2}','match','once');
newinfo=horzcat(newldsids,newldsdates);
newinfo=array2table(newinfo);
newinfo.Properties.VariableNames={'ID','FBBPET_Date'};

%% Let's save the tdb_fbb list also for the longitudinal ref processing
% we will use this to replicate the approach of reading last database and
% add on top also for the long database

tdb_fbb_longFBB=tdb_fbb;

%% Read the last extraction database, this will be used to know what's new %%
%% This command reads all the files in the extraction folder and then reads the latest one %%
%% The final aim is to cleanup our list to only extract from new scans. If no extraction files are found then it will extract from all the scans

olddb=dir('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/FBB_ROI_Extraction*');

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

    checkinfo=isequal(newinfo(:,{'ID','FBBPET_Date'}),oldinfo(:,{'ID','FBBPET_Date'}));% Check if the b-ids columns in the new and old tables is unchanged

        if checkinfo==0

        Indexa=not(ismember(newinfo(:,{'ID','FBBPET_Date'}),oldinfo(:,{'ID','FBBPET_Date'}))); %% logical indexing for the new cases. Both are based on BID+date to include longitudinal scans
        Indexb=not(ismember(oldinfo(:,{'ID','FBBPET_Date'}),newinfo(:,{'ID','FBBPET_Date'}))); %% logical indexing for the removed cases

        newcases=newinfo(Indexa,:); %% Creates a table with the New FBB scans
        remcases=oldinfo(Indexb,:); %% Creates a table with eventually removed FBB scans

            if size(newcases,1)>0

                fprintf(2,'**New FBB-PET scans were found for Cross-Sectional Processing:\n\n');
                disp(newcases);
                tdb_fbb=tdb_fbb(Indexa,:);

            end % end if condition there are new scans available

            if size(remcases,1)>0

                fprintf(2,'Warning! %s cases were in the last FBB database and now are not found:\n\n', num2str(size(remcases,1)));

                choice = input(['\n\n Choose an action:' ,...
                    '\n     [1] Delete removed cases from present database',...
                    '\n     [2] Keep removed cases, I will check later'
                    '\n     --> ']);

                if choice==1

                    filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/Removed_Cases_FBB_ROI_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
                    writetable(remcases,filename,'WriteRowNames',true)
                    oldinfo(Indexb,:)=[];

                elseif choice==2

                        fprintf(1,'Fine! I will still save a csv file to keep track of that.\n\n');
                        filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/RemovedbutKept_Cases_FBB_ROI_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
                        writetable(remcases,filename,'WriteRowNames',true)
                end % end if condition choice on removed cases
                
                clear choice filename

            end % end if condition existence of removed cases
            
            elseif checkinfo==1
                
                fprintf(2,'No new FBB scans needed extraction!\n');
                tdb_fbb={};
                
         end % end if condition there is a difference between old extraction and new extraction list
        
    end % end if condition there is a previous FBB extraction database available
    
% To be safe, let's separately scan also the longitudinal database

olddb_longFBB=dir('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/FBB_CompWM_Extraction*');

    if size(olddb_longFBB,1)>0
    
    olddb_longFBB=struct2cell(olddb_longFBB)';
    olddb_longFBB=olddb_longFBB(:,[1 3]);
    olddb_longFBB=array2table(olddb_longFBB);
    olddb_longFBB.Properties.VariableNames(1) = cellstr(strcat('Filename'));
    olddb_longFBB.Properties.VariableNames(2) = cellstr(strcat('DateCreated'));
    olddb_longFBB.DateCreated=datetime(olddb_longFBB.DateCreated);
    filt=max(olddb_longFBB.DateCreated); % Store most recent date
    olddb_longFBB=olddb_longFBB(olddb_longFBB.DateCreated==filt,:); % Subset to select the most recent file
    oldinfo_longFBB=readtable(strcat('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/',char(table2cell(olddb_longFBB(1,1))))); %% ready to check differences now

    newinfo=sortrows(newinfo);
    oldinfo_longFBB=sortrows(oldinfo_longFBB);

    checkinfo=isequal(newinfo(:,{'ID','FBBPET_Date'}),oldinfo_longFBB(:,{'ID','FBBPET_Date'}));% Check if the b-ids columns in the new and old tables is unchanged

        if checkinfo==0

        Indexa=not(ismember(newinfo(:,{'ID','FBBPET_Date'}),oldinfo_longFBB(:,{'ID','FBBPET_Date'}))); %% logical indexing for the new cases. Both are based on BID+date to include longitudinal scans
        Indexb=not(ismember(oldinfo_longFBB(:,{'ID','FBBPET_Date'}),newinfo(:,{'ID','FBBPET_Date'}))); %% logical indexing for the removed cases

        newcases_longFBB=newinfo(Indexa,:); %% Creates a table with the New FBB scans
        remcases_longFBB=oldinfo_longFBB(Indexb,:); %% Creates a table with eventually removed FBB scans

            if size(newcases_longFBB,1)>0

                fprintf(2,'**New FBB-PET scans were found for Longitudinal Processing:\n\n');
                disp(newcases_longFBB);
                tdb_fbb_longFBB=tdb_fbb_longFBB(Indexa,:);

            end % end if condition there are new scans available

            if size(remcases_longFBB,1)>0

                fprintf(2,'Warning! %s cases were in the last Longitudinal FBB database and now are not found:\n\n', num2str(size(remcases,1)));

                choice = input(['\n\n Choose an action:' ,...
                    '\n     [1] Delete removed cases from present database',...
                    '\n     [2] Keep removed cases, I will check later'
                    '\n     --> ']);

                if choice==1

                    filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/Removed_Cases_FBB_CompWM_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
                    writetable(remcases_longFBB,filename,'WriteRowNames',true)
                    oldinfo_longFBB(Indexb,:)=[];

                elseif choice==2

                        fprintf(1,'Fine! I will still save a csv file to keep track of that.\n\n');
                        filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/RemovedbutKept_Cases_FBB_CompWM_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
                        writetable(remcases_longFBB,filename,'WriteRowNames',true)
                        
                end % end if condition choice on removed cases
                
                clear choice filename

            end % end if condition existence of removed cases
                        
            elseif checkinfo==1
                
                fprintf(2,'No new FBB scans needed longitudinal extraction!\n');
                tdb_fbb_longFBB={};
                
         end % end if condition there is a difference between old extraction and new extraction list
        
    end % end if condition there is a previous FBB Longitudinal extraction database available
    


    if size(tdb_fbb,1)>0
    
%%% Time to extract with lists we have

% 1. Cross-sectional Module
    
%%% Create the list of regions we are interested in, plus creating labels

fsids=[2; 4;	5;	7;	8;	10;	11;	12;	13;	14;	15;	16;	17;	18;	24;	26;	28;	30;	31;	41;	43;	44;	46;	47;	49;	50;	51;	52;	53;	54;	58;	60;	62;	63;	72;	77;	80;	85;	251;	252;	253;	254;	255;	1000;	1001;	1002;	1003;	1005;	1006;	1007;	1008;	1009;	1010;	1011;	1012;	1013;	1014;	1015;	1016;	1017;	1018;	1019;	1020;	1021;	1022;	1023;	1024;	1025;	1026;	1027;	1028;	1029;	1030;	1031;	1032;	1033;	1034;	1035;	2000;	2001;	2002;	2003;	2005;	2006;	2007;	2008;	2009;	2010;	2011;	2012;	2013;	2014;	2015;	2016;	2017;	2018;	2019;	2020;	2021;	2022;	2023;	2024;	2025;	2026;	2027;	2028;	2029;	2030;	2031;	2032;	2033;	2034;	2035];
fslabs={'Left_Cerebral_White_Matter'; 'Left_Lateral_Ventricle';	'Left_Inf_Lat_Vent';	'Left_Cerebellum_White_Matter';	'Left_Cerebellum_Cortex';	'Left_Thalamus_Proper';	'Left_Caudate';	'Left_Putamen';	'Left_Pallidum';	'Third_Ventricle';	'Fourth_Ventricle';	'Brain_Stem';	'Left_Hippocampus';	'Left_Amygdala';	'CSF';	'Left_Accumbens_area';	'Left_VentralDC';	'Left_vessel';	'Left_choroid_plexus';	'Right_Cerebral_White_Matter'; 'Right_Lateral_Ventricle';	'Right_Inf_Lat_Vent';	'Right_Cerebellum_White_Matter';	'Right_Cerebellum_Cortex';	'Right_Thalamus_Proper';	'Right_Caudate';	'Right_Putamen';	'Right_Pallidum';	'Right_Hippocampus';	'Right_Amygdala';	'Right_Accumbens_area';	'Right_VentralDC';	'Right_vessel';	'Right_choroid_plexus';	'Fifth_Ventricle';	'WM_hypointensities';	'non_WM_hypointensities';	'Optic_Chiasm';	'CC_Posterior';	'CC_Mid_Posterior';	'CC_Central';	'CC_Mid_Anterior';	'CC_Anterior';	'ctx_lh_unknown';	'ctx_lh_bankssts';	'ctx_lh_caudalanteriorcingulate';	'ctx_lh_caudalmiddlefrontal';	'ctx_lh_cuneus';	'ctx_lh_entorhinal';	'ctx_lh_fusiform';	'ctx_lh_inferiorparietal';	'ctx_lh_inferiortemporal';	'ctx_lh_isthmuscingulate';	'ctx_lh_lateraloccipital';	'ctx_lh_lateralorbitofrontal';	'ctx_lh_lingual';	'ctx_lh_medialorbitofrontal';	'ctx_lh_middletemporal';	'ctx_lh_parahippocampal';	'ctx_lh_paracentral';	'ctx_lh_parsopercularis';	'ctx_lh_parsorbitalis';	'ctx_lh_parstriangularis';	'ctx_lh_pericalcarine';	'ctx_lh_postcentral';	'ctx_lh_posteriorcingulate';	'ctx_lh_precentral';	'ctx_lh_precuneus';	'ctx_lh_rostralanteriorcingulate';	'ctx_lh_rostralmiddlefrontal';	'ctx_lh_superiorfrontal';	'ctx_lh_superiorparietal';	'ctx_lh_superiortemporal';	'ctx_lh_supramarginal';	'ctx_lh_frontalpole';	'ctx_lh_temporalpole';	'ctx_lh_transversetemporal';	'ctx_lh_insula';	'ctx_rh_unknown';	'ctx_rh_bankssts';	'ctx_rh_caudalanteriorcingulate';	'ctx_rh_caudalmiddlefrontal';	'ctx_rh_cuneus';	'ctx_rh_entorhinal';	'ctx_rh_fusiform';	'ctx_rh_inferiorparietal';	'ctx_rh_inferiortemporal';	'ctx_rh_isthmuscingulate';	'ctx_rh_lateraloccipital';	'ctx_rh_lateralorbitofrontal';	'ctx_rh_lingual';	'ctx_rh_medialorbitofrontal';	'ctx_rh_middletemporal';	'ctx_rh_parahippocampal';	'ctx_rh_paracentral';	'ctx_rh_parsopercularis';	'ctx_rh_parsorbitalis';	'ctx_rh_parstriangularis';	'ctx_rh_pericalcarine';	'ctx_rh_postcentral';	'ctx_rh_posteriorcingulate';	'ctx_rh_precentral';	'ctx_rh_precuneus';	'ctx_rh_rostralanteriorcingulate';	'ctx_rh_rostralmiddlefrontal';	'ctx_rh_superiorfrontal';	'ctx_rh_superiorparietal';	'ctx_rh_superiortemporal';	'ctx_rh_supramarginal';	'ctx_rh_frontalpole';	'ctx_rh_temporalpole';	'ctx_rh_transversetemporal';	'ctx_rh_insula'};

macroslabsext={'GlobalSUVR_Type2', 'GlobalSUVR','ACC_PCC','Frontal','Temporal','Parietal','WholeCerebellum'};
macroslabs={'ACC_PCC','Frontal','Temporal','Parietal','WholeCerebellum'};

% Generate list of images I will work with

vols=char(tdb_fbb(:,1));
vols_spm=tdb_fbb(:,1);
numv1= size(vols_spm,1);

aparcs = char(tdb_fbb(:,2));
aparc_spm = tdb_fbb(:,2);

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

Maccpcc = zeros(1,size(ind_accpcc,1)); % create empty matrix in which to store values from all the ROIs
Maccpcc_sz = zeros(1,size(ind_accpcc,1));

for f = 1:size(ind_accpcc,1) % Looping through the ROIs
      reg = ind_accpcc(f);  % Selecting the ROI from the loop
      mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
      img_mask=img1.*mask; % Creating the masked image
      temp=nonzeros(img_mask);
      temp=temp(~isnan(temp));
      ext_val=mean(temp);
      vec_amy(f)=ext_val; % Store values for the individual ROIs across the images
      vec_amy_sz(f)=size(temp,1);             
end
      Maccpcc(1,:)=vec_amy; % save the values ROI-wise (column-wise)
      Maccpcc_sz(1,:)=vec_amy_sz;  % save the roi size
clear vec_amy;
clear vec_amy_sz;

Mfront = zeros(1,size(ind_front,1)); % create empty matrix in which to store values from all the ROIs
Mfront_sz = zeros(1,size(ind_front,1));

for f = 1:size(ind_front,1) % Looping through the ROIs
      reg = ind_front(f);  % Selecting the ROI from the loop
      mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
      img_mask=img1.*mask; % Creating the masked image
      temp=nonzeros(img_mask);
      temp=temp(~isnan(temp));
      ext_val=mean(temp);
      vec_amy(f)=ext_val; % Store values for the individual ROIs across the images
      vec_amy_sz(f)=size(temp,1);             
end
      Mfront(1,:)=vec_amy; % save the values ROI-wise (column-wise)
      Mfront_sz(1,:)=vec_amy_sz;  % save the roi size
clear vec_amy;
clear vec_amy_sz;

Mtemp = zeros(1,size(ind_temp,1)); % create empty matrix in which to store values from all the ROIs
Mtemp_sz = zeros(1,size(ind_temp,1));

for f = 1:size(ind_temp,1) % Looping through the ROIs
      reg = ind_temp(f);  % Selecting the ROI from the loop
      mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
      img_mask=img1.*mask; % Creating the masked image
      temp=nonzeros(img_mask);
      temp=temp(~isnan(temp));
      ext_val=mean(temp);
      vec_amy(f)=ext_val; % Store values for the individual ROIs across the images
      vec_amy_sz(f)=size(temp,1);             
end
      Mtemp(1,:)=vec_amy; % save the values ROI-wise (column-wise)
      Mtemp_sz(1,:)=vec_amy_sz;  % save the roi size
clear vec_amy;
clear vec_amy_sz;

Mpariet = zeros(1,size(ind_pariet,1)); % create empty matrix in which to store values from all the ROIs
Mpariet_sz = zeros(1,size(ind_pariet,1));

for f = 1:size(ind_pariet,1) % Looping through the ROIs
      reg = ind_pariet(f);  % Selecting the ROI from the loop
      mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
      img_mask=img1.*mask; % Creating the masked image
      temp=nonzeros(img_mask);
      temp=temp(~isnan(temp));
      ext_val=mean(temp);
      vec_amy(f)=ext_val; % Store values for the individual ROIs across the images
      vec_amy_sz(f)=size(temp,1);             
end
      Mpariet(1,:)=vec_amy; % save the values ROI-wise (column-wise)
      Mpariet_sz(1,:)=vec_amy_sz;  % save the roi size

clear vec_amy;
clear vec_amy_sz;

Mcbl = zeros(1,size(ind_cbl,1)); % create empty matrix in which to store values from all the ROIs
Mcbl_sz = zeros(1,size(ind_cbl,1));

for f = 1:size(ind_cbl,1) % Looping through the ROIs
      reg = ind_cbl(f);  % Selecting the ROI from the loop
      mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
      img_mask=img1.*mask; % Creating the masked image
      temp=nonzeros(img_mask);
      temp=temp(~isnan(temp));
      ext_val=mean(temp);
      vec_amy(f)=ext_val; % Store values for the individual ROIs across the images
      vec_amy_sz(f)=size(temp,1);             
end
      Mcbl(1,:)=vec_amy; % save the values ROI-wise (column-wise)
      Mcbl_sz(1,:)=vec_amy_sz;  % save the roi size

    clear vec_amy;
    clear vec_amy_sz;

    %%% Get weighted mean values %%%
    wmean_accpcc = sum(Maccpcc_sz.*Maccpcc,2)./sum(Maccpcc_sz,2);
    wmean_front = sum(Mfront_sz.*Mfront,2)./sum(Mfront_sz,2);
    wmean_temp = sum(Mtemp_sz.*Mtemp,2)./sum(Mtemp_sz,2);
    wmean_pariet = sum(Mpariet_sz.*Mpariet,2)./sum(Mpariet_sz,2);
    wmean_cbl = sum(Mcbl_sz.*Mcbl,2)./sum(Mcbl_sz,2);
    size_accpcc=sum(Maccpcc_sz,2);
    size_front=sum(Mfront_sz,2);
    size_temp=sum(Mtemp_sz,2);
    size_pariet=sum(Mpariet_sz,2);
    size_cbl=sum(Mcbl_sz,2);

allmacros_comp=[wmean_accpcc wmean_front wmean_temp wmean_pariet];
compsuvr=(mean(allmacros_comp))/wmean_cbl;

%%% Module for new global score 

[sz1,sz2,sz3]=size(img1);
rsuvr=reshape(img1,sz1*sz2*sz3,1);
rparc=reshape(r1,sz1*sz2*sz3,1);

ind2=find((rparc==1003 | rparc==1012 | rparc==1014 ...
    | rparc==1018 | rparc==1019 | rparc==1020 ...
    | rparc==1027 | rparc==1028 | rparc==1032 ...
    | rparc==1009 | rparc==1015 | rparc==1030 ...
    | rparc==1008 | rparc==1025 | rparc==1029 ...
    | rparc==1031 | rparc==1002 | rparc==1023 | rparc==1010 ...
    | rparc==1026 | rparc==2003 | rparc==2012 | rparc==2014 ...
    | rparc==2018 | rparc==2019 | rparc==2020 ...
    | rparc==2027 | rparc==2028 | rparc==2032 ...
    | rparc==2009 | rparc==2015 | rparc==2030 ...
    | rparc==2008 | rparc==2025 | rparc==2029 ...
    | rparc==2031 | rparc==2002 | rparc==2023 | rparc==2010 ...
    | rparc==2026) & rsuvr>0);

compsuvr_t2=mean(rsuvr(ind2));

%%% 

allmacros(v,:)=[compsuvr_t2 compsuvr wmean_accpcc wmean_front wmean_temp wmean_pariet wmean_cbl];
allmacros_sz(v,:)=[size_accpcc size_front size_temp size_pariet size_cbl];

% Now extract from all aparc parcellation

for f = 1:numf1 % Looping through the ROIs
  reg = regs(f);  % Selecting the ROI from the loop
  mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
  img_mask=img1.*mask; % Creating the masked image
  temp=nonzeros(img_mask);
  temp=temp(~isnan(temp));
  ext_val=mean(temp);
  vec(f)=ext_val; % Store values for the individual ROIs across the images
  vec_sz(f)=size(temp,1);
end % end for loop for each ROI in aparc+aseg

% store subrois info

M(v,:)=vec; % save the values ROI-wise (column-wise)
M_sz(v,:)=vec_sz;  % save the roi size

%%%% Small module to get reference region
%%%% values 

templdsid=regexp(vols_spm(v),'LDS\d{7}','match','once');
tempfbbdate=regexp(vols_spm(v),'\d{4}-\d{2}-\d{2}','match','once');
tempfbbfold=regexp(vols_spm(v),'\FBB_\d{4}-\d{2}-\d{2}','match','once');
temptimepoint=regexp(vols_spm(v),'\Timepoint\d{1}','match','once');

refreg=strcat('/mnt/coredata/Projects/LEADS/data_f7p1/processed/',templdsid,'/',temptimepoint,'/',tempfbbfold,'/wholecbl_ref_mask.nii');
tempsuv=strcat('/mnt/coredata/Projects/LEADS/data_f7p1/processed/',templdsid,'/',temptimepoint,'/',tempfbbfold,'/r',templdsid,'_FBB_',tempfbbdate,'.nii');

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

clear r1n r1 img img1 p f e Maccpcc Maccpcc_sz Mfront Mfront_sz Mtemp Mtemp_sz Mpariet Mpariet_sz Mcbl Mcbl_sz wmean_accpcc wmean_front wmean_temp wmean_pariet wmean_cbl size_accpcc size_front size_temp size_pariet size_cbl allmacros_comp compsuvr rsuvr raparc ind2 compsuvr_t2 templdsid tempfbbdate tempfbbfold temptimepoint refreg tempsuv refregn regreg1 suv suv1 mask suv_mask temp refregval refreg_sz 
             
end

% creating the tables that will be merged together and concatenated to the
% final one

T_refreg=array2table(M_refreg);
T_refreg_sz=array2table(M_refreg_sz);
T_refreg.Properties.VariableNames={'ScalingFactor_WholeCereb'};
T_refreg_sz.Properties.VariableNames={'ScalingFactor_WholeCereb_ClustSize'};

T_macros=array2table(allmacros);
T_macros_sz=array2table(allmacros_sz);
T_macros.Properties.VariableNames = macroslabsext;
T_macros_sz.Properties.VariableNames = cellstr(strcat(macroslabs,'_ClustSize'));

T = array2table(M); % create Tables for value reports
T_sz=array2table(M_sz);
T.Properties.VariableNames = fslabs;
T_sz.Properties.VariableNames = cellstr(strcat(fslabs,'_ClustSize'));

templdsids=regexp(vols_spm,'LDS\d{7}','match','once');
tempdates=regexp(vols_spm,'\d{4}-\d{2}-\d{2}','match','once');
tempmridates=regexp(aparc_spm,'(?<=MRI_T1_)\d{4}-\d{2}-\d{2}','match','once');

meta=horzcat(templdsids, tempdates, tempmridates);
T_meta=array2table(meta);
T_meta.Properties.VariableNames={'ID','FBBPET_Date','MRI_Date'};

qc=horzcat(vols_spm, aparc_spm);
T_qc=array2table(qc);
T_qc.Properties.VariableNames={'FBBPET_path','APARC_path'};

T=[T_meta T_qc T_refreg T_macros T T_refreg_sz T_macros_sz T_sz]; 

if size(newcases,1)>0
T=vertcat(oldinfo,T);
end % end if condition this is a batch of new images that have to be appended to the existing one

% export database
filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/FBB_ROI_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
writetable(T,filename,'WriteRowNames',true)
copyfile(filename,'/shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/'); % copy in the shared drive for report generation in R            
 
    end % end if condition new images need extraction

% 2. Longitudinal Module

if size(tdb_fbb_longFBB,1)>0

fprintf(1,'Now Starting longitudinal module for FBB...\n');

fbbscans_longFBB=char(tdb_fbb_longFBB(:,1)); fbbscans_longFBB_spm=cellstr(fbbscans_longFBB);
aparcscans_longFBB=char(tdb_fbb_longFBB(:,2)); aparcscans_longFBB_spm=cellstr(aparcscans_longFBB);

M_longFBB = zeros(size(fbbscans_longFBB,1),1); % create empty matrix in which to store values from all the ROIs
M_suvr_longFBB=zeros(size(fbbscans_longFBB,1),1);
M_sz_longFBB = zeros(size(fbbscans_longFBB,1),1);

for i=1:size(fbbscans_longFBB,1)
   
    fprintf(1,'Now starting Longitudinal processing for %s\n',fbbscans_longFBB_spm{i});
    
    tempimg=fbbscans_longFBB_spm{i};
    tempaparc=aparcscans_longFBB_spm{i};
    
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
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
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
    matlabbatch{3}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch); clear matlabbatch;
    
    % Send QC multislice to the shared petcore
    
    tempbs_wm_slovname=quickmultislice3(tempnu, char(swmaparcfnamethres),'axial','nih.lut','0.5 1.5','-30 6 58');
    copyfile(tempbs_wm_slovname,path_qccompwm);
    
    % 2. Isolate Brainstem 
    
    bsaparcfname=strcat(appath,'/',apfname,'_bs.nii');
    
    spm('defaults','PET');
    matlabbatch{1}.spm.util.imcalc.input = cellstr(tempaparc);
    matlabbatch{1}.spm.util.imcalc.output = char(bsaparcfname);
    matlabbatch{1}.spm.util.imcalc.outdir = {''};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1==16';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 0;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch); clear matlabbatch;
    
    % 3. combine WM, bs and wcbl that has to be already there
    
    compwmfname=strcat(appath,'/compwm_ref_mask.nii');
    
    compwmfs=vertcat(swmaparcfnamethres,bsaparcfname,cellstr(strcat(appath,'/wholecbl_ref_mask.nii')));
    
    spm('defaults','PET');
    matlabbatch{1}.spm.util.imcalc.input = compwmfs;
    matlabbatch{1}.spm.util.imcalc.output = char(compwmfname);
    matlabbatch{1}.spm.util.imcalc.outdir = {''};
    matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3)>0';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 0;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch); clear matlabbatch;
    
    % Grab the rSCAN to get the scaling factor to be saved in the database
    
    [imgpath,imgfname,~]=spm_fileparts(char(tempimg));
    
    rimgfname=strcat(imgpath,'/r',imgfname(1:end-9),'.nii');
    
    % Extract value
    
    roi=spm_vol(char(strcat(appath,'/compwm_ref_mask.nii'))); 
    roi1=spm_read_vols(roi);
    
    img=spm_vol(rimgfname); % Reading img header from the loop
  img1=spm_read_vols(img); % Reading img values 
  img_mask=img1.*roi1; % Creating the masked image
  img_mask=nonzeros(img_mask);
  img_mask=img_mask(~isnan(img_mask));
  ext_val=mean(img_mask);
  
  M_longFBB(i,1)=ext_val;
  
  % Extract value from suvr image
  
  imgsuvr=spm_vol(char(tempimg)); % Reading imgsuvr header from the loop
  imgsuvr1=spm_read_vols(imgsuvr); % Reading imgsuvr values 
  imgsuvr_mask=imgsuvr1.*roi1; % Creating the masked image
  imgsuvr_mask=nonzeros(imgsuvr_mask);
  imgsuvr_mask=imgsuvr_mask(~isnan(imgsuvr_mask));
  ext_val_suvr=mean(imgsuvr_mask);
  M_suvr_longFBB(i,1)=ext_val_suvr;
  
  M_sz_longFBB(i,1)=size(img_mask,1);
  
  % Create the new suvr image, last thing I want
      
    exp=char(strcat('i1/',num2str(ext_val)));
    newfname=char(strcat(imgpath,'/',imgfname(1:end-9),'_suvr_compWM.nii'));

    spm('defaults','PET');
    matlabbatch{1}.spm.util.imcalc.input = cellstr(rimgfname);
    matlabbatch{1}.spm.util.imcalc.output = newfname;
    matlabbatch{1}.spm.util.imcalc.outdir = {''};
    matlabbatch{1}.spm.util.imcalc.expression = exp;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch); clear matlabbatch;

    clear apfname appath tempnu bsaparcfname compwmfname compwmfs exp ext_val ext_val_suvr img img1 img_mask imgfname imgpath newfname rimgfname roi roi1 swmaparcfname swmaparcfnamethres tempaparc tempimg wmaparcfname 
    
end % end for each new image to be processed longitudinally

T_longFBB=array2table(M_longFBB);
Tsuvr_longFBB=array2table(M_suvr_longFBB);
T_sz_longFBB=array2table(M_sz_longFBB);
T_longFBB.Properties.VariableNames={'ScalingFactor_CompWM'};
Tsuvr_longFBB.Properties.VariableNames={'SUVR_CompWM'};
T_sz_longFBB.Properties.VariableNames={'ScalingFactor_CompWM_ClustSize'};

templdsids_longFBB=regexp(fbbscans_longFBB_spm,'LDS\d{7}','match','once');
tempdates_longFBB=regexp(fbbscans_longFBB_spm,'\d{4}-\d{2}-\d{2}','match','once');
tempmridates_longFBB=regexp(aparcscans_longFBB_spm,'(?<=MRI_T1_)\d{4}-\d{2}-\d{2}','match','once');

meta_longFBB=horzcat(templdsids_longFBB, tempdates_longFBB, tempmridates_longFBB);
T_meta_longFBB=array2table(meta_longFBB);
T_meta_longFBB.Properties.VariableNames={'ID','FBBPET_Date','MRI_Date'};

qc_longFBB=horzcat(fbbscans_longFBB_spm, aparcscans_longFBB_spm);
T_qc_longFBB=array2table(qc_longFBB);
T_qc_longFBB.Properties.VariableNames={'FBBPET_path','APARC_path'};

T_longFBB=[T_meta_longFBB T_qc_longFBB T_longFBB Tsuvr_longFBB T_sz_longFBB]; 

if size(newcases_longFBB,1)>0
T_longFBB=vertcat(oldinfo_longFBB,T_longFBB);
end % end if condition new longitudinal cases exist

filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/FBB_CompWM_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
writetable(T_longFBB,filename,'WriteRowNames',true)
copyfile(filename,'/shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/');

end % end if condition existence of files in need of longitudinal extraction


% Done, wrapping up - clear everything before proceeding to the FTP ROI
% Extraction
            
clearvars -except tdb_ftp tdb_fdg path_processed path_extraction path_qccompwm
fprintf(1,'**Completed ROI and Longitudinal extraction for all the FBB SUVR images available!\n');
