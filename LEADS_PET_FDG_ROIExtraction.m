% Create a temp table with the info from the "updated" list of scans
% available in the processed folder
% Sep 2021: Code cleanup

% Feedback for the user

fprintf(1,'Checking what FDG scans are new...\n');

% Generate the entire list of images now available at all timepoints in the
% processed folder

newldsids=regexp(tdb_fdg(:,1),'LDS\d{7}','match','once');
newldsdates=regexp(tdb_fdg(:,1),'\d{4}-\d{2}-\d{2}','match','once');
newinfo=horzcat(newldsids,newldsdates);
newinfo=array2table(newinfo);
newinfo.Properties.VariableNames={'ID','FDGPET_Date'};

%% Read the last extraction database, this will be used to know what's new %%
%% This command reads all the files in the extraction folder and then reads the latest one %%
%% The final aim is to cleanup our list to only extract from new scans. If no extraction files are found then it will extract from all the scans

olddb=dir('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/FDG_ROI_Extraction*');

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

    checkinfo=isequal(newinfo(:,{'ID','FDGPET_Date'}),oldinfo(:,{'ID','FDGPET_Date'}));% Check if the b-ids columns in the new and old tables is unchanged

        if checkinfo==0

        Indexa=not(ismember(newinfo(:,{'ID','FDGPET_Date'}),oldinfo(:,{'ID','FDGPET_Date'}))); %% logical indexing for the new cases. Both are based on BID+date to include longitudinal scans
        Indexb=not(ismember(oldinfo(:,{'ID','FDGPET_Date'}),newinfo(:,{'ID','FDGPET_Date'}))); %% logical indexing for the removed cases

        newcases=newinfo(Indexa,:); %% Creates a table with the New FDG scans
        remcases=oldinfo(Indexb,:); %% Creates a table with eventually removed FDG scans

            if size(newcases,1)>0

                fprintf(2,'**New FDG-PET scans were found for Cross-Sectional Processing:\n\n');
                disp(newcases);
                tdb_fdg=tdb_fdg(Indexa,:);

            end % end if condition there are new scans available

            if size(remcases,1)>0

                fprintf(2,'Warning! %s cases were in the last FDG database and now are not found:\n\n', num2str(size(remcases,1)));

                choice = input(['\n\n Choose an action:' ,...
                    '\n     [1] Delete removed cases from present database',...
                    '\n     [2] Keep removed cases, I will check later'
                    '\n     --> ']);

                if choice==1

                    filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/Removed_Cases_FDG_ROI_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
                    writetable(remcases,filename,'WriteRowNames',true)
                    oldinfo(Indexb,:)=[];

                elseif choice==2

                        fprintf(1,'Fine! I will still save a csv file to keep track of that.\n\n');
                        filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/RemovedbutKept_Cases_FDG_ROI_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
                        writetable(remcases,filename,'WriteRowNames',true)
                end % end if condition choice on removed cases
                
                clear choice filename

            end % end if condition existence of removed cases
            
            elseif checkinfo==1
                
                fprintf(2,'No new FDG scans needed extraction!\n');
                tdb_fdg={};
                
         end % end if condition there is a difference between old extraction and new extraction list
        
    end % end if condition there is a previous FDG extraction database available
    
    if size(tdb_fdg,1)>0
    
% 1. Cross-sectional Module

%%% Create the list of regions we are interested in, plus creating labels

fsids=[2;  4;	5;	7;	8;	10;	11;	12;	13;	14;	15;	16;	17;	18;	24;	26;	28;	30;	31;	41;	43;	44;	46;	47;	49;	50;	51;	52;	53;	54;	58;	60;	62;	63;	72;	77;	80;	85;	251;	252;	253;	254;	255;	1000;	1001;	1002;	1003;	1005;	1006;	1007;	1008;	1009;	1010;	1011;	1012;	1013;	1014;	1015;	1016;	1017;	1018;	1019;	1020;	1021;	1022;	1023;	1024;	1025;	1026;	1027;	1028;	1029;	1030;	1031;	1032;	1033;	1034;	1035;	2000;	2001;	2002;	2003;	2005;	2006;	2007;	2008;	2009;	2010;	2011;	2012;	2013;	2014;	2015;	2016;	2017;	2018;	2019;	2020;	2021;	2022;	2023;	2024;	2025;	2026;	2027;	2028;	2029;	2030;	2031;	2032;	2033;	2034;	2035];
fslabs={'Left_Cerebral_White_Matter'; 'Left_Lateral_Ventricle';	'Left_Inf_Lat_Vent';	'Left_Cerebellum_White_Matter';	'Left_Cerebellum_Cortex';	'Left_Thalamus_Proper';	'Left_Caudate';	'Left_Putamen';	'Left_Pallidum';	'Third_Ventricle';	'Fourth_Ventricle';	'Brain_Stem';	'Left_Hippocampus';	'Left_Amygdala';	'CSF';	'Left_Accumbens_area';	'Left_VentralDC';	'Left_vessel';	'Left_choroid_plexus';	'Right_Cerebral_White_Matter'; 'Right_Lateral_Ventricle';	'Right_Inf_Lat_Vent';	'Right_Cerebellum_White_Matter';	'Right_Cerebellum_Cortex';	'Right_Thalamus_Proper';	'Right_Caudate';	'Right_Putamen';	'Right_Pallidum';	'Right_Hippocampus';	'Right_Amygdala';	'Right_Accumbens_area';	'Right_VentralDC';	'Right_vessel';	'Right_choroid_plexus';	'Fifth_Ventricle';	'WM_hypointensities';	'non_WM_hypointensities';	'Optic_Chiasm';	'CC_Posterior';	'CC_Mid_Posterior';	'CC_Central';	'CC_Mid_Anterior';	'CC_Anterior';	'ctx_lh_unknown';	'ctx_lh_bankssts';	'ctx_lh_caudalanteriorcingulate';	'ctx_lh_caudalmiddlefrontal';	'ctx_lh_cuneus';	'ctx_lh_entorhinal';	'ctx_lh_fusiform';	'ctx_lh_inferiorparietal';	'ctx_lh_inferiortemporal';	'ctx_lh_isthmuscingulate';	'ctx_lh_lateraloccipital';	'ctx_lh_lateralorbitofrontal';	'ctx_lh_lingual';	'ctx_lh_medialorbitofrontal';	'ctx_lh_middletemporal';	'ctx_lh_parahippocampal';	'ctx_lh_paracentral';	'ctx_lh_parsopercularis';	'ctx_lh_parsorbitalis';	'ctx_lh_parstriangularis';	'ctx_lh_pericalcarine';	'ctx_lh_postcentral';	'ctx_lh_posteriorcingulate';	'ctx_lh_precentral';	'ctx_lh_precuneus';	'ctx_lh_rostralanteriorcingulate';	'ctx_lh_rostralmiddlefrontal';	'ctx_lh_superiorfrontal';	'ctx_lh_superiorparietal';	'ctx_lh_superiortemporal';	'ctx_lh_supramarginal';	'ctx_lh_frontalpole';	'ctx_lh_temporalpole';	'ctx_lh_transversetemporal';	'ctx_lh_insula';	'ctx_rh_unknown';	'ctx_rh_bankssts';	'ctx_rh_caudalanteriorcingulate';	'ctx_rh_caudalmiddlefrontal';	'ctx_rh_cuneus';	'ctx_rh_entorhinal';	'ctx_rh_fusiform';	'ctx_rh_inferiorparietal';	'ctx_rh_inferiortemporal';	'ctx_rh_isthmuscingulate';	'ctx_rh_lateraloccipital';	'ctx_rh_lateralorbitofrontal';	'ctx_rh_lingual';	'ctx_rh_medialorbitofrontal';	'ctx_rh_middletemporal';	'ctx_rh_parahippocampal';	'ctx_rh_paracentral';	'ctx_rh_parsopercularis';	'ctx_rh_parsorbitalis';	'ctx_rh_parstriangularis';	'ctx_rh_pericalcarine';	'ctx_rh_postcentral';	'ctx_rh_posteriorcingulate';	'ctx_rh_precentral';	'ctx_rh_precuneus';	'ctx_rh_rostralanteriorcingulate';	'ctx_rh_rostralmiddlefrontal';	'ctx_rh_superiorfrontal';	'ctx_rh_superiorparietal';	'ctx_rh_superiortemporal';	'ctx_rh_supramarginal';	'ctx_rh_frontalpole';	'ctx_rh_temporalpole';	'ctx_rh_transversetemporal';	'ctx_rh_insula'};

% Generate list of images I will work with

vols=char(tdb_fdg(:,1));
vols_spm=tdb_fdg(:,1);
numv1= size(vols_spm,1);

aparcs = char(tdb_fdg(:,2));
aparc_spm = tdb_fdg(:,2);

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
     tempfdgdate=regexp(vols_spm(v),'\d{4}-\d{2}-\d{2}','match','once');
     tempfdgfold=regexp(vols_spm(v),'\FDG_\d{4}-\d{2}-\d{2}','match','once');
     temptimepoint=regexp(vols_spm(v),'\Timepoint\d{1}','match','once');

     refreg=strcat('/mnt/coredata/Projects/LEADS/data_f7p1/processed/',templdsid,'/',temptimepoint,'/',tempfdgfold,'/pons_ref_mask.nii');
     tempsuv=strcat('/mnt/coredata/Projects/LEADS/data_f7p1/processed/',templdsid,'/',temptimepoint,'/',tempfdgfold,'/r',templdsid,'_FDG_',tempfdgdate,'.nii');

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

     clear r1n r1 img img1 p f e templdsid tempftpdate tempftpfold temptimepoint refreg tempsuv refregn regreg1 suv suv1 mask suv_mask temp refregval refreg_sz

     end % end for loop for each image

T_refreg=array2table(M_refreg);
T_refreg_sz=array2table(M_refreg_sz);
T_refreg.Properties.VariableNames={'ScalingFactor_Pons'};
T_refreg_sz.Properties.VariableNames={'ScalingFactor_Pons_ClustSize'};

T = array2table(M); % create Tables for value reports
T_sz=array2table(M_sz);
T.Properties.VariableNames = fslabs;
T_sz.Properties.VariableNames = cellstr(strcat(fslabs,'_ClustSize'));

templdsids=regexp(vols_spm,'LDS\d{7}','match','once');
tempdates=regexp(vols_spm,'\d{4}-\d{2}-\d{2}','match','once');
tempmridates=regexp(aparc_spm,'(?<=MRI_T1_)\d{4}-\d{2}-\d{2}','match','once');

meta=horzcat(templdsids, tempdates, tempmridates);
T_meta=array2table(meta);
T_meta.Properties.VariableNames={'ID','FDGPET_Date','MRI_Date'};

qc=horzcat(vols_spm, aparc_spm);
T_qc=array2table(qc);
T_qc.Properties.VariableNames={'FDGPET_path','APARC_path'};

T=[T_meta T_qc T_refreg T T_refreg_sz T_sz]; 

    if size(newcases,1)>0
    T=vertcat(oldinfo,T);
    end % end if condition there are new cases

filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/FDG_ROI_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
writetable(T,filename,'WriteRowNames',true)
copyfile(filename,'/shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/');

    end % end if condition new FDG scans needed extractions
