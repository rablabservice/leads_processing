#!/usr/bin/env python

"""
Parse ROI extraction files in processed scan directories and compile
them in the style of LEADS ROI extraction CSV files (one for each
PET tracer: FBB, FTP, and FDG).
"""

import argparse
import os.path as op
import sys
import time
import warnings

import pandas as pd

utils_dir = op.join(op.dirname(__file__), "..", "utils")
if utils_dir not in sys.path:
    sys.path.append(utils_dir)
import utilities as uts

# Define globals
TODAY = pd.Timestamp.today().strftime("%Y-%m-%d")


class XReport:
    """Class to manage the creation of LEADS ROI extraction files."""

    def __init__(self, proj_dir="/mnt/coredata/processing/leads"):
        """Initialize the XReport object.

        Parameters
        ----------
        proj_dir : str
            Path to the top-level project directory
        report_date : pd.Timestamp | str like "YYYY-MM-DD"
            Date to use for the report period
        """
        self.get_paths(proj_dir)

        # Define each (tracer, ref. region) pair that we will save an
        # ROI extraction spreadsheet for
        self.extraction_keys = (
            ("fbb", "wcbl"),
            ("fbb", "compwm"),
            ("ftp", "infcblgm"),
            ("ftp", "eroded-subcortwm"),
            ("fdg", "pons"),
        )

        # Define amyloid tracers
        self.amyloid_tracers = ["fbb"]

    def __repr__(self):
        return f"XReport(proj_dir={self.paths['proj']})"

    def compile(self, save_output=True, overwrite=True):
        self.load_data()
        self.create_extraction_files(save_output, overwrite)
        print("-" * 80, "Finished compiling LEADS ROI extraction files!", sep="\n")

    def get_paths(self, proj_dir):
        """Get a dictionary of important project paths.

        Parameters
        ----------
        proj_dir : str
            Path to the top-level project directory

        Creates
        -------
        self.paths : dict
        """
        self.paths = {"proj": proj_dir}
        self.paths["code"] = op.join(self.paths["proj"], "code")
        self.paths["metadata"] = op.join(self.paths["proj"], "metadata")
        self.paths["atri"] = op.join(self.paths["metadata"], "atri")
        self.paths["loni"] = op.join(self.paths["metadata"], "loni")
        self.paths["qc"] = op.join(self.paths["metadata"], "qc")
        self.paths["scans_to_process"] = op.join(
            self.paths["metadata"], "scans_to_process"
        )
        self.paths["data"] = op.join(self.paths["proj"], "data")
        self.paths["extraction"] = op.join(self.paths["data"], "extraction")
        self.paths["qreport"] = op.join(
            self.paths["extraction"], "quarterly_report_files"
        )
        self.paths["rois"] = op.join(self.paths["extraction"], "internal_roi_files")
        self.paths["proc"] = op.join(self.paths["data"], "processed")

    def load_data(self):
        """Load all data needed for the XReport."""
        self.load_processed_pet_index()
        self.load_cohort()
        self.load_eligibility_list()
        self.load_subject_demo()
        self.load_apoe_genotype()
        self.load_anti_amyloid_treatment()
        self.load_clinical_baseline_data()
        self.load_ref_region_dat()
        self.load_roi_dat()
        self.load_centiloid_dat()
        self.load_umich_qc()
        self.load_ucsf_qc()

    def load_processed_pet_index(self):
        """Load the dataframe with all processed PET scans.

        This is the raw_PET_index*.csv file created by
        `select_scans_to_process.py` in the Setup Module of the
        processing pipeline.

        Creates
        -------
        self.tracers : list
            List of PET tracers that we are working with
        self.pet_idx : dict
            Dictionary of processed PET scan dataframes, one per tracer
        """
        # Load the latest PET index CSV file from scans_to_process
        keep_cols = [
            "subj",
            "tracer",
            "pet_date",
            "pet_image_id",
            "mri_date",
            "mri_image_id",
            "days_mri_to_pet",
            "pet_proc_dir",
        ]
        pet_scan_idx = pd.read_csv(
            uts.glob_sort_mtime(
                op.join(self.paths["scans_to_process"], "raw_PET_index*.csv")
            )[0]
        )

        # Remove PET scans that have not been fully processed
        pet_scan_idx = (
            pet_scan_idx.loc[pet_scan_idx["pet_processing_complete"] == 1, keep_cols]
            .reset_index(drop=True)
            .copy()
        )

        # Rename columns
        pet_scan_idx = pet_scan_idx.rename(columns={"subj": "subject_id"})

        # Modify columns
        pet_scan_idx["tracer"] = pet_scan_idx["tracer"].str.lower()

        # Separate each PET tracer into its own dataframe
        n_scans = len(pet_scan_idx)
        n_subjs = pet_scan_idx["subject_id"].nunique()
        print(
            "",
            "-" * 80,
            f"Loading latest PET scan index from {self.paths['scans_to_process']}",
            sep="\n",
        )
        print(
            f"- {n_scans:,} fully processed PET scans from {n_subjs:,} subjects",
            end="\n" * 2,
        )
        self.tracers = sorted(pet_scan_idx["tracer"].unique())
        self.pet_idx = {}
        for tracer, grp in pet_scan_idx.groupby("tracer"):
            self.pet_idx[tracer] = grp.reset_index(drop=True).copy()

    def set_pet_subjs(self):
        """Set the list of PET subjects from the processed PET index."""
        self.pet_subjs = set()
        for tracer in self.tracers:
            self.pet_subjs.update(self.pet_idx[tracer]["subject_id"].unique())

    def load_cohort(self):
        """Load the subject index with cohort assignment info.

        Creates
        -------
        self.cohort : DataFrame
            DataFrame with subject IDs and cohort assignments
        """
        self.set_pet_subjs()

        # Initialize the subjects dataframe
        self.cohort = pd.DataFrame(self.pet_subjs, columns=["subject_id"])

        # ------------------------------
        # Load the study_group dataframe
        study_group = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_subject.csv")
        )

        # Rename columns
        study_group = study_group.rename(
            columns={"subject_label": "subject_id", "ptcoh": "study_group"}
        )

        # Map integer values to string labels for the study_group column
        study_group["study_group"] = study_group["study_group"].map({1: "CN", 2: "PT"})

        # Select needed columns
        keep_cols = ["subject_id", "study_group"]
        study_group = study_group[keep_cols]

        # Merge into the main subject dataframe
        self.cohort = self.cohort.merge(
            study_group,
            on="subject_id",
            how="left",
        )

        # --------------------------------------
        # Load the amyloid eligibility dataframe
        amyelg = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_amyelg.csv")
        )

        # Rename columns
        amyelg = amyelg.rename(
            columns={
                "subject_label": "subject_id",
                "amyelg": "CohortAssgn",
            }
        )

        # Select rows from Screening visits
        amyelg = amyelg[amyelg["event_code"] == "sc"]

        # Select neeeded columns
        keep_cols = [
            "subject_id",
            "CohortAssgn",
        ]
        amyelg = amyelg[keep_cols]

        # Map integer values to string labels
        amyelg["CohortAssgn"] = amyelg["CohortAssgn"].map({0: "EOnonAD", 1: "EOAD"})

        # Merge into the main subject dataframe
        self.cohort = self.cohort.merge(
            amyelg,
            on="subject_id",
            how="left",
        )

        # Combine the two diagnosis columns
        self.cohort["CohortAssgn"] = self.cohort.apply(
            lambda x: (
                x["CohortAssgn"] if pd.notnull(x["CohortAssgn"]) else x["study_group"]
            ),
            axis=1,
        )
        self.cohort = self.cohort.drop(columns=["study_group"])

    def load_eligibility_list(self):
        """Load the LEADS eligibility list.

        Creates
        -------
        self.eligibility : DataFrame
        """
        eligibility = pd.read_csv(
            op.join(self.paths["atri"], "LEADS_Eligibility_list.csv")
        )
        self.eligible_subjects = set(eligibility["subject_label"])

    def load_subject_demo(self):
        """Load LEADS subject demographic data.

        Creates
        -------
        self.subj_demo : DataFrame
        """
        self.subj_demo = pd.read_csv(
            uts.glob_sort_mtime(op.join(self.paths["loni"], "LEADS_PTDEMOG*.csv"))[0]
        )

        # Rename columns
        self.subj_demo = self.subj_demo.rename(
            columns={
                "subject_code": "subject_id",
                "ptdob": "dob",
                "ptgender": "sex",
                "ptraccat": "race",
                "ptethcat": "ethnicity",
                "pteducat": "years_education",
            }
        )

        # Make subject IDs upper case
        self.subj_demo["subject_id"] = self.subj_demo["subject_id"].str.upper()

        # Map sex values
        self.subj_demo["sex"] = self.subj_demo["sex"].map({1: "Male", 2: "Female"})

        # Map race values
        self.subj_demo.loc[
            self.subj_demo["race"].astype(str).str.contains("\|"), "race"
        ] = 6
        self.subj_demo["race"] = self.subj_demo["race"].map(
            {
                1: "American Indian or Alaskan Native",
                2: "Asian",
                3: "Native Hawaiian or Other Pacific Islander",
                4: "Black or African American",
                5: "White",
                6: "More than one race",
                7: "Unknown",
            }
        )

        # Map ethnicity values
        self.subj_demo["ethnicity"] = self.subj_demo["ethnicity"].map(
            {1: "Hispanic or Latino", 2: "Not Hispanic or Latino", 3: "Unknown"}
        )

        # Select only needed columns
        keep_cols = ["subject_id", "dob", "sex", "race", "ethnicity", "years_education"]
        self.subj_demo = self.subj_demo[keep_cols]

    def load_apoe_genotype(self):
        """Load LEADS APOE genotype data.

        Creates
        -------
        self.apoe_geno : DataFrame
        """
        self.apoe_geno = pd.read_csv(
            uts.glob_sort_mtime(
                op.join(self.paths["loni"], "Biospecimen_Analysis_Results*.csv")
            )[0]
        )

        # Filter to only include APOE Genotype results
        self.apoe_geno = self.apoe_geno.query("TESTNAME == 'APOE Genotype'")

        # Rename columns
        self.apoe_geno = self.apoe_geno.rename(
            columns={"SUBJECT_CODE": "subject_id", "TESTVALUE": "apoe_genotype"}
        )

        # Add column counting number of APOE4 alleles
        self.apoe_geno["apoe4_alleles"] = self.apoe_geno["apoe_genotype"].apply(
            lambda x: str(x).count("4")
        )

        # Select only needed columns
        keep_cols = ["subject_id", "apoe_genotype", "apoe4_alleles"]
        self.apoe_geno = self.apoe_geno[keep_cols]

    def load_anti_amyloid_treatment(self):
        """Load LEADS anti-amyloid treatment data.

        Creates
        -------
        self.antiamy_tx : DataFrame
        """
        self.antiamy_tx = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_antiamytx.csv")
        )

        # Rename columns
        self.antiamy_tx = self.antiamy_tx.rename(
            columns={
                "subject_label": "subject_id",
                "txrec": "antiamy_tx_received",
            }
        )

        # Remap treatment received values from 1/2 to 1/0
        self.antiamy_tx["antiamy_tx_received"] = self.antiamy_tx[
            "antiamy_tx_received"
        ].map(
            {1: 1, 2: 0}  # 1 = Yes, 0 = No
        )

        # Record whether a patient has ever received anti-amyloid
        # treatment, according to any of their visit forms
        self.antiamy_tx = (
            self.antiamy_tx.groupby("subject_id")["antiamy_tx_received"]
            .any()
            .astype("Int64")
            .reset_index()
        )

    def load_clinical_baseline_data(self):
        """Load LEADS subject baseline clinical data.

        Calls subfunctions for individual assessment dataframes,
        all cognitive testing and clinical characteristics at baseline.
        """
        self.load_preliminary_dx()
        self.load_mmse_baseline()
        self.load_cdr_baseline()
        self.load_dx_pca()
        self.load_dx_lvppa()

    def load_preliminary_dx(self):
        """Load LEADS preliminary diagnosis data.

        Categorizes patients as MCI or Dementia at baseline.

        Creates
        -------
        self.prelim_dx : DataFrame
        """
        self.prelim_dx = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_prelimdx.csv")
        )

        # Filter to only include screening visits
        self.prelim_dx = self.prelim_dx.query("event_code == 'sc'")

        # Rename columns
        self.prelim_dx = self.prelim_dx.rename(
            columns={
                "subject_label": "subject_id",
                "dementia": "clinical_severity_baseline",
            }
        )

        # Map clinical severity values
        self.prelim_dx["clinical_severity_baseline"] = self.prelim_dx[
            "clinical_severity_baseline"
        ].map({0: "MCI", 1: "Dementia"})

        # Select only needed columns
        keep_cols = ["subject_id", "clinical_severity_baseline"]
        self.prelim_dx = self.prelim_dx[keep_cols]

    def load_mmse_baseline(self):
        """Load LEADS MMSE baseline data.

        Creates
        -------
        self.mmse_baseline : DataFrame
        """
        self.mmse_baseline = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_mmse.csv")
        )

        # Filter to only include screening visits
        self.mmse_baseline = self.mmse_baseline.query("event_code == 'sc'")

        # Rename columns
        self.mmse_baseline = self.mmse_baseline.rename(
            columns={
                "subject_label": "subject_id",
                "examdate": "mmse_date_baseline",
                "mmscore": "mmse_baseline",
            }
        )

        # Select only needed columns
        keep_cols = ["subject_id", "mmse_baseline"]
        self.mmse_baseline = self.mmse_baseline[keep_cols]

    def load_cdr_baseline(self):
        """Load LEADS CDR baseline data.

        Creates
        -------
        self.cdr_baseline : DataFrame
        """
        self.cdr_baseline = pd.read_csv(
            uts.glob_sort_mtime(
                op.join(self.paths["loni"], "Clinical_Dementia_Rating*.csv")
            )[0]
        )

        # Filter to only include screening visits
        self.cdr_baseline = self.cdr_baseline.query("LEADS_SCREENING_VISIT == 'sc'")

        # Create cdr_date_baseline column from year, month, day columns
        self.cdr_baseline["cdr_date_baseline"] = pd.to_datetime(
            {
                "year": self.cdr_baseline["C2VISITYR"],
                "month": self.cdr_baseline["C2VISITMO"],
                "day": self.cdr_baseline["C2VISITDAY"],
            }
        ).dt.strftime("%Y-%m-%d")

        # Rename columns
        self.cdr_baseline = self.cdr_baseline.rename(
            columns={
                "LEADS_ID": "subject_id",
                "CDRGLOB": "cdr_global_baseline",
                "CDRSUM": "cdr_sb_baseline",
            }
        )

        # Capitalize subject_id values
        self.cdr_baseline["subject_id"] = self.cdr_baseline["subject_id"].str.upper()

        # Select only needed columns
        keep_cols = [
            "subject_id",
            "cdr_date_baseline",
            "cdr_global_baseline",
            "cdr_sb_baseline",
        ]
        self.cdr_baseline = self.cdr_baseline[keep_cols]

    def load_dx_pca(self):
        """Load LEADS PCA diagnosis data.

        Creates
        -------
        self.dx_pca : DataFrame
        """
        self.dx_pca = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_pcadx.csv")
        )

        # Filter to only include screening visits
        self.dx_pca = self.dx_pca.query("event_code == 'sc'")

        # Rename columns
        self.dx_pca = self.dx_pca.rename(
            columns={"subject_label": "subject_id", "pcaformal": "pca_formal"}
        )

        # Select only needed columns
        keep_cols = ["subject_id", "pca_formal"]
        self.dx_pca = self.dx_pca[keep_cols]

    def load_dx_lvppa(self):
        """Load LEADS lvPPA diagnosis data.

        Creates
        -------
        self.dx_lvppa : DataFrame
        """
        self.dx_lvppa = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_ppadx.csv")
        )

        # Filter to only include screening visits
        self.dx_lvppa = self.dx_lvppa.query("event_code == 'sc'")

        # Rename columns
        self.dx_lvppa = self.dx_lvppa.rename(
            columns={"subject_label": "subject_id", "lvppaformal": "lvppa_formal"}
        )

        # Select only needed columns
        keep_cols = ["subject_id", "lvppa_formal"]
        self.dx_lvppa = self.dx_lvppa[keep_cols]

    def get_roi_name_from_file(self, filepath, prefix):
        """Return the ROI name from a filepath."""
        roi_name = op.basename(filepath)
        i = roi_name.find(prefix)
        if i != -1:
            roi_name = roi_name[i + len(prefix) :]
        roi_name = ".".join(roi_name.split(".")[:-1])
        return roi_name

    def load_ref_region_dat(self):
        """Load reference region scaling factors for each tracer.

        Creates
        -------
        self.ref_region_dat : dict
            Dictionary with one dataframe of reference region scaling
            factors for each tracer
        """

        def get_ref_region_means(pet_proc_dir):
            def find_ref_region_means_file(pet_proc_dir):
                """Return the filepath to the ROI mean CSV file."""
                subj, tracer, pet_date = uts.parse_scan_tag(
                    uts.get_scan_tag(pet_proc_dir)
                )
                filepath = op.join(
                    pet_proc_dir,
                    f"r{subj}_{tracer}_{pet_date}_ref-region-means.csv",
                )
                if op.isfile(filepath):
                    return filepath
                else:
                    warnings.warn(f"File not found: {filepath}")

            def load_ref_region_means_file(filepath):
                """Load and format the ROI extractions CSV file."""

                def scrape_ref_regions(mask_file):
                    ref_regions = []
                    for roi in mask_file.split(";"):
                        ref_regions.append(self.get_roi_name_from_file(roi, "mask-"))
                    if len(ref_regions) > 1:
                        return "compwm"
                    else:
                        return ref_regions[0]

                # Load the CSV file
                df = pd.read_csv(filepath)

                # Parse the PET proc dir
                subj, _, pet_date = uts.parse_scan_tag(uts.get_scan_tag(filepath))

                # Add subject_id and pet_date columns
                df.insert(0, "subject_id", subj)
                df.insert(1, "pet_date", pet_date)
                df.insert(
                    2,
                    "ref_region",
                    df["mask_file"].apply(scrape_ref_regions),
                )

                return df

            def format_ref_region_means(df):
                """Format ROI extractions dataframe for LEADS ROI extractions."""
                # Format ROI names
                df["ref_region"] = df["ref_region"].astype(
                    pd.CategoricalDtype(df["ref_region"], ordered=True)
                )

                # Remove unnecessary columns
                df = df.drop(columns=["image_file", "mask_file", "voxel_count"])

                # Rename columns
                df = df.rename(
                    columns={
                        "mean": "ScalingFactor",
                    }
                )

                # Pivot the dataframe from long to wide format
                df = df.set_index(["subject_id", "pet_date", "ref_region"]).unstack(
                    "ref_region"
                )

                # Flatten the column index
                df.columns = ["_".join(col).strip() for col in df.columns.values]

                # Reset the index
                df = df.reset_index()

                return df

            try:
                rr_dat = format_ref_region_means(
                    load_ref_region_means_file(find_ref_region_means_file(pet_proc_dir))
                )
                return rr_dat
            except Exception as e:
                print(e)
                return None

        # Load reference region scaling factors for each tracer
        self.ref_region_dat = {}
        for tracer in self.tracers:
            self.ref_region_dat[tracer] = pd.concat(
                list(
                    self.pet_idx[tracer].apply(
                        lambda x: get_ref_region_means(x["pet_proc_dir"]), axis=1
                    )
                ),
                ignore_index=True,
            )

    def load_roi_dat(self):
        """Load ROI extraction means for each tracer.

        Creates
        -------
        self.roi_dat : dict
            Dictionary of dataframes, one per (tracer, ref. region) pair
        """

        def get_roi_extractions(pet_proc_dir, ref_region):
            def find_roi_extractions_file(pet_proc_dir, ref_region):
                """Return the filepath to the ROI mean CSV file."""
                subj, tracer, pet_date = uts.parse_scan_tag(
                    uts.get_scan_tag(pet_proc_dir)
                )
                filepath = op.join(
                    pet_proc_dir,
                    f"r{subj}_{tracer}_{pet_date}_suvr-{ref_region}_roi-extractions.csv",
                )
                if op.isfile(filepath):
                    return filepath
                else:
                    warnings.warn(f"File not found: {filepath}")

            def load_roi_extractions_file(filepath):
                """Load and format the ROI extractions CSV file."""
                # Load the CSV file
                df = pd.read_csv(filepath)

                # Add subject_id and pet_date columns
                df.insert(
                    0,
                    "subject_id",
                    df["image_file"].apply(lambda x: op.basename(x).split("_")[0][1:]),
                )
                df.insert(
                    1,
                    "pet_date",
                    df["image_file"].apply(lambda x: op.basename(x).split("_")[2]),
                )

                return df

            def format_roi_extractions(df):
                """Format ROI extractions dataframe for LEADS ROI extractions."""
                # Format ROI names
                df["roi"] = df["roi"].astype(
                    pd.CategoricalDtype(df["roi"], ordered=True)
                )

                # Remove unnecessary columns
                df = df.drop(columns=["image_file", "roi_file"])

                # Rename columns
                df = df.rename(
                    columns={
                        "mean": "MRIBASED_SUVR",
                        "voxel_count": "ClustSize",
                    }
                )

                # Pivot the dataframe from long to wide format
                df = df.set_index(["subject_id", "pet_date", "roi"]).unstack("roi")

                # Flatten the column index
                df.columns = ["_".join(col[::-1]).strip() for col in df.columns.values]

                # Reset the index
                df = df.reset_index()

                return df

            try:
                roi_dat = format_roi_extractions(
                    load_roi_extractions_file(
                        find_roi_extractions_file(pet_proc_dir, ref_region)
                    )
                )
                return roi_dat
            except Exception as e:
                print(e)
                return None

        # Load ROI extractions for each tracer
        self.roi_dat = {}
        for key in self.extraction_keys:
            tracer, ref_region = key
            self.roi_dat[key] = pd.concat(
                list(
                    self.pet_idx[tracer].apply(
                        lambda x: get_roi_extractions(x["pet_proc_dir"], ref_region),
                        axis=1,
                    )
                ),
                ignore_index=True,
            )

    def load_centiloid_dat(self):
        """Load Centiloid values for amyloid PET.

        Creates
        -------
        self.centiloid_dat : dict
            Dictionary with one dataframe of Centiloid values per tracer
        """

        def get_centiloids(pet_proc_dir):
            def find_centiloid_file(pet_proc_dir):
                """Return the filepath to the ROI mean CSV file."""
                subj, tracer, pet_date = uts.parse_scan_tag(
                    uts.get_scan_tag(pet_proc_dir)
                )
                filepath = op.join(
                    pet_proc_dir,
                    f"r{subj}_{tracer}_{pet_date}_amyloid-cortical-summary.csv",
                )
                if op.isfile(filepath):
                    return filepath
                else:
                    warnings.warn(f"File not found: {filepath}")

            def load_centiloid_file(filepath):
                """Load and format the ROI extractions CSV file."""

                def _scrape_ref_region(image_file):
                    ref_region = (
                        op.basename(image_file)
                        .split(".")[0]
                        .split("_")[3]
                        .split("-")[1]
                    )
                    return ref_region

                # Load the CSV file
                df = pd.read_csv(filepath)

                # Parse the PET proc dir
                subj, _, pet_date = uts.parse_scan_tag(
                    uts.get_scan_tag(op.dirname(filepath))
                )

                # Add subject_id and pet_date columns
                df.insert(0, "subject_id", subj)
                df.insert(1, "pet_date", pet_date)
                df.insert(
                    2,
                    "ref_region",
                    df["image_file"].apply(
                        lambda x: self.get_roi_name_from_file(x, "suvr-")
                    ),
                )

                return df

            def format_centiloid_dat(df):
                """Format ROI extractions dataframe for LEADS ROI extractions."""
                # Format ROI names
                df["ref_region"] = df["ref_region"].astype(
                    pd.CategoricalDtype(df["ref_region"], ordered=True)
                )

                # Remove unnecessary columns
                df = df.drop(columns=["image_file", "mask_file", "mean_suvr"])

                # Pivot the dataframe from long to wide format
                df = df.set_index(["subject_id", "pet_date", "ref_region"]).unstack(
                    "ref_region"
                )

                # Flatten the column index
                df.columns = ["_".join(col) for col in df.columns.values]

                # Reset the index
                df = df.reset_index()

                return df

            try:
                rr_dat = format_centiloid_dat(
                    load_centiloid_file(find_centiloid_file(pet_proc_dir))
                )
                return rr_dat
            except Exception as e:
                print(e)
                return None

        # Load Centiloid values
        self.centiloid_dat = {}
        for tracer in self.amyloid_tracers:
            self.centiloid_dat[tracer] = pd.concat(
                list(
                    self.pet_idx[tracer].apply(
                        lambda x: get_centiloids(x["pet_proc_dir"]), axis=1
                    )
                ),
                ignore_index=True,
            )

    def load_umich_qc(self):
        """Filter PET scans and retain rows from scans that passed UMich QC.

        Creates
        -------
        self.qc_mich : dict
            Dictionary of dataframes, one per tracer
        """

        # Load the UMich QC spreadsheets
        self.qc_mich = {
            "fbb": pd.read_csv(
                op.join(self.paths["atri"], "leads_codebook_study_data_amyqc.csv")
            ),
            "ftp": pd.read_csv(
                op.join(self.paths["atri"], "leads_codebook_study_data_tauqc.csv")
            ),
            "fdg": pd.read_csv(
                op.join(self.paths["atri"], "leads_codebook_study_data_fdgqc.csv")
            ),
        }

        # Rename columns
        for tracer in self.tracers:
            if tracer == "ftp":
                self.qc_mich[tracer] = self.qc_mich[tracer].rename(
                    columns={"scanqlty": "scanqltya"}
                )
            self.qc_mich[tracer] = self.qc_mich[tracer].rename(
                columns={
                    "subject_label": "subject_id",
                    "scandate": "pet_date",
                    "scanqltya": "qc_umich_passed",
                }
            )

        # Retain only the columns we need
        keep_cols = ["subject_id", "pet_date", "qc_umich_passed"]
        for tracer in self.tracers:
            self.qc_mich[tracer] = self.qc_mich[tracer][keep_cols]

        # Fix scan date mismatches between our processing and UMich CSVs
        idx = (
            self.qc_mich["ftp"]
            .query("(subject_id=='LDS0360283') & (pet_date=='2020-11-11')")
            .index
        )
        self.qc_mich["ftp"].loc[idx, "pet_date"] = "2020-11-12"
        idx = (
            self.qc_mich["ftp"]
            .query("(subject_id=='LDS0370672') & (pet_date=='2023-09-13')")
            .index
        )
        self.qc_mich["ftp"].loc[idx, "pet_date"] = "2023-09-14"
        idx = (
            self.qc_mich["fdg"]
            .query("(subject_id=='LDS0370012') & (pet_date=='2020-12-17')")
            .index
        )
        self.qc_mich["fdg"].loc[idx, "pet_date"] = "2020-12-15"

    def load_ucsf_qc(self):
        """Filter PET scans and retain rows from scans that passed UCSF QC.

        Creates
        -------
        self.qc_ucsf : dict
            Dictionary of dataframes, one per tracer
        """
        # Load UCSF QC spreadsheets
        self.qc_ucsf = {
            "mri": pd.read_csv(
                uts.glob_sort(
                    op.join(self.paths["qc"], "processed_MRI-T1_qc-evals_*.csv")
                )[-1]
            )
        }
        for tracer in self.tracers:
            self.qc_ucsf[tracer] = pd.read_csv(
                uts.glob_sort(
                    op.join(
                        self.paths["qc"], f"processed_{tracer.upper()}_qc-evals_*.csv"
                    )
                )[-1]
            )

        # Rename columns
        for scan_type in self.qc_ucsf:
            if scan_type == "mri":
                self.qc_ucsf[scan_type] = self.qc_ucsf[scan_type].rename(
                    columns={
                        "subj": "subject_id",
                        "scan_date": "mri_date",
                    }
                )
            else:
                self.qc_ucsf[scan_type] = self.qc_ucsf[scan_type].rename(
                    columns={
                        "subj": "subject_id",
                        "scan_date": "pet_date",
                    }
                )

            # Prefix all columns except subject_id and scan_date with "qc_ucsf_"
            cols = self.qc_ucsf[scan_type].columns
            self.qc_ucsf[scan_type].columns = [
                f"qc_ucsf_{x}" if x not in ["subject_id", "mri_date", "pet_date"] else x
                for x in cols
            ]

        # Drop unnecessary columns
        for scan_type in self.qc_ucsf:
            if scan_type == "mri":
                drop_cols = [
                    "qc_ucsf_rater",
                    "qc_ucsf_spm_seg_ok",
                    "qc_ucsf_affine_nu_ok",
                    "qc_ucsf_warped_nu_ok",
                    "qc_ucsf_notes",
                ]
                self.qc_ucsf[scan_type] = self.qc_ucsf[scan_type].drop(
                    columns=drop_cols
                )
            else:
                drop_cols = [
                    "qc_ucsf_rater",
                    "qc_ucsf_affine_pet_ok",
                    "qc_ucsf_warped_pet_ok",
                    "qc_ucsf_notes",
                ]
                self.qc_ucsf[scan_type] = self.qc_ucsf[scan_type].drop(
                    columns=drop_cols
                )

    def create_extraction_files(self, save_output=True, overwrite=True):
        """Merge spreadsheets and save each extraction file."""
        self.extraction = {}

        # Prepare the final dataframes for each tracer
        for key in self.extraction_keys:
            tracer, ref_region = key
            print("-" * 80)
            print(
                f"Building ROI extraction spreadsheet for {tracer.upper()} - {ref_region.upper()}"
            )

            # Construct the main dataframe for each extraction file.
            # Start by copying the processed scan list
            self.extraction[key] = self.pet_idx[tracer].copy()

            # Sort the data by subject_id and PET date
            self.extraction[key] = (
                self.extraction[key]
                .sort_values(["subject_id", "pet_date"])
                .reset_index(drop=True)
            )

            # Add cohort assignment
            self.extraction[key] = self.extraction[key].merge(
                self.cohort,
                on="subject_id",
                how="left",
            )

            # Add a column for whether subject is in the eligibility list
            self.extraction[key]["eligible_subject"] = (
                self.extraction[key]["subject_id"]
                .isin(self.eligible_subjects)
                .astype("Int64")
            )

            # Merge subject demographic and clinical data into the main dataframe.
            merge_dfs = [
                self.subj_demo,
                self.antiamy_tx,
                self.apoe_geno,
                self.prelim_dx,
                self.mmse_baseline,
                self.cdr_baseline,
                self.dx_pca,
                self.dx_lvppa,
            ]
            for df in merge_dfs:
                self.extraction[key] = self.extraction[key].merge(
                    df, on="subject_id", how="left"
                )

            # Insert age at PET in place of date of birth
            mask = (
                self.extraction[key]["dob"].notna()
                & self.extraction[key]["pet_date"].notna()
            )
            self.extraction[key].insert(
                self.extraction[key].columns.tolist().index("dob"), "age_at_pet", pd.NA
            )
            self.extraction[key].loc[mask, "age_at_pet"] = round(
                (
                    pd.to_datetime(self.extraction[key].loc[mask, "pet_date"])
                    - pd.to_datetime(self.extraction[key].loc[mask, "dob"])
                ).dt.days
                / 365.25,
                1,
            )
            self.extraction[key] = self.extraction[key].drop(columns=["dob"])

            # Merge reference region scaling factors into the main dataframe
            keep_cols = ["subject_id", "pet_date", f"ScalingFactor_{ref_region}"]
            self.extraction[key] = self.extraction[key].merge(
                self.ref_region_dat[tracer][keep_cols],
                on=["subject_id", "pet_date"],
                how="left",
            )

            # Merge Centiloid data into the main dataframe
            keep_cols = ["subject_id", "pet_date", f"centiloids_{ref_region}"]
            if tracer in self.amyloid_tracers:
                self.extraction[key] = self.extraction[key].merge(
                    self.centiloid_dat[tracer][keep_cols],
                    on=["subject_id", "pet_date"],
                    how="left",
                )

            # Merge ROI mean SUVR and volume data into the main dataframe
            self.extraction[key] = self.extraction[key].merge(
                self.roi_dat[key], on=["subject_id", "pet_date"], how="left"
            )

            # Merge UMich QC columns into the main dataframe
            self.extraction[key] = self.extraction[key].merge(
                self.qc_mich[tracer],
                on=["subject_id", "pet_date"],
                how="left",
            )

            # Merge UCSF QC columns into the main dataframe
            self.extraction[key] = self.extraction[key].merge(
                self.qc_ucsf[tracer],
                on=["subject_id", "pet_date"],
                how="left",
            )

            # Save the output dataframe
            output_file = op.join(
                self.paths["rois"],
                f"LEADS_{tracer.upper()}-ROI-means_{TODAY}.csv",
            )
            if save_output:
                if overwrite or not op.isfile(output_file):
                    self.extraction[key].to_csv(output_file, index=False)
                    print(f"Saved {output_file}")

            # Print final dataframe shape
            n_scans = len(self.extraction[key])
            n_subjs = self.extraction[key]["subject_id"].nunique()
            print(
                "",
                f"Final {tracer.upper()} ROI extraction dataframe:",
                f"- {n_scans:,} {tracer.upper()} scans from {n_subjs:,} subjects",
                f"- Dataframe shape: {self.extraction[key].shape}",
                sep="\n",
                end="\n" * 2,
            )


def _parse_args():
    """Parse and return command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Save LEADS PET ROI extraction files\n\n"
            + "Overview\n--------\n"
            + "This program aggregates ROI extraction files in LEADS processed PET scan\n"
            + "directories, creating 3 CSV files (one each for FBB, FTP, and FDG) that serve\n"
            + "as internal reference spreadsheets containing all standard ROI extraction data\n"
            + "for all processed PET scans, along with PET and MRI scan QC results and basic\n"
            + "subject-level demographic and clinical information."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        exit_on_error=False,
    )
    parser.add_argument(
        "-p",
        "--proj-dir",
        default="/mnt/coredata/processing/leads",
        help=(
            "Full path to the top-level project directory, from which subdirectory\n"
            "paths are inferred. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--no-save",
        action="store_false",
        dest="save",
        help="Generate ROI extraction dataframes but do not save them as CSVs",
    )
    parser.add_argument(
        "--no-overwrite",
        action="store_false",
        dest="overwrite",
        help="Do not overwrite existing ROI extraction CSV files",
    )

    return parser.parse_args()


def main():
    """Main function to run the XReport class."""
    start_time = time.time()

    # Get command line arguments.
    args = _parse_args()

    # Instantiate the XReport class and run the main function.
    xr = XReport(args.proj_dir)
    xr.compile(args.save, args.overwrite)

    # Report elapsed time.
    elapsed = time.time() - start_time
    minutes, seconds = divmod(elapsed, 60)
    print(f"Elapsed time: {int(minutes)} min, {int(seconds)} s", end="\n" * 2)

    # Exit successfully
    sys.exit(0)


if __name__ == "__main__":
    main()
