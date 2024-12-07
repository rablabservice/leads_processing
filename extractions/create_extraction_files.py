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

import numpy as np
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

    def __repr__(self):
        return (
            f"XReport(proj_dir={self.paths['proj']}, report_date={self.report_period})"
        )

    def compile(self, save_output=True, overwrite=True):
        self.load_processed_pet_index()
        self.filter_pet_by_umich_qc()
        self.filter_pet_by_ucsf_qc()
        self.load_subject_index()
        self.load_eligibility_list()
        self.load_ref_region_dat()
        self.load_roi_dat()
        self.load_centiloid_dat()
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

    def load_processed_pet_index(self):
        """Load the dataframe with all processed PET scans.

        This is the raw_PET_index*.csv file created by
        `select_scans_to_process.py` in the Setup Module of the
        processing pipeline.

        Creates
        -------
        self.tracers : list
            List of PET tracers that are working with
        self.pet_idx : dict
            Dictionary of processed PET scan dataframes, one per tracer
        """
        # Load the latest PET index CSV file from scans_to_process
        keep_cols = [
            "subj",
            "tracer",
            "pet_date",
            "pet_image_id",
            "pet_scan_number",
            "n_pet_scans",
            "days_from_baseline_pet",
            "days_from_last_pet",
            "pet_res",
            "mri_date",
            "mri_image_id",
            "days_mri_to_pet",
            "abs_days_mri_to_pet",
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
        print(f"- {n_scans:,} fully processed PET scans from {n_subjs:,} subjects")
        self.tracers = sorted(pet_scan_idx["tracer"].unique())
        self.pet_idx = {}
        for tracer, grp in pet_scan_idx.groupby("tracer"):
            n_scans = len(grp)
            n_subjs = grp["subject_id"].nunique()
            self.pet_idx[tracer] = grp.reset_index(drop=True).copy()
            print(f"  * {n_scans:>5,} {tracer.upper()}, {n_subjs:>3,} subjects")
        print()

    def filter_pet_by_umich_qc(self, drop_failed_scans=False):
        """Filter PET scans and retain rows from scans that passed UMich QC.

        Parameters
        ----------
        drop_failed_scans : bool
            Whether to drop scans that failed UMich QC.
            - If False, all scans in the input dataframes are retained,
            and a Boolean "passed_umich_qc" column is added to the output
            dataframes.
            - If True, only scans that passed QC are retained, and the
            "passed_umich_qc" column is dropped from the output
            dataframes.

        Updates
        -------
        self.pet_idx : dict
            Dictionary of processed PET scan dataframes, one per tracer
        """

        # Load the UMich QC spreadsheets
        qc_mich = {
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
        qc_mich["fbb"] = qc_mich["fbb"].rename(
            columns={
                "subject_label": "subject_id",
                "scandate": "pet_date",
                "scanqltya": "passed_umich_qc",
            }
        )
        qc_mich["ftp"] = qc_mich["ftp"].rename(
            columns={
                "subject_label": "subject_id",
                "scandate": "pet_date",
                "scanqlty": "passed_umich_qc",
            }
        )
        qc_mich["fdg"] = qc_mich["fdg"].rename(
            columns={
                "subject_label": "subject_id",
                "scandate": "pet_date",
                "scanqltya": "passed_umich_qc",
            }
        )

        # Retain only the columns we need
        keep_cols = ["subject_id", "pet_date", "passed_umich_qc"]
        for tracer in self.tracers:
            qc_mich[tracer] = qc_mich[tracer][keep_cols]

        # Select rows that passed QC, making sure each scan has only one row
        for tracer in self.tracers:
            qc_mich[tracer] = qc_mich[tracer].query("(passed_umich_qc==1)")
            assert len(qc_mich[tracer]) == len(
                qc_mich[tracer].drop_duplicates(["subject_id", "pet_date"])
            )

        # Fix scan date mismatches between our processing and UMich CSVs
        idx = (
            qc_mich["ftp"]
            .query("(subject_id=='LDS0360283') & (pet_date=='2020-11-11')")
            .index
        )
        qc_mich["ftp"].loc[idx, "pet_date"] = "2020-11-12"
        idx = (
            qc_mich["ftp"]
            .query("(subject_id=='LDS0370672') & (pet_date=='2023-09-13')")
            .index
        )
        qc_mich["ftp"].loc[idx, "pet_date"] = "2023-09-14"
        idx = (
            qc_mich["fdg"]
            .query("(subject_id=='LDS0370012') & (pet_date=='2020-12-17')")
            .index
        )
        qc_mich["fdg"].loc[idx, "pet_date"] = "2020-12-15"

        # Merge into the main PET scans dataframe
        for tracer in self.tracers:
            self.pet_idx[tracer] = self.pet_idx[tracer].merge(
                qc_mich[tracer],
                on=["subject_id", "pet_date"],
                how="left",
            )

        # Select scans that passed QC
        print("-" * 80, "Checking PET scans against UMich QC", sep="\n")
        drop_msg = "Removing" if drop_failed_scans else "These are the"
        for tracer in self.tracers:
            n_scans = len(self.pet_idx[tracer])
            n_passed = int(self.pet_idx[tracer]["passed_umich_qc"].sum())
            print(f"- {n_passed:,}/{n_scans:,} {tracer.upper()} scans passed UMich QC")
            if n_scans > n_passed:
                print(
                    f"- {drop_msg} {n_scans - n_passed:,} {tracer.upper()} scans that did not pass UMich QC:"
                )
                print(
                    "    "
                    + "\n    ".join(
                        self.pet_idx[tracer]
                        .query("(passed_umich_qc!=1)")[["subject_id", "pet_date"]]
                        .to_markdown(index=False, tablefmt="rst")
                        .split("\n")
                    )
                )
                print()
                if drop_failed_scans:
                    self.pet_idx[tracer] = (
                        self.pet_idx[tracer]
                        .query("(passed_umich_qc==1)")
                        .reset_index(drop=True)
                    )
                    self.pet_idx[tracer] = self.pet_idx[tracer].drop(
                        columns=["passed_umich_qc"]
                    )
        print()

    def filter_pet_by_ucsf_qc(self, drop_failed_scans=False):
        """Filter PET scans and retain rows from scans that passed UCSF QC.

        Parameters
        ----------
        drop_failed_scans : bool
            Whether to drop scans that failed UCSF QC.
            - If False, all scans in the input dataframes are retained,
            and a Boolean "ucsf_qc_pass" column is added to the output
            dataframes.
            - If True, only scans that passed QC are retained, and the
            "ucsf_qc_pass" column is dropped from the output
            dataframes.

        Updates
        -------
        self.pet_idx : dict
            Dictionary of processed PET scan dataframes, one per tracer
        """
        # Load UCSF QC spreadsheets
        qc_ucsf = {
            "mri": pd.read_csv(
                uts.glob_sort(
                    op.join(self.paths["qc"], "processed_MRI-T1_qc-evals_*.csv")
                )[-1]
            ),
            "fbb": pd.read_csv(
                uts.glob_sort(
                    op.join(self.paths["qc"], "processed_FBB_qc-evals_*.csv")
                )[-1]
            ),
            "fdg": pd.read_csv(
                uts.glob_sort(
                    op.join(self.paths["qc"], "processed_FDG_qc-evals_*.csv")
                )[-1]
            ),
            "ftp": pd.read_csv(
                uts.glob_sort(
                    op.join(self.paths["qc"], "processed_FTP_qc-evals_*.csv")
                )[-1]
            ),
        }

        # Rename columns
        for scan_type in qc_ucsf:
            if scan_type == "mri":
                qc_ucsf[scan_type] = qc_ucsf[scan_type].rename(
                    columns={
                        "subj": "subject_id",
                        "scan_date": "mri_date",
                        "notes": "mri_qc_notes",
                    }
                )
            else:
                qc_ucsf[scan_type] = qc_ucsf[scan_type].rename(
                    columns={
                        "subj": "subject_id",
                        "scan_date": "pet_date",
                        "notes": "pet_qc_notes",
                    }
                )

        # Drop unnecessary columns
        for scan_type in qc_ucsf:
            if scan_type == "mri":
                drop_cols = [
                    "rater",
                    "manual_freesurfer_edits",
                    "spm_seg_ok",
                    "affine_nu_ok",
                ]
                qc_ucsf[scan_type] = qc_ucsf[scan_type].drop(columns=drop_cols)
            else:
                drop_cols = [
                    "rater",
                    "affine_pet_ok",
                ]
                qc_ucsf[scan_type] = qc_ucsf[scan_type].drop(columns=drop_cols)

        # Add QC pass columns
        scan_type = "mri"
        qc_ucsf[scan_type].insert(
            2,
            "mri_qc_pass",
            qc_ucsf[scan_type].apply(
                lambda x: np.all(
                    (x["native_nu_rating"] > 0, x["aparc_rating"] > 0)
                ).astype(float),
                axis=1,
            ),
        )

        scan_type = "fbb"
        qc_ucsf[scan_type].insert(
            2,
            "pet_qc_pass",
            qc_ucsf[scan_type].apply(
                lambda x: np.all(
                    (
                        x["native_pet_ok"] > 0,
                        x["pet_to_mri_coreg_ok"] > 0,
                        x["wcbl_mask_ok"] > 0,
                    )
                ).astype(float),
                axis=1,
            ),
        )

        scan_type = "ftp"
        qc_ucsf[scan_type].insert(
            2,
            "pet_qc_pass",
            qc_ucsf[scan_type].apply(
                lambda x: np.all(
                    (
                        x["native_pet_ok"] > 0,
                        x["pet_to_mri_coreg_ok"] > 0,
                        x["infcblgm_mask_ok"] > 0,
                    )
                ).astype(float),
                axis=1,
            ),
        )

        scan_type = "fdg"
        qc_ucsf[scan_type].insert(
            2,
            "pet_qc_pass",
            qc_ucsf[scan_type].apply(
                lambda x: np.all(
                    (
                        x["native_pet_ok"] > 0,
                        x["pet_to_mri_coreg_ok"] > 0,
                        x["pons_mask_ok"] > 0,
                    )
                ).astype(float),
                axis=1,
            ),
        )

        # Merge the UCSF QC data into the PET scans dataframes
        for tracer in self.tracers:
            self.pet_idx[tracer] = self.pet_idx[tracer].merge(
                qc_ucsf[tracer],
                on=["subject_id", "pet_date"],
                how="left",
            )
            self.pet_idx[tracer] = self.pet_idx[tracer].merge(
                qc_ucsf["mri"],
                on=["subject_id", "mri_date"],
                how="left",
            )

        # Add a ucsf_qc_pass column combining MRI and PET QC fields
        for tracer in self.tracers:
            self.pet_idx[tracer]["ucsf_qc_pass"] = self.pet_idx[tracer].apply(
                lambda x: np.all(
                    (
                        x["mri_qc_pass"] > 0,
                        x["pet_qc_pass"] > 0,
                    )
                ).astype(float),
                axis=1,
            )

        # Select scans that passed QC
        print("-" * 80, "Checking PET scans against UCSF PET/MRI QC", sep="\n")
        drop_msg = "Removing" if drop_failed_scans else "These are the"
        for tracer in self.tracers:
            n_scans = len(self.pet_idx[tracer])
            n_passed = int(self.pet_idx[tracer]["ucsf_qc_pass"].sum())
            print(f"- {n_passed:,}/{n_scans:,} {tracer.upper()} scans passed UCSF QC")
            if n_scans > n_passed:
                print(
                    f"- {drop_msg} {n_scans - n_passed:,} {tracer.upper()} scans that did not pass UCSF QC:"
                )
                print(
                    "    "
                    + "\n    ".join(
                        self.pet_idx[tracer]
                        .query("(ucsf_qc_pass!=1)")[["subject_id", "pet_date"]]
                        .to_markdown(index=False, tablefmt="rst")
                        .split("\n")
                    )
                )
                print()
                if drop_failed_scans:
                    self.pet_idx[tracer] = (
                        self.pet_idx[tracer]
                        .query("(ucsf_qc_pass==1)")
                        .reset_index(drop=True)
                    )
                    self.pet_idx[tracer] = self.pet_idx[tracer].drop(
                        columns=["ucsf_qc_pass"]
                    )
        print()

    def set_pet_subjs(self):
        """Set the list of PET subjects from the processed PET index."""
        self.pet_subjs = set()
        for tracer in self.tracers:
            self.pet_subjs.update(self.pet_idx[tracer]["subject_id"].unique())

    # ### Cohort assignment
    # Determine each participant's enrollment status as a control or patient.
    # For patients, determine their diagnostic assignment into EOAD or EOnonAD
    # groups based on visual read of the baseline amyloid PET scan.
    #
    # - `atri/leads_codebook_study_data_subject.csv`
    # - Join on columns
    #     - `subject_label` -> `subject_id`
    # - Include columns
    #     - `ptcoh` -> `study_cohort`
    #     - 1 = "Control" (cognitively normal)
    #     - 2 = "Patient" (cognitively impaired)
    #
    # - `atri/leads_codebook_study_data_amyelg.csv`
    # - Filter by columns
    #     - `event_code` == "sc"
    # - Join on columns
    #     - `subject_label` -> `subject_id`
    # - Include columns
    #     - `amyelg` -> `CohortAssgn`
    #         - 0 = "EOnonAD"
    #         - 1 = "EOAD"
    #     - `suvr` -> `Screening_PETONLY_Composite_SUVR`
    #     - `outcome` -> `Screening_PETONLY_VisualRead`
    #     - `consensres` -> `Screening_PETONLY_Final_Read`
    # - Add columns
    #     - `Screening_PETONLY_AmyPos_Quantification_1p18`
    #     - `Screening_PETONLY_Disagreement`
    def load_subject_index(self):
        """Load the subject index with cohort assignment info.

        Creates
        -------
        self.subj_idx : DataFrame
            DataFrame with subject IDs and cohort assignments
        self.screening : DataFrame
            Same as self.subj_idx but with amyloid visual read columns
            from the screening visit retained
        """
        self.set_pet_subjs()

        # Initialize the subjects dataframe
        self.subj_idx = pd.DataFrame(self.pet_subjs, columns=["subject_id"])

        # Load the cohort dataframe
        cohort = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_subject.csv")
        )

        # Rename columns
        cohort = cohort.rename(
            columns={"subject_label": "subject_id", "ptcoh": "study_group"}
        )

        # Map integer values to string labels for the study_group column
        cohort["study_group"] = cohort["study_group"].map({1: "CN", 2: "PT"})

        # Select needed columns
        keep_cols = ["subject_id", "study_group"]
        cohort = cohort[keep_cols]

        # Merge into the main subject dataframe
        self.subj_idx = self.subj_idx.merge(
            cohort,
            on="subject_id",
            how="left",
        )

        # Load the amyloid eligibility dataframe
        amyelg = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_amyelg.csv")
        )

        # Rename columns
        amyelg = amyelg.rename(
            columns={
                "subject_label": "subject_id",
                "suvr": "Screening_PETONLY_Composite_SUVR",
                "outcome": "Screening_PETONLY_VisualRead",
                "consensres": "Screening_PETONLY_Final_Read",
                "amyelg": "CohortAssgn",
            }
        )

        # Select rows from Screening visits
        amyelg = amyelg[amyelg["event_code"] == "sc"]

        # Add additional columns
        FBB_PETONLY_THRESH = 1.18
        amyelg["Screening_PETONLY_AmyPos_Quantification_1p18"] = (
            amyelg["Screening_PETONLY_Composite_SUVR"] > FBB_PETONLY_THRESH
        ).astype(int)
        amyelg["Screening_PETONLY_Disagreement"] = np.logical_xor(
            amyelg["Screening_PETONLY_AmyPos_Quantification_1p18"],
            amyelg["Screening_PETONLY_VisualRead"],
        ).astype(int)

        # Replace nans in Screening_PETONLY_Final_Read with
        # Screening_PETONLY_VisualRead values
        amyelg["Screening_PETONLY_Final_Read"] = amyelg.apply(
            lambda x: (
                x["Screening_PETONLY_Final_Read"]
                if pd.notnull(x["Screening_PETONLY_Final_Read"])
                else x["Screening_PETONLY_VisualRead"]
            ),
            axis=1,
        )

        # Map integer values to string labels
        amyelg["CohortAssgn"] = amyelg["CohortAssgn"].map({0: "EOnonAD", 1: "EOAD"})

        # Select neeeded columns
        keep_cols = [
            "subject_id",
            "Screening_PETONLY_Composite_SUVR",
            "Screening_PETONLY_AmyPos_Quantification_1p18",
            "Screening_PETONLY_VisualRead",
            "Screening_PETONLY_Disagreement",
            "Screening_PETONLY_Final_Read",
            "CohortAssgn",
        ]
        amyelg = amyelg[keep_cols]

        # Merge into the main subject dataframe
        self.subj_idx = self.subj_idx.merge(
            amyelg,
            on="subject_id",
            how="left",
        )

        # Combine the two diagnosis columns
        self.subj_idx["CohortAssgn"] = self.subj_idx.apply(
            lambda x: (
                x["CohortAssgn"] if pd.notnull(x["CohortAssgn"]) else x["study_group"]
            ),
            axis=1,
        )
        self.subj_idx = self.subj_idx.drop(columns=["study_group"])

        # Create the screening dataframe and drop unnecessary columns from the
        # subjects dataframe
        self.screening = self.subj_idx.copy()
        keep_cols = ["subject_id", "CohortAssgn"]
        self.subj_idx = self.subj_idx[keep_cols]

    # ### Participant demographics
    # Determine each participant's date of birth, sex, race and ethnicity, and years of education
    #
    # - `loni/LEADS_PTDEMOG*.csv`
    # - Join on columns
    #     - `subject_code` -> `subject_id`
    # - Include columns
    #     - `ptdob` -> `dob`
    #     - `ptgender` -> `sex`
    #     - 1 = "Male"
    #     - 2 = "Female"
    #     - `ptraccat` -> `race`
    #     - 1 = "American Indian or Alaskan Native"
    #     - 2 = "Asian"
    #     - 3 = "Native Hawaiian or Other Pacific Islander"
    #     - 4 = "Black or African American"
    #     - 5 = "White"
    #     - 6 = "More than one race"
    #     - 7 = "Unknown"
    #     - `ptethcat` -> `ethnicity`
    #     - 1 = "Hispanic or Latino"
    #     - 2 = "Not Hispanic or Latino"
    #     - 3 = "Unknown"
    #     - `pteducat` -> `years_education`
    #     - Years of education
    def load_demographics(self):
        """Load LEADS subject demographics data.

        Creates
        -------
        self.pt_demo : DataFrame
        """
        self.pt_demographic_dat = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_subject.csv")
        )
        # TODO - WRITE THE REST OF THIS FUNCTION

    # ### Baseline clinical characteristics
    # Determine each patient's clinical severity (MCI or Dementia), cognitive
    # test scores (MMSE, CDR), and clinical phenotype at baseline.
    #
    # - `atri/leads_codebook_study_data_prelimdx.csv`
    # - Filter by columns
    #     - `event_code` == "sc"
    # - Join on columns
    #     - `subject_label` -> `subject_id`
    # - Include columns
    #     - `dementia` -> `clinical_severity_baseline`
    #     - 0 = "MCI"
    #     - 1 = "Dementia"
    #
    # - `atri/leads_codebook_study_data_mmse.csv`
    # - Filter by columns
    #     - `event_code` == "sc"
    # - Join on columns
    #     - `subject_label` -> `subject_id`
    # - Include columns
    #     - `mmscore` -> `mmse_baseline`
    #
    # - `loni/Clinical_Dementia_Rating*.csv`
    # - Filter by columns
    #     - `LEADS_SCREENING_VISIT` == "sc"
    # - Join on columns
    #     - `LEADS_ID` -> `subject_id`
    # - Include columns
    #     - `C2VISITYR`
    #     - `C2VISITMO`
    #     - `C2VISITDAY`
    #     - `CDRGLOB` -> `cdr_global_baseline`
    #     - `CDRSUM` -> `cdr_sb_baseline`
    # - Add columns
    #     - `cdr_date`
    #
    # - `atri/leads_codebook_study_data_pcadx.csv`
    # - Filter by columns
    #     - `event_code` == "sc"
    # - Join on columns
    #     - `subject_label` -> `subject_id`
    # - Include columns
    #     - `pcaformal` -> `pca_formal`
    #         - 0 = "Does not meet formal clinical diagnostic criteria for PCA"
    #         - 1 = "Meets formal clinical diagnostic criteria for PCA"
    #
    # - `atri/leads_codebook_study_data_ppadx.csv`
    # - Filter by columns
    #     - `event_code` == "sc"
    # - Join on columns
    #     - `subject_label` -> `subject_id`
    # - Include columns
    #     - `lvppaformal` -> `lvppa_formal`
    #         - 0 = "Does not meet formal clinical diagnostic criteria for lvPPA"
    #         - 1 = "Meets formal clinical diagnostic criteria for lvPPA"
    def load_baseline_clinical(self):
        """Load LEADS baseline clinical data.

        Creates
        -------
        self.clinical_dat : DataFrame
        """
        self.clinical_dat = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_prelimdx.csv")
        )
        # TODO - WRITE THE REST OF THIS FUNCTION

    # ### APOE genotype
    # Determine each participant's APOE genotype
    #
    # - `loni/Biospecimen_Analysis_Results*.csv`
    # - Join on columns
    #     - `SUBJECT_CODE` -> `subject_id`
    # - Filter on columns
    #     - `TESTNAME` == "APOE Genotype"
    # - Include columns
    #     - `TESTVALUE` -> `apoe_genotype`
    # - Add columns
    #     - `apoe4_alleles`
    #     - 0
    #     - 1
    #     - 2
    def load_apoe_genotype(self):
        """Load LEADS APOE genotype data.

        Creates
        -------
        self.apoe_dat : DataFrame
        """
        self.apoe_dat = pd.read_csv(
            op.join(self.paths["loni"], "Biospecimen_Analysis_Results*.csv")
        )
        # TODO - WRITE THE REST OF THIS FUNCTION

    ### Anti-amyloid treatment
    def load_anti_amyloid_treatment(self):
        """Load LEADS anti-amyloid treatment data.

        Creates
        -------
        self.treatment_dat : DataFrame
        """
        self.treatment_dat = pd.read_csv(
            op.join(self.paths["atri"], "leads_codebook_study_data_antiamy.csv")
        )
        # TODO - WRITE THE REST OF THIS FUNCTION

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
                    ref_regions = [
                        op.basename(x).split(".")[0].split("_")[3].split("-")[1]
                        for x in mask_file.split(";")
                    ]
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
                df["ref_region"] = df["ref_region"].apply(lambda x: x.replace("-", "_"))
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
            Dictionary with one dataframe of reference region scaling
            factors for each tracer
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
                df["roi"] = df["roi"].apply(lambda x: x.replace("-", "_"))
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

        # Define reference regions for each tracer
        ref_regions = {
            "fbb": "wcbl",
            "ftp": "infcblgm",
            "fdg": "pons",
        }

        # Load ROI extractions for each tracer
        self.roi_dat = {}
        for tracer, ref_region in ref_regions.items():
            self.roi_dat[tracer] = pd.concat(
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
        self.centiloid_dat : DataFrame
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
                    df["image_file"].apply(_scrape_ref_region),
                )

                return df

            def format_centiloid_dat(df):
                """Format ROI extractions dataframe for LEADS ROI extractions."""
                # Format ROI names
                df["ref_region"] = df["ref_region"].apply(lambda x: x.replace("-", "_"))
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

        # Load Centiloid values for each reference region
        tracer = "fbb"
        self.centiloid_dat = pd.concat(
            list(
                self.pet_idx[tracer].apply(
                    lambda x: get_centiloids(x["pet_proc_dir"]), axis=1
                )
            ),
            ignore_index=True,
        )

    def create_extraction_files(self, save_output=True, overwrite=True):

        # Prepare the final dataframes for each tracer
        for tracer in self.tracers:
            print("-" * 80)
            print(f"Building ROI extraction spreadsheet for {tracer.upper()}")

            # Merge reference region scaling factors into the main dataframe
            self.pet_idx[tracer] = self.pet_idx[tracer].merge(
                self.ref_region_dat[tracer], on=["subject_id", "pet_date"], how="left"
            )
            if tracer == "fbb":
                # Merge amyloid visual read info into the main dataframe
                self.pet_idx[tracer] = self.pet_idx[tracer].merge(
                    self.screening, on="subject_id", how="left"
                )
                # Merge Centiloid data into the main dataframe
                self.pet_idx[tracer] = self.pet_idx[tracer].merge(
                    self.centiloid_dat, on=["subject_id", "pet_date"], how="left"
                )
            else:
                # Merge cohort assignment into the main dataframe
                self.pet_idx[tracer] = self.pet_idx[tracer].merge(
                    self.screening[["subject_id", "CohortAssgn"]],
                    on="subject_id",
                    how="left",
                )

            # Merge ROI mean SUVR and volume data into the main dataframe
            self.pet_idx[tracer] = self.pet_idx[tracer].merge(
                self.roi_dat[tracer], on=["subject_id", "pet_date"], how="left"
            )

            # Sort the data by subject_id and PET date
            self.pet_idx[tracer] = (
                self.pet_idx[tracer]
                .sort_values(["subject_id", "pet_date"])
                .reset_index(drop=True)
            )

            # Add a column for whether subject is in the eligibility list
            print("Checking eligibility list")
            self.pet_idx[tracer]["eligible_subject"] = (
                self.pet_idx[tracer]["subject_id"]
                .isin(self.eligible_subjects)
                .astype(int)
            )

            print(
                "{:,}/{:,} ({:.1%}) {} scans are on the eligibility list".format(
                    self.pet_idx[tracer]["eligible_subject"].sum(),
                    len(self.pet_idx[tracer]),
                    self.pet_idx[tracer]["eligible_subject"].mean(),
                    tracer.upper(),
                ),
                end="\n" * 2,
            )

            # Save the output dataframe
            output_file = op.join(
                self.paths["rois"],
                f"LEADS_{tracer.upper()}-ROI-means_{TODAY}.csv",
            )
            if save_output:
                if overwrite or not op.isfile(output_file):
                    self.pet_idx[tracer].to_csv(output_file, index=False)
                    print(f"Saved {output_file}")

            # Print final dataframe shape
            n_scans = len(self.pet_idx[tracer])
            n_subjs = self.pet_idx[tracer]["ID"].nunique()
            print(
                "",
                f"Final {tracer.upper()} ROI extraction dataframe:",
                f"- {n_scans:,} {tracer.upper()} scans from {n_subjs:,} subjects",
                f"- Dataframe shape: {self.pet_idx[tracer].shape}",
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


if __name__ == "__main__":
    start_time = time.time()

    # Get command line arguments.
    args = _parse_args()

    # Instantiate the QRReport class and run the main function.
    qr = XReport(args.proj_dir)
    qr.compile(args.save, args.overwrite)

    # Report elapsed time.
    elapsed = time.time() - start_time
    minutes, seconds = divmod(elapsed, 60)
    print(f"Elapsed time: {int(minutes)} min, {int(seconds)} s", end="\n" * 2)

    # Exit successfully
    sys.exit(0)
