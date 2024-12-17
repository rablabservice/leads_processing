#!/usr/bin/env python

"""
Parse ROI extraction files in processed scan directories and compile
them in the style of LEADS quarterly report CSV files (one for each
PET tracer: FBB, FTP, and FDG).
"""

import argparse
import os.path as op
import re
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


class QReport:
    """Class to manage the creation of LEADS quarterly report files."""

    def __init__(self, proj_dir="/mnt/coredata/processing/leads", report_date=None):
        """Initialize the QReport object.

        Parameters
        ----------
        proj_dir : str
            Path to the top-level project directory
        report_date : pd.Timestamp | str like "YYYY-MM-DD"
            Date to use for the report period
        """
        self.get_paths(proj_dir)
        self.get_report_period(report_date)
        self.get_stop_date()

    def __repr__(self):
        return (
            f"QReport(proj_dir={self.paths['proj']}, report_date={self.report_period})"
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
        self.create_qreport_files(save_output, overwrite)
        print("-" * 80, "Finished compiling LEADS quarterly report files!", sep="\n")

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

    def get_report_period(self, report_date=None):
        """Get the report period (e.g. '2024Q3') for a given date.

        The report period is determined by the month of the input date
        and will use today's date by default. If report_date is a valid
        report period string (e.g. '2024Q3'), it will be used as is.

        Parameters
        ----------
        report_date : pd.Timestamp | str like "YYYY-MM-DD" | None | "YYYYQ[1-4]"

        Quarters
        --------
        - Q1: January 1 - March 31
        - Q2: April 1 - June 30
        - Q3: July 1 - September 30
        - Q4: October 1 - December 31

        Creates
        -------
        self.report_period : str
        """
        if isinstance(report_date, str) and re.match(
            r"^(20[0-9]{2})-Q[1-4]$", report_date
        ):
            self.report_period = report_date
        else:
            report_date = (
                pd.Timestamp.today()
                if report_date is None
                else pd.Timestamp(report_date)
            )
            self.report_period = f"{report_date.year}-Q{report_date.quarter}"

    def get_stop_date(self):
        """Get the last date for inclusion in the current report period.

        LEADS PET data are reported at a one-quarter lag, so the stop date
        should match the end of the previous quarter.

        Creates
        -------
        self.stop_date : pd.Timestamp
        """
        this_year, this_quarter = self.report_period.split("-")
        this_year = int(this_year)
        this_quarter = int(this_quarter[1])

        if this_quarter == 1:
            stop_year = this_year - 1
            stop_month = 12
            stop_day = 31
        elif this_quarter == 2:
            stop_year = this_year
            stop_month = 3
            stop_day = 31
        elif this_quarter == 3:
            stop_year = this_year
            stop_month = 6
            stop_day = 30
        elif this_quarter == 4:
            stop_year = this_year
            stop_month = 9
            stop_day = 30

        self.stop_date = pd.Timestamp(year=stop_year, month=stop_month, day=stop_day)

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

    def filter_pet_by_umich_qc(self, drop_failed_scans=True):
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

    def filter_pet_by_ucsf_qc(self, drop_failed_scans=True):
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
                """Format ROI extractions dataframe for LEADS quarterly reports."""
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
                """Format ROI extractions dataframe for LEADS quarterly reports."""
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
                """Format ROI extractions dataframe for LEADS quarterly reports."""
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

    def create_qreport_files(self, save_output=True, overwrite=True):
        def format_tracer(self, tracer):
            def format_fbb(self):
                # Drop unnecessary columns
                drop_cols = [
                    "tracer",
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
                    "pet_qc_pass",
                    "native_pet_ok",
                    "pet_to_mri_coreg_ok",
                    "wcbl_mask_ok",
                    "eroded_wm_and_brainstem_masks_ok",
                    "warped_pet_ok",
                    "pet_qc_notes",
                    "mri_qc_notes",
                    "centiloids_compwm",
                    "brainstem_MRIBASED_SUVR",
                    "eroded_subcortwm_MRIBASED_SUVR",
                    "ctx_desikan_MRIBASED_SUVR",
                    "brainstem_ClustSize",
                    "eroded_subcortwm_ClustSize",
                    "amyloid_cortical_summary_ClustSize",
                    "ctx_desikan_ClustSize",
                ]
                self.qrep_dat[tracer] = self.qrep_dat[tracer].drop(columns=drop_cols)

                # Rename columns
                self.qrep_dat[tracer] = self.qrep_dat[tracer].rename(
                    columns={
                        "subject_id": "ID",
                        "pet_date": "FBBPET_Date",
                        "pet_image_id": "ImageID",
                        "ScalingFactor_wcbl": "ScalingFactor_WholeCereb",
                        "ScalingFactor_compwm": "ScalingFactor_CompositeWM",
                        "centiloids_wcbl": "MRIBASED_Composite_Centiloids",
                        "wcbl_MRIBASED_SUVR": "WholeCerebellum_MRIBASED_SUVR",
                        "amyloid_cortical_summary_MRIBASED_SUVR": "MRIBASED_Composite_SUVR",
                        "wcbl_ClustSize": "WholeCerebellum_ClustSize",
                    }
                )

                # Add empty columns
                add_cols = [
                    "ACC_PCC_MRIBASED_SUVR",
                    "Frontal_MRIBASED_SUVR",
                    "Temporal_MRIBASED_SUVR",
                    "Parietal_MRIBASED_SUVR",
                    "CompositeWM_MRIBASED_SUVR",
                    "Left_vessel_MRIBASED_SUVR",
                    "Right_vessel_MRIBASED_SUVR",
                    "Fifth_Ventricle_MRIBASED_SUVR",
                    "non_WM_hypointensities_MRIBASED_SUVR",
                    "ctx_lh_unknown_MRIBASED_SUVR",
                    "ctx_rh_unknown_MRIBASED_SUVR",
                    "ScalingFactor_WholeCereb_ClustSize",
                    "ScalingFactor_CompositeWM_ClustSize",
                    "ACC_PCC_ClustSize",
                    "Frontal_ClustSize",
                    "Temporal_ClustSize",
                    "Parietal_ClustSize",
                    "CompositeWM_ClustSize",
                    "Left_vessel_ClustSize",
                    "Right_vessel_ClustSize",
                    "Fifth_Ventricle_ClustSize",
                    "non_WM_hypointensities_ClustSize",
                    "ctx_lh_unknown_ClustSize",
                    "ctx_rh_unknown_ClustSize",
                ]
                for col in add_cols:
                    self.qrep_dat[tracer][col] = np.nan

                # Convert the scan date column to a string formatted like '%m/%d/%y'
                self.qrep_dat[tracer]["FBBPET_Date"] = self.qrep_dat[tracer][
                    "FBBPET_Date"
                ].apply(lambda x: pd.to_datetime(x).strftime("%m/%d/%y"))

                # Convert binary screening columsn to string (0 becomes "0" and 1 becomes "1")
                fmt_cols = [
                    "Screening_PETONLY_AmyPos_Quantification_1p18",
                    "Screening_PETONLY_VisualRead",
                    "Screening_PETONLY_Disagreement",
                    "Screening_PETONLY_Final_Read",
                ]
                for col in fmt_cols:
                    self.qrep_dat[tracer][col] = self.qrep_dat[tracer][col].apply(
                        lambda x: "" if pd.isna(x) else str(int(x))
                    )

                # Round scaling factor columns to 4 decimals
                round_cols = [
                    col
                    for col in self.qrep_dat[tracer].columns
                    if col.startswith("ScalingFactor")
                ]
                self.qrep_dat[tracer][round_cols] = self.qrep_dat[tracer][
                    round_cols
                ].round(4)

                # Round MRIBASED_SUVR columns to 7 decimals
                round_cols = [
                    col
                    for col in self.qrep_dat[tracer].columns
                    if col.endswith("MRIBASED_SUVR")
                ]
                round_cols += [
                    "MRIBASED_Composite_SUVR",
                    "MRIBASED_Composite_Centiloids",
                ]
                self.qrep_dat[tracer][round_cols] = self.qrep_dat[tracer][
                    round_cols
                ].round(7)

                # Reorder columns
                cols_in_order = [
                    "ID",
                    "FBBPET_Date",
                    "ImageID",
                    "Screening_PETONLY_Composite_SUVR",
                    "Screening_PETONLY_AmyPos_Quantification_1p18",
                    "Screening_PETONLY_VisualRead",
                    "Screening_PETONLY_Disagreement",
                    "Screening_PETONLY_Final_Read",
                    "CohortAssgn",
                    "ScalingFactor_WholeCereb",
                    "ScalingFactor_CompositeWM",
                    "MRIBASED_Composite_SUVR",
                    "MRIBASED_Composite_Centiloids",
                    "ACC_PCC_MRIBASED_SUVR",
                    "Frontal_MRIBASED_SUVR",
                    "Temporal_MRIBASED_SUVR",
                    "Parietal_MRIBASED_SUVR",
                    "WholeCerebellum_MRIBASED_SUVR",
                    "CompositeWM_MRIBASED_SUVR",
                    "Left_Cerebral_White_Matter_MRIBASED_SUVR",
                    "Left_Lateral_Ventricle_MRIBASED_SUVR",
                    "Left_Inf_Lat_Vent_MRIBASED_SUVR",
                    "Left_Cerebellum_White_Matter_MRIBASED_SUVR",
                    "Left_Cerebellum_Cortex_MRIBASED_SUVR",
                    "Left_Thalamus_Proper_MRIBASED_SUVR",
                    "Left_Caudate_MRIBASED_SUVR",
                    "Left_Putamen_MRIBASED_SUVR",
                    "Left_Pallidum_MRIBASED_SUVR",
                    "Third_Ventricle_MRIBASED_SUVR",
                    "Fourth_Ventricle_MRIBASED_SUVR",
                    "Brain_Stem_MRIBASED_SUVR",
                    "Left_Hippocampus_MRIBASED_SUVR",
                    "Left_Amygdala_MRIBASED_SUVR",
                    "CSF_MRIBASED_SUVR",
                    "Left_Accumbens_area_MRIBASED_SUVR",
                    "Left_VentralDC_MRIBASED_SUVR",
                    "Left_vessel_MRIBASED_SUVR",
                    "Left_choroid_plexus_MRIBASED_SUVR",
                    "Right_Cerebral_White_Matter_MRIBASED_SUVR",
                    "Right_Lateral_Ventricle_MRIBASED_SUVR",
                    "Right_Inf_Lat_Vent_MRIBASED_SUVR",
                    "Right_Cerebellum_White_Matter_MRIBASED_SUVR",
                    "Right_Cerebellum_Cortex_MRIBASED_SUVR",
                    "Right_Thalamus_Proper_MRIBASED_SUVR",
                    "Right_Caudate_MRIBASED_SUVR",
                    "Right_Putamen_MRIBASED_SUVR",
                    "Right_Pallidum_MRIBASED_SUVR",
                    "Right_Hippocampus_MRIBASED_SUVR",
                    "Right_Amygdala_MRIBASED_SUVR",
                    "Right_Accumbens_area_MRIBASED_SUVR",
                    "Right_VentralDC_MRIBASED_SUVR",
                    "Right_vessel_MRIBASED_SUVR",
                    "Right_choroid_plexus_MRIBASED_SUVR",
                    "Fifth_Ventricle_MRIBASED_SUVR",
                    "WM_hypointensities_MRIBASED_SUVR",
                    "non_WM_hypointensities_MRIBASED_SUVR",
                    "Optic_Chiasm_MRIBASED_SUVR",
                    "CC_Posterior_MRIBASED_SUVR",
                    "CC_Mid_Posterior_MRIBASED_SUVR",
                    "CC_Central_MRIBASED_SUVR",
                    "CC_Mid_Anterior_MRIBASED_SUVR",
                    "CC_Anterior_MRIBASED_SUVR",
                    "ctx_lh_unknown_MRIBASED_SUVR",
                    "ctx_lh_bankssts_MRIBASED_SUVR",
                    "ctx_lh_caudalanteriorcingulate_MRIBASED_SUVR",
                    "ctx_lh_caudalmiddlefrontal_MRIBASED_SUVR",
                    "ctx_lh_cuneus_MRIBASED_SUVR",
                    "ctx_lh_entorhinal_MRIBASED_SUVR",
                    "ctx_lh_fusiform_MRIBASED_SUVR",
                    "ctx_lh_inferiorparietal_MRIBASED_SUVR",
                    "ctx_lh_inferiortemporal_MRIBASED_SUVR",
                    "ctx_lh_isthmuscingulate_MRIBASED_SUVR",
                    "ctx_lh_lateraloccipital_MRIBASED_SUVR",
                    "ctx_lh_lateralorbitofrontal_MRIBASED_SUVR",
                    "ctx_lh_lingual_MRIBASED_SUVR",
                    "ctx_lh_medialorbitofrontal_MRIBASED_SUVR",
                    "ctx_lh_middletemporal_MRIBASED_SUVR",
                    "ctx_lh_parahippocampal_MRIBASED_SUVR",
                    "ctx_lh_paracentral_MRIBASED_SUVR",
                    "ctx_lh_parsopercularis_MRIBASED_SUVR",
                    "ctx_lh_parsorbitalis_MRIBASED_SUVR",
                    "ctx_lh_parstriangularis_MRIBASED_SUVR",
                    "ctx_lh_pericalcarine_MRIBASED_SUVR",
                    "ctx_lh_postcentral_MRIBASED_SUVR",
                    "ctx_lh_posteriorcingulate_MRIBASED_SUVR",
                    "ctx_lh_precentral_MRIBASED_SUVR",
                    "ctx_lh_precuneus_MRIBASED_SUVR",
                    "ctx_lh_rostralanteriorcingulate_MRIBASED_SUVR",
                    "ctx_lh_rostralmiddlefrontal_MRIBASED_SUVR",
                    "ctx_lh_superiorfrontal_MRIBASED_SUVR",
                    "ctx_lh_superiorparietal_MRIBASED_SUVR",
                    "ctx_lh_superiortemporal_MRIBASED_SUVR",
                    "ctx_lh_supramarginal_MRIBASED_SUVR",
                    "ctx_lh_frontalpole_MRIBASED_SUVR",
                    "ctx_lh_temporalpole_MRIBASED_SUVR",
                    "ctx_lh_transversetemporal_MRIBASED_SUVR",
                    "ctx_lh_insula_MRIBASED_SUVR",
                    "ctx_rh_unknown_MRIBASED_SUVR",
                    "ctx_rh_bankssts_MRIBASED_SUVR",
                    "ctx_rh_caudalanteriorcingulate_MRIBASED_SUVR",
                    "ctx_rh_caudalmiddlefrontal_MRIBASED_SUVR",
                    "ctx_rh_cuneus_MRIBASED_SUVR",
                    "ctx_rh_entorhinal_MRIBASED_SUVR",
                    "ctx_rh_fusiform_MRIBASED_SUVR",
                    "ctx_rh_inferiorparietal_MRIBASED_SUVR",
                    "ctx_rh_inferiortemporal_MRIBASED_SUVR",
                    "ctx_rh_isthmuscingulate_MRIBASED_SUVR",
                    "ctx_rh_lateraloccipital_MRIBASED_SUVR",
                    "ctx_rh_lateralorbitofrontal_MRIBASED_SUVR",
                    "ctx_rh_lingual_MRIBASED_SUVR",
                    "ctx_rh_medialorbitofrontal_MRIBASED_SUVR",
                    "ctx_rh_middletemporal_MRIBASED_SUVR",
                    "ctx_rh_parahippocampal_MRIBASED_SUVR",
                    "ctx_rh_paracentral_MRIBASED_SUVR",
                    "ctx_rh_parsopercularis_MRIBASED_SUVR",
                    "ctx_rh_parsorbitalis_MRIBASED_SUVR",
                    "ctx_rh_parstriangularis_MRIBASED_SUVR",
                    "ctx_rh_pericalcarine_MRIBASED_SUVR",
                    "ctx_rh_postcentral_MRIBASED_SUVR",
                    "ctx_rh_posteriorcingulate_MRIBASED_SUVR",
                    "ctx_rh_precentral_MRIBASED_SUVR",
                    "ctx_rh_precuneus_MRIBASED_SUVR",
                    "ctx_rh_rostralanteriorcingulate_MRIBASED_SUVR",
                    "ctx_rh_rostralmiddlefrontal_MRIBASED_SUVR",
                    "ctx_rh_superiorfrontal_MRIBASED_SUVR",
                    "ctx_rh_superiorparietal_MRIBASED_SUVR",
                    "ctx_rh_superiortemporal_MRIBASED_SUVR",
                    "ctx_rh_supramarginal_MRIBASED_SUVR",
                    "ctx_rh_frontalpole_MRIBASED_SUVR",
                    "ctx_rh_temporalpole_MRIBASED_SUVR",
                    "ctx_rh_transversetemporal_MRIBASED_SUVR",
                    "ctx_rh_insula_MRIBASED_SUVR",
                    "ScalingFactor_WholeCereb_ClustSize",
                    "ScalingFactor_CompositeWM_ClustSize",
                    "ACC_PCC_ClustSize",
                    "Frontal_ClustSize",
                    "Temporal_ClustSize",
                    "Parietal_ClustSize",
                    "WholeCerebellum_ClustSize",
                    "CompositeWM_ClustSize",
                    "Left_Cerebral_White_Matter_ClustSize",
                    "Left_Lateral_Ventricle_ClustSize",
                    "Left_Inf_Lat_Vent_ClustSize",
                    "Left_Cerebellum_White_Matter_ClustSize",
                    "Left_Cerebellum_Cortex_ClustSize",
                    "Left_Thalamus_Proper_ClustSize",
                    "Left_Caudate_ClustSize",
                    "Left_Putamen_ClustSize",
                    "Left_Pallidum_ClustSize",
                    "Third_Ventricle_ClustSize",
                    "Fourth_Ventricle_ClustSize",
                    "Brain_Stem_ClustSize",
                    "Left_Hippocampus_ClustSize",
                    "Left_Amygdala_ClustSize",
                    "CSF_ClustSize",
                    "Left_Accumbens_area_ClustSize",
                    "Left_VentralDC_ClustSize",
                    "Left_vessel_ClustSize",
                    "Left_choroid_plexus_ClustSize",
                    "Right_Cerebral_White_Matter_ClustSize",
                    "Right_Lateral_Ventricle_ClustSize",
                    "Right_Inf_Lat_Vent_ClustSize",
                    "Right_Cerebellum_White_Matter_ClustSize",
                    "Right_Cerebellum_Cortex_ClustSize",
                    "Right_Thalamus_Proper_ClustSize",
                    "Right_Caudate_ClustSize",
                    "Right_Putamen_ClustSize",
                    "Right_Pallidum_ClustSize",
                    "Right_Hippocampus_ClustSize",
                    "Right_Amygdala_ClustSize",
                    "Right_Accumbens_area_ClustSize",
                    "Right_VentralDC_ClustSize",
                    "Right_vessel_ClustSize",
                    "Right_choroid_plexus_ClustSize",
                    "Fifth_Ventricle_ClustSize",
                    "WM_hypointensities_ClustSize",
                    "non_WM_hypointensities_ClustSize",
                    "Optic_Chiasm_ClustSize",
                    "CC_Posterior_ClustSize",
                    "CC_Mid_Posterior_ClustSize",
                    "CC_Central_ClustSize",
                    "CC_Mid_Anterior_ClustSize",
                    "CC_Anterior_ClustSize",
                    "ctx_lh_unknown_ClustSize",
                    "ctx_lh_bankssts_ClustSize",
                    "ctx_lh_caudalanteriorcingulate_ClustSize",
                    "ctx_lh_caudalmiddlefrontal_ClustSize",
                    "ctx_lh_cuneus_ClustSize",
                    "ctx_lh_entorhinal_ClustSize",
                    "ctx_lh_fusiform_ClustSize",
                    "ctx_lh_inferiorparietal_ClustSize",
                    "ctx_lh_inferiortemporal_ClustSize",
                    "ctx_lh_isthmuscingulate_ClustSize",
                    "ctx_lh_lateraloccipital_ClustSize",
                    "ctx_lh_lateralorbitofrontal_ClustSize",
                    "ctx_lh_lingual_ClustSize",
                    "ctx_lh_medialorbitofrontal_ClustSize",
                    "ctx_lh_middletemporal_ClustSize",
                    "ctx_lh_parahippocampal_ClustSize",
                    "ctx_lh_paracentral_ClustSize",
                    "ctx_lh_parsopercularis_ClustSize",
                    "ctx_lh_parsorbitalis_ClustSize",
                    "ctx_lh_parstriangularis_ClustSize",
                    "ctx_lh_pericalcarine_ClustSize",
                    "ctx_lh_postcentral_ClustSize",
                    "ctx_lh_posteriorcingulate_ClustSize",
                    "ctx_lh_precentral_ClustSize",
                    "ctx_lh_precuneus_ClustSize",
                    "ctx_lh_rostralanteriorcingulate_ClustSize",
                    "ctx_lh_rostralmiddlefrontal_ClustSize",
                    "ctx_lh_superiorfrontal_ClustSize",
                    "ctx_lh_superiorparietal_ClustSize",
                    "ctx_lh_superiortemporal_ClustSize",
                    "ctx_lh_supramarginal_ClustSize",
                    "ctx_lh_frontalpole_ClustSize",
                    "ctx_lh_temporalpole_ClustSize",
                    "ctx_lh_transversetemporal_ClustSize",
                    "ctx_lh_insula_ClustSize",
                    "ctx_rh_unknown_ClustSize",
                    "ctx_rh_bankssts_ClustSize",
                    "ctx_rh_caudalanteriorcingulate_ClustSize",
                    "ctx_rh_caudalmiddlefrontal_ClustSize",
                    "ctx_rh_cuneus_ClustSize",
                    "ctx_rh_entorhinal_ClustSize",
                    "ctx_rh_fusiform_ClustSize",
                    "ctx_rh_inferiorparietal_ClustSize",
                    "ctx_rh_inferiortemporal_ClustSize",
                    "ctx_rh_isthmuscingulate_ClustSize",
                    "ctx_rh_lateraloccipital_ClustSize",
                    "ctx_rh_lateralorbitofrontal_ClustSize",
                    "ctx_rh_lingual_ClustSize",
                    "ctx_rh_medialorbitofrontal_ClustSize",
                    "ctx_rh_middletemporal_ClustSize",
                    "ctx_rh_parahippocampal_ClustSize",
                    "ctx_rh_paracentral_ClustSize",
                    "ctx_rh_parsopercularis_ClustSize",
                    "ctx_rh_parsorbitalis_ClustSize",
                    "ctx_rh_parstriangularis_ClustSize",
                    "ctx_rh_pericalcarine_ClustSize",
                    "ctx_rh_postcentral_ClustSize",
                    "ctx_rh_posteriorcingulate_ClustSize",
                    "ctx_rh_precentral_ClustSize",
                    "ctx_rh_precuneus_ClustSize",
                    "ctx_rh_rostralanteriorcingulate_ClustSize",
                    "ctx_rh_rostralmiddlefrontal_ClustSize",
                    "ctx_rh_superiorfrontal_ClustSize",
                    "ctx_rh_superiorparietal_ClustSize",
                    "ctx_rh_superiortemporal_ClustSize",
                    "ctx_rh_supramarginal_ClustSize",
                    "ctx_rh_frontalpole_ClustSize",
                    "ctx_rh_temporalpole_ClustSize",
                    "ctx_rh_transversetemporal_ClustSize",
                    "ctx_rh_insula_ClustSize",
                ]
                self.qrep_dat[tracer] = self.qrep_dat[tracer][cols_in_order]

                # Save the output dataframe
                output_file = op.join(
                    self.paths["qreport"],
                    f"LEADS-PETCore-quarterly-report_{self.report_period}_{tracer.upper()}-ROI-means_{TODAY}.csv",
                )
                if save_output:
                    if overwrite or not op.isfile(output_file):
                        self.qrep_dat[tracer].to_csv(output_file, index=False)
                        print(f"Saved {output_file}")

            def format_fdg(self):
                # Drop unnecessary columns
                drop_cols = [
                    "tracer",
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
                    "pet_qc_pass",
                    "native_pet_ok",
                    "pet_to_mri_coreg_ok",
                    "pons_mask_ok",
                    "warped_pet_ok",
                    "pet_qc_notes",
                    "mri_qc_pass",
                    "native_nu_rating",
                    "aparc_rating",
                    "warped_nu_ok",
                    "mri_qc_notes",
                    "pons_MRIBASED_SUVR",
                    "meta_temporal_MRIBASED_SUVR",
                    "mtl_no_hippocampus_MRIBASED_SUVR",
                    "basolateral_temporal_MRIBASED_SUVR",
                    "temporoparietal_MRIBASED_SUVR",
                    "ctx_desikan_MRIBASED_SUVR",
                    "pons_ClustSize",
                    "meta_temporal_ClustSize",
                    "mtl_no_hippocampus_ClustSize",
                    "basolateral_temporal_ClustSize",
                    "temporoparietal_ClustSize",
                    "ctx_desikan_ClustSize",
                ]
                self.qrep_dat[tracer] = self.qrep_dat[tracer].drop(columns=drop_cols)

                # Rename columns
                self.qrep_dat[tracer] = self.qrep_dat[tracer].rename(
                    columns={
                        "subject_id": "ID",
                        "pet_date": "FDGPET_Date",
                        "pet_image_id": "ImageID",
                        "ScalingFactor_pons": "ScalingFactor_Pons",
                    }
                )

                # Add empty columns
                add_cols = [
                    "Left_vessel_MRIBASED_SUVR",
                    "Right_vessel_MRIBASED_SUVR",
                    "Fifth_Ventricle_MRIBASED_SUVR",
                    "non_WM_hypointensities_MRIBASED_SUVR",
                    "ctx_lh_unknown_MRIBASED_SUVR",
                    "ctx_rh_unknown_MRIBASED_SUVR",
                    "ScalingFactor_Pons_ClustSize",
                    "Left_vessel_ClustSize",
                    "Right_vessel_ClustSize",
                    "Fifth_Ventricle_ClustSize",
                    "non_WM_hypointensities_ClustSize",
                    "ctx_lh_unknown_ClustSize",
                    "ctx_rh_unknown_ClustSize",
                ]
                for col in add_cols:
                    self.qrep_dat[tracer][col] = np.nan

                # Round scaling factor columns to 4 decimals
                round_cols = [
                    col
                    for col in self.qrep_dat[tracer].columns
                    if col.startswith("ScalingFactor")
                ]
                self.qrep_dat[tracer][round_cols] = self.qrep_dat[tracer][
                    round_cols
                ].round(4)

                # Round MRIBASED_SUVR columns to 7 decimals
                round_cols = [
                    col
                    for col in self.qrep_dat[tracer].columns
                    if col.endswith("MRIBASED_SUVR")
                ]
                self.qrep_dat[tracer][round_cols] = self.qrep_dat[tracer][
                    round_cols
                ].round(7)

                # Reorder columns
                cols_in_order = [
                    "ID",
                    "FDGPET_Date",
                    "ImageID",
                    "CohortAssgn",
                    "ScalingFactor_Pons",
                    "Left_Cerebral_White_Matter_MRIBASED_SUVR",
                    "Left_Lateral_Ventricle_MRIBASED_SUVR",
                    "Left_Inf_Lat_Vent_MRIBASED_SUVR",
                    "Left_Cerebellum_White_Matter_MRIBASED_SUVR",
                    "Left_Cerebellum_Cortex_MRIBASED_SUVR",
                    "Left_Thalamus_Proper_MRIBASED_SUVR",
                    "Left_Caudate_MRIBASED_SUVR",
                    "Left_Putamen_MRIBASED_SUVR",
                    "Left_Pallidum_MRIBASED_SUVR",
                    "Third_Ventricle_MRIBASED_SUVR",
                    "Fourth_Ventricle_MRIBASED_SUVR",
                    "Brain_Stem_MRIBASED_SUVR",
                    "Left_Hippocampus_MRIBASED_SUVR",
                    "Left_Amygdala_MRIBASED_SUVR",
                    "CSF_MRIBASED_SUVR",
                    "Left_Accumbens_area_MRIBASED_SUVR",
                    "Left_VentralDC_MRIBASED_SUVR",
                    "Left_vessel_MRIBASED_SUVR",
                    "Left_choroid_plexus_MRIBASED_SUVR",
                    "Right_Cerebral_White_Matter_MRIBASED_SUVR",
                    "Right_Lateral_Ventricle_MRIBASED_SUVR",
                    "Right_Inf_Lat_Vent_MRIBASED_SUVR",
                    "Right_Cerebellum_White_Matter_MRIBASED_SUVR",
                    "Right_Cerebellum_Cortex_MRIBASED_SUVR",
                    "Right_Thalamus_Proper_MRIBASED_SUVR",
                    "Right_Caudate_MRIBASED_SUVR",
                    "Right_Putamen_MRIBASED_SUVR",
                    "Right_Pallidum_MRIBASED_SUVR",
                    "Right_Hippocampus_MRIBASED_SUVR",
                    "Right_Amygdala_MRIBASED_SUVR",
                    "Right_Accumbens_area_MRIBASED_SUVR",
                    "Right_VentralDC_MRIBASED_SUVR",
                    "Right_vessel_MRIBASED_SUVR",
                    "Right_choroid_plexus_MRIBASED_SUVR",
                    "Fifth_Ventricle_MRIBASED_SUVR",
                    "WM_hypointensities_MRIBASED_SUVR",
                    "non_WM_hypointensities_MRIBASED_SUVR",
                    "Optic_Chiasm_MRIBASED_SUVR",
                    "CC_Posterior_MRIBASED_SUVR",
                    "CC_Mid_Posterior_MRIBASED_SUVR",
                    "CC_Central_MRIBASED_SUVR",
                    "CC_Mid_Anterior_MRIBASED_SUVR",
                    "CC_Anterior_MRIBASED_SUVR",
                    "ctx_lh_unknown_MRIBASED_SUVR",
                    "ctx_lh_bankssts_MRIBASED_SUVR",
                    "ctx_lh_caudalanteriorcingulate_MRIBASED_SUVR",
                    "ctx_lh_caudalmiddlefrontal_MRIBASED_SUVR",
                    "ctx_lh_cuneus_MRIBASED_SUVR",
                    "ctx_lh_entorhinal_MRIBASED_SUVR",
                    "ctx_lh_fusiform_MRIBASED_SUVR",
                    "ctx_lh_inferiorparietal_MRIBASED_SUVR",
                    "ctx_lh_inferiortemporal_MRIBASED_SUVR",
                    "ctx_lh_isthmuscingulate_MRIBASED_SUVR",
                    "ctx_lh_lateraloccipital_MRIBASED_SUVR",
                    "ctx_lh_lateralorbitofrontal_MRIBASED_SUVR",
                    "ctx_lh_lingual_MRIBASED_SUVR",
                    "ctx_lh_medialorbitofrontal_MRIBASED_SUVR",
                    "ctx_lh_middletemporal_MRIBASED_SUVR",
                    "ctx_lh_parahippocampal_MRIBASED_SUVR",
                    "ctx_lh_paracentral_MRIBASED_SUVR",
                    "ctx_lh_parsopercularis_MRIBASED_SUVR",
                    "ctx_lh_parsorbitalis_MRIBASED_SUVR",
                    "ctx_lh_parstriangularis_MRIBASED_SUVR",
                    "ctx_lh_pericalcarine_MRIBASED_SUVR",
                    "ctx_lh_postcentral_MRIBASED_SUVR",
                    "ctx_lh_posteriorcingulate_MRIBASED_SUVR",
                    "ctx_lh_precentral_MRIBASED_SUVR",
                    "ctx_lh_precuneus_MRIBASED_SUVR",
                    "ctx_lh_rostralanteriorcingulate_MRIBASED_SUVR",
                    "ctx_lh_rostralmiddlefrontal_MRIBASED_SUVR",
                    "ctx_lh_superiorfrontal_MRIBASED_SUVR",
                    "ctx_lh_superiorparietal_MRIBASED_SUVR",
                    "ctx_lh_superiortemporal_MRIBASED_SUVR",
                    "ctx_lh_supramarginal_MRIBASED_SUVR",
                    "ctx_lh_frontalpole_MRIBASED_SUVR",
                    "ctx_lh_temporalpole_MRIBASED_SUVR",
                    "ctx_lh_transversetemporal_MRIBASED_SUVR",
                    "ctx_lh_insula_MRIBASED_SUVR",
                    "ctx_rh_unknown_MRIBASED_SUVR",
                    "ctx_rh_bankssts_MRIBASED_SUVR",
                    "ctx_rh_caudalanteriorcingulate_MRIBASED_SUVR",
                    "ctx_rh_caudalmiddlefrontal_MRIBASED_SUVR",
                    "ctx_rh_cuneus_MRIBASED_SUVR",
                    "ctx_rh_entorhinal_MRIBASED_SUVR",
                    "ctx_rh_fusiform_MRIBASED_SUVR",
                    "ctx_rh_inferiorparietal_MRIBASED_SUVR",
                    "ctx_rh_inferiortemporal_MRIBASED_SUVR",
                    "ctx_rh_isthmuscingulate_MRIBASED_SUVR",
                    "ctx_rh_lateraloccipital_MRIBASED_SUVR",
                    "ctx_rh_lateralorbitofrontal_MRIBASED_SUVR",
                    "ctx_rh_lingual_MRIBASED_SUVR",
                    "ctx_rh_medialorbitofrontal_MRIBASED_SUVR",
                    "ctx_rh_middletemporal_MRIBASED_SUVR",
                    "ctx_rh_parahippocampal_MRIBASED_SUVR",
                    "ctx_rh_paracentral_MRIBASED_SUVR",
                    "ctx_rh_parsopercularis_MRIBASED_SUVR",
                    "ctx_rh_parsorbitalis_MRIBASED_SUVR",
                    "ctx_rh_parstriangularis_MRIBASED_SUVR",
                    "ctx_rh_pericalcarine_MRIBASED_SUVR",
                    "ctx_rh_postcentral_MRIBASED_SUVR",
                    "ctx_rh_posteriorcingulate_MRIBASED_SUVR",
                    "ctx_rh_precentral_MRIBASED_SUVR",
                    "ctx_rh_precuneus_MRIBASED_SUVR",
                    "ctx_rh_rostralanteriorcingulate_MRIBASED_SUVR",
                    "ctx_rh_rostralmiddlefrontal_MRIBASED_SUVR",
                    "ctx_rh_superiorfrontal_MRIBASED_SUVR",
                    "ctx_rh_superiorparietal_MRIBASED_SUVR",
                    "ctx_rh_superiortemporal_MRIBASED_SUVR",
                    "ctx_rh_supramarginal_MRIBASED_SUVR",
                    "ctx_rh_frontalpole_MRIBASED_SUVR",
                    "ctx_rh_temporalpole_MRIBASED_SUVR",
                    "ctx_rh_transversetemporal_MRIBASED_SUVR",
                    "ctx_rh_insula_MRIBASED_SUVR",
                    "ScalingFactor_Pons_ClustSize",
                    "Left_Cerebral_White_Matter_ClustSize",
                    "Left_Lateral_Ventricle_ClustSize",
                    "Left_Inf_Lat_Vent_ClustSize",
                    "Left_Cerebellum_White_Matter_ClustSize",
                    "Left_Cerebellum_Cortex_ClustSize",
                    "Left_Thalamus_Proper_ClustSize",
                    "Left_Caudate_ClustSize",
                    "Left_Putamen_ClustSize",
                    "Left_Pallidum_ClustSize",
                    "Third_Ventricle_ClustSize",
                    "Fourth_Ventricle_ClustSize",
                    "Brain_Stem_ClustSize",
                    "Left_Hippocampus_ClustSize",
                    "Left_Amygdala_ClustSize",
                    "CSF_ClustSize",
                    "Left_Accumbens_area_ClustSize",
                    "Left_VentralDC_ClustSize",
                    "Left_vessel_ClustSize",
                    "Left_choroid_plexus_ClustSize",
                    "Right_Cerebral_White_Matter_ClustSize",
                    "Right_Lateral_Ventricle_ClustSize",
                    "Right_Inf_Lat_Vent_ClustSize",
                    "Right_Cerebellum_White_Matter_ClustSize",
                    "Right_Cerebellum_Cortex_ClustSize",
                    "Right_Thalamus_Proper_ClustSize",
                    "Right_Caudate_ClustSize",
                    "Right_Putamen_ClustSize",
                    "Right_Pallidum_ClustSize",
                    "Right_Hippocampus_ClustSize",
                    "Right_Amygdala_ClustSize",
                    "Right_Accumbens_area_ClustSize",
                    "Right_VentralDC_ClustSize",
                    "Right_vessel_ClustSize",
                    "Right_choroid_plexus_ClustSize",
                    "Fifth_Ventricle_ClustSize",
                    "WM_hypointensities_ClustSize",
                    "non_WM_hypointensities_ClustSize",
                    "Optic_Chiasm_ClustSize",
                    "CC_Posterior_ClustSize",
                    "CC_Mid_Posterior_ClustSize",
                    "CC_Central_ClustSize",
                    "CC_Mid_Anterior_ClustSize",
                    "CC_Anterior_ClustSize",
                    "ctx_lh_unknown_ClustSize",
                    "ctx_lh_bankssts_ClustSize",
                    "ctx_lh_caudalanteriorcingulate_ClustSize",
                    "ctx_lh_caudalmiddlefrontal_ClustSize",
                    "ctx_lh_cuneus_ClustSize",
                    "ctx_lh_entorhinal_ClustSize",
                    "ctx_lh_fusiform_ClustSize",
                    "ctx_lh_inferiorparietal_ClustSize",
                    "ctx_lh_inferiortemporal_ClustSize",
                    "ctx_lh_isthmuscingulate_ClustSize",
                    "ctx_lh_lateraloccipital_ClustSize",
                    "ctx_lh_lateralorbitofrontal_ClustSize",
                    "ctx_lh_lingual_ClustSize",
                    "ctx_lh_medialorbitofrontal_ClustSize",
                    "ctx_lh_middletemporal_ClustSize",
                    "ctx_lh_parahippocampal_ClustSize",
                    "ctx_lh_paracentral_ClustSize",
                    "ctx_lh_parsopercularis_ClustSize",
                    "ctx_lh_parsorbitalis_ClustSize",
                    "ctx_lh_parstriangularis_ClustSize",
                    "ctx_lh_pericalcarine_ClustSize",
                    "ctx_lh_postcentral_ClustSize",
                    "ctx_lh_posteriorcingulate_ClustSize",
                    "ctx_lh_precentral_ClustSize",
                    "ctx_lh_precuneus_ClustSize",
                    "ctx_lh_rostralanteriorcingulate_ClustSize",
                    "ctx_lh_rostralmiddlefrontal_ClustSize",
                    "ctx_lh_superiorfrontal_ClustSize",
                    "ctx_lh_superiorparietal_ClustSize",
                    "ctx_lh_superiortemporal_ClustSize",
                    "ctx_lh_supramarginal_ClustSize",
                    "ctx_lh_frontalpole_ClustSize",
                    "ctx_lh_temporalpole_ClustSize",
                    "ctx_lh_transversetemporal_ClustSize",
                    "ctx_lh_insula_ClustSize",
                    "ctx_rh_unknown_ClustSize",
                    "ctx_rh_bankssts_ClustSize",
                    "ctx_rh_caudalanteriorcingulate_ClustSize",
                    "ctx_rh_caudalmiddlefrontal_ClustSize",
                    "ctx_rh_cuneus_ClustSize",
                    "ctx_rh_entorhinal_ClustSize",
                    "ctx_rh_fusiform_ClustSize",
                    "ctx_rh_inferiorparietal_ClustSize",
                    "ctx_rh_inferiortemporal_ClustSize",
                    "ctx_rh_isthmuscingulate_ClustSize",
                    "ctx_rh_lateraloccipital_ClustSize",
                    "ctx_rh_lateralorbitofrontal_ClustSize",
                    "ctx_rh_lingual_ClustSize",
                    "ctx_rh_medialorbitofrontal_ClustSize",
                    "ctx_rh_middletemporal_ClustSize",
                    "ctx_rh_parahippocampal_ClustSize",
                    "ctx_rh_paracentral_ClustSize",
                    "ctx_rh_parsopercularis_ClustSize",
                    "ctx_rh_parsorbitalis_ClustSize",
                    "ctx_rh_parstriangularis_ClustSize",
                    "ctx_rh_pericalcarine_ClustSize",
                    "ctx_rh_postcentral_ClustSize",
                    "ctx_rh_posteriorcingulate_ClustSize",
                    "ctx_rh_precentral_ClustSize",
                    "ctx_rh_precuneus_ClustSize",
                    "ctx_rh_rostralanteriorcingulate_ClustSize",
                    "ctx_rh_rostralmiddlefrontal_ClustSize",
                    "ctx_rh_superiorfrontal_ClustSize",
                    "ctx_rh_superiorparietal_ClustSize",
                    "ctx_rh_superiortemporal_ClustSize",
                    "ctx_rh_supramarginal_ClustSize",
                    "ctx_rh_frontalpole_ClustSize",
                    "ctx_rh_temporalpole_ClustSize",
                    "ctx_rh_transversetemporal_ClustSize",
                    "ctx_rh_insula_ClustSize",
                ]
                self.qrep_dat[tracer] = self.qrep_dat[tracer][cols_in_order]

                # Save the output dataframe
                output_file = op.join(
                    self.paths["qreport"],
                    f"LEADS-PETCore-quarterly-report_{self.report_period}_{tracer.upper()}-ROI-means_{TODAY}.csv",
                )
                if save_output:
                    if overwrite or not op.isfile(output_file):
                        self.qrep_dat[tracer].to_csv(output_file, index=False)
                        print(f"Saved {output_file}")

            def format_ftp(self):
                # Drop unnecessary columns
                drop_cols = [
                    "tracer",
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
                    "pet_qc_pass",
                    "native_pet_ok",
                    "pet_to_mri_coreg_ok",
                    "infcblgm_mask_ok",
                    "erodedwm_mask_ok",
                    "warped_pet_ok",
                    "pet_qc_notes",
                    "mri_qc_pass",
                    "native_nu_rating",
                    "aparc_rating",
                    "warped_nu_ok",
                    "mri_qc_notes",
                    "infcblgm_MRIBASED_SUVR",
                    "mtl_no_hippocampus_MRIBASED_SUVR",
                    "basolateral_temporal_MRIBASED_SUVR",
                    "temporoparietal_MRIBASED_SUVR",
                    "ctx_desikan_MRIBASED_SUVR",
                    "mtl_no_hippocampus_ClustSize",
                    "basolateral_temporal_ClustSize",
                    "temporoparietal_ClustSize",
                    "ctx_desikan_ClustSize",
                ]
                self.qrep_dat[tracer] = self.qrep_dat[tracer].drop(columns=drop_cols)

                # Rename columns
                self.qrep_dat[tracer] = self.qrep_dat[tracer].rename(
                    columns={
                        "subject_id": "ID",
                        "pet_date": "FTPPET_Date",
                        "pet_image_id": "ImageID",
                        "ScalingFactor_infcblgm": "ScalingFactor_InfCerebGray",
                        "ScalingFactor_eroded": "ScalingFactor_ErodedWM",
                        "eroded_subcortwm_MRIBASED_SUVR": "ErodedWM_MRIBASED_SUVR",
                        "meta_temporal_MRIBASED_SUVR": "MetaROI_MRIBASED_SUVR",
                        "eroded_subcortwm_ClustSize": "ScalingFactor_ErodedWM_ClustSize",
                        "infcblgm_ClustSize": "ScalingFactor_InfCerebGray_ClustSize",
                        "meta_temporal_ClustSize": "MetaROI_ClustSize",
                    }
                )

                # Add empty columns
                add_cols = [
                    "Assigned_MRIBASED_MetaROI_ADNIcutoff_1p2",
                    "Braak_1_MRIBASED_SUVR",
                    "Braak_12_MRIBASED_SUVR",
                    "Braak_34_MRIBASED_SUVR",
                    "Braak_56_MRIBASED_SUVR",
                    "Left_vessel_MRIBASED_SUVR",
                    "Right_vessel_MRIBASED_SUVR",
                    "Fifth_Ventricle_MRIBASED_SUVR",
                    "non_WM_hypointensities_MRIBASED_SUVR",
                    "ctx_lh_unknown_MRIBASED_SUVR",
                    "ctx_rh_unknown_MRIBASED_SUVR",
                    "Braak_1_ClustSize",
                    "Braak_12_ClustSize",
                    "Braak_34_ClustSize",
                    "Braak_56_ClustSize",
                    "ErodedWM_ClustSize",
                    "Left_vessel_ClustSize",
                    "Right_vessel_ClustSize",
                    "Fifth_Ventricle_ClustSize",
                    "non_WM_hypointensities_ClustSize",
                    "ctx_lh_unknown_ClustSize",
                    "ctx_rh_unknown_ClustSize",
                ]
                for col in add_cols:
                    self.qrep_dat[tracer][col] = np.nan

                # Convert the scan date column to a string formatted like '%m/%d/%y'
                self.qrep_dat[tracer]["FTPPET_Date"] = self.qrep_dat[tracer][
                    "FTPPET_Date"
                ].apply(lambda x: pd.to_datetime(x).strftime("%m/%d/%y"))

                # Round scaling factor columns to 4 decimals
                round_cols = [
                    col
                    for col in self.qrep_dat[tracer].columns
                    if col.startswith("ScalingFactor")
                ]
                self.qrep_dat[tracer][round_cols] = self.qrep_dat[tracer][
                    round_cols
                ].round(4)

                # Round MRIBASED_SUVR columns to 7 decimals
                round_cols = [
                    col
                    for col in self.qrep_dat[tracer].columns
                    if col.endswith("MRIBASED_SUVR")
                ]
                self.qrep_dat[tracer][round_cols] = self.qrep_dat[tracer][
                    round_cols
                ].round(7)

                # Reorder columns
                cols_in_order = [
                    "ID",
                    "FTPPET_Date",
                    "ImageID",
                    "CohortAssgn",
                    "ScalingFactor_InfCerebGray",
                    "ScalingFactor_ErodedWM",
                    "Assigned_MRIBASED_MetaROI_ADNIcutoff_1p2",
                    "MetaROI_MRIBASED_SUVR",
                    "Braak_1_MRIBASED_SUVR",
                    "Braak_12_MRIBASED_SUVR",
                    "Braak_34_MRIBASED_SUVR",
                    "Braak_56_MRIBASED_SUVR",
                    "ErodedWM_MRIBASED_SUVR",
                    "Left_Cerebral_White_Matter_MRIBASED_SUVR",
                    "Left_Lateral_Ventricle_MRIBASED_SUVR",
                    "Left_Inf_Lat_Vent_MRIBASED_SUVR",
                    "Left_Cerebellum_White_Matter_MRIBASED_SUVR",
                    "Left_Cerebellum_Cortex_MRIBASED_SUVR",
                    "Left_Thalamus_Proper_MRIBASED_SUVR",
                    "Left_Caudate_MRIBASED_SUVR",
                    "Left_Putamen_MRIBASED_SUVR",
                    "Left_Pallidum_MRIBASED_SUVR",
                    "Third_Ventricle_MRIBASED_SUVR",
                    "Fourth_Ventricle_MRIBASED_SUVR",
                    "Brain_Stem_MRIBASED_SUVR",
                    "Left_Hippocampus_MRIBASED_SUVR",
                    "Left_Amygdala_MRIBASED_SUVR",
                    "CSF_MRIBASED_SUVR",
                    "Left_Accumbens_area_MRIBASED_SUVR",
                    "Left_VentralDC_MRIBASED_SUVR",
                    "Left_vessel_MRIBASED_SUVR",
                    "Left_choroid_plexus_MRIBASED_SUVR",
                    "Right_Cerebral_White_Matter_MRIBASED_SUVR",
                    "Right_Lateral_Ventricle_MRIBASED_SUVR",
                    "Right_Inf_Lat_Vent_MRIBASED_SUVR",
                    "Right_Cerebellum_White_Matter_MRIBASED_SUVR",
                    "Right_Cerebellum_Cortex_MRIBASED_SUVR",
                    "Right_Thalamus_Proper_MRIBASED_SUVR",
                    "Right_Caudate_MRIBASED_SUVR",
                    "Right_Putamen_MRIBASED_SUVR",
                    "Right_Pallidum_MRIBASED_SUVR",
                    "Right_Hippocampus_MRIBASED_SUVR",
                    "Right_Amygdala_MRIBASED_SUVR",
                    "Right_Accumbens_area_MRIBASED_SUVR",
                    "Right_VentralDC_MRIBASED_SUVR",
                    "Right_vessel_MRIBASED_SUVR",
                    "Right_choroid_plexus_MRIBASED_SUVR",
                    "Fifth_Ventricle_MRIBASED_SUVR",
                    "WM_hypointensities_MRIBASED_SUVR",
                    "non_WM_hypointensities_MRIBASED_SUVR",
                    "Optic_Chiasm_MRIBASED_SUVR",
                    "CC_Posterior_MRIBASED_SUVR",
                    "CC_Mid_Posterior_MRIBASED_SUVR",
                    "CC_Central_MRIBASED_SUVR",
                    "CC_Mid_Anterior_MRIBASED_SUVR",
                    "CC_Anterior_MRIBASED_SUVR",
                    "ctx_lh_unknown_MRIBASED_SUVR",
                    "ctx_lh_bankssts_MRIBASED_SUVR",
                    "ctx_lh_caudalanteriorcingulate_MRIBASED_SUVR",
                    "ctx_lh_caudalmiddlefrontal_MRIBASED_SUVR",
                    "ctx_lh_cuneus_MRIBASED_SUVR",
                    "ctx_lh_entorhinal_MRIBASED_SUVR",
                    "ctx_lh_fusiform_MRIBASED_SUVR",
                    "ctx_lh_inferiorparietal_MRIBASED_SUVR",
                    "ctx_lh_inferiortemporal_MRIBASED_SUVR",
                    "ctx_lh_isthmuscingulate_MRIBASED_SUVR",
                    "ctx_lh_lateraloccipital_MRIBASED_SUVR",
                    "ctx_lh_lateralorbitofrontal_MRIBASED_SUVR",
                    "ctx_lh_lingual_MRIBASED_SUVR",
                    "ctx_lh_medialorbitofrontal_MRIBASED_SUVR",
                    "ctx_lh_middletemporal_MRIBASED_SUVR",
                    "ctx_lh_parahippocampal_MRIBASED_SUVR",
                    "ctx_lh_paracentral_MRIBASED_SUVR",
                    "ctx_lh_parsopercularis_MRIBASED_SUVR",
                    "ctx_lh_parsorbitalis_MRIBASED_SUVR",
                    "ctx_lh_parstriangularis_MRIBASED_SUVR",
                    "ctx_lh_pericalcarine_MRIBASED_SUVR",
                    "ctx_lh_postcentral_MRIBASED_SUVR",
                    "ctx_lh_posteriorcingulate_MRIBASED_SUVR",
                    "ctx_lh_precentral_MRIBASED_SUVR",
                    "ctx_lh_precuneus_MRIBASED_SUVR",
                    "ctx_lh_rostralanteriorcingulate_MRIBASED_SUVR",
                    "ctx_lh_rostralmiddlefrontal_MRIBASED_SUVR",
                    "ctx_lh_superiorfrontal_MRIBASED_SUVR",
                    "ctx_lh_superiorparietal_MRIBASED_SUVR",
                    "ctx_lh_superiortemporal_MRIBASED_SUVR",
                    "ctx_lh_supramarginal_MRIBASED_SUVR",
                    "ctx_lh_frontalpole_MRIBASED_SUVR",
                    "ctx_lh_temporalpole_MRIBASED_SUVR",
                    "ctx_lh_transversetemporal_MRIBASED_SUVR",
                    "ctx_lh_insula_MRIBASED_SUVR",
                    "ctx_rh_unknown_MRIBASED_SUVR",
                    "ctx_rh_bankssts_MRIBASED_SUVR",
                    "ctx_rh_caudalanteriorcingulate_MRIBASED_SUVR",
                    "ctx_rh_caudalmiddlefrontal_MRIBASED_SUVR",
                    "ctx_rh_cuneus_MRIBASED_SUVR",
                    "ctx_rh_entorhinal_MRIBASED_SUVR",
                    "ctx_rh_fusiform_MRIBASED_SUVR",
                    "ctx_rh_inferiorparietal_MRIBASED_SUVR",
                    "ctx_rh_inferiortemporal_MRIBASED_SUVR",
                    "ctx_rh_isthmuscingulate_MRIBASED_SUVR",
                    "ctx_rh_lateraloccipital_MRIBASED_SUVR",
                    "ctx_rh_lateralorbitofrontal_MRIBASED_SUVR",
                    "ctx_rh_lingual_MRIBASED_SUVR",
                    "ctx_rh_medialorbitofrontal_MRIBASED_SUVR",
                    "ctx_rh_middletemporal_MRIBASED_SUVR",
                    "ctx_rh_parahippocampal_MRIBASED_SUVR",
                    "ctx_rh_paracentral_MRIBASED_SUVR",
                    "ctx_rh_parsopercularis_MRIBASED_SUVR",
                    "ctx_rh_parsorbitalis_MRIBASED_SUVR",
                    "ctx_rh_parstriangularis_MRIBASED_SUVR",
                    "ctx_rh_pericalcarine_MRIBASED_SUVR",
                    "ctx_rh_postcentral_MRIBASED_SUVR",
                    "ctx_rh_posteriorcingulate_MRIBASED_SUVR",
                    "ctx_rh_precentral_MRIBASED_SUVR",
                    "ctx_rh_precuneus_MRIBASED_SUVR",
                    "ctx_rh_rostralanteriorcingulate_MRIBASED_SUVR",
                    "ctx_rh_rostralmiddlefrontal_MRIBASED_SUVR",
                    "ctx_rh_superiorfrontal_MRIBASED_SUVR",
                    "ctx_rh_superiorparietal_MRIBASED_SUVR",
                    "ctx_rh_superiortemporal_MRIBASED_SUVR",
                    "ctx_rh_supramarginal_MRIBASED_SUVR",
                    "ctx_rh_frontalpole_MRIBASED_SUVR",
                    "ctx_rh_temporalpole_MRIBASED_SUVR",
                    "ctx_rh_transversetemporal_MRIBASED_SUVR",
                    "ctx_rh_insula_MRIBASED_SUVR",
                    "ScalingFactor_InfCerebGray_ClustSize",
                    "ScalingFactor_ErodedWM_ClustSize",
                    "MetaROI_ClustSize",
                    "Braak_1_ClustSize",
                    "Braak_12_ClustSize",
                    "Braak_34_ClustSize",
                    "Braak_56_ClustSize",
                    "ErodedWM_ClustSize",
                    "Left_Cerebral_White_Matter_ClustSize",
                    "Left_Lateral_Ventricle_ClustSize",
                    "Left_Inf_Lat_Vent_ClustSize",
                    "Left_Cerebellum_White_Matter_ClustSize",
                    "Left_Cerebellum_Cortex_ClustSize",
                    "Left_Thalamus_Proper_ClustSize",
                    "Left_Caudate_ClustSize",
                    "Left_Putamen_ClustSize",
                    "Left_Pallidum_ClustSize",
                    "Third_Ventricle_ClustSize",
                    "Fourth_Ventricle_ClustSize",
                    "Brain_Stem_ClustSize",
                    "Left_Hippocampus_ClustSize",
                    "Left_Amygdala_ClustSize",
                    "CSF_ClustSize",
                    "Left_Accumbens_area_ClustSize",
                    "Left_VentralDC_ClustSize",
                    "Left_vessel_ClustSize",
                    "Left_choroid_plexus_ClustSize",
                    "Right_Cerebral_White_Matter_ClustSize",
                    "Right_Lateral_Ventricle_ClustSize",
                    "Right_Inf_Lat_Vent_ClustSize",
                    "Right_Cerebellum_White_Matter_ClustSize",
                    "Right_Cerebellum_Cortex_ClustSize",
                    "Right_Thalamus_Proper_ClustSize",
                    "Right_Caudate_ClustSize",
                    "Right_Putamen_ClustSize",
                    "Right_Pallidum_ClustSize",
                    "Right_Hippocampus_ClustSize",
                    "Right_Amygdala_ClustSize",
                    "Right_Accumbens_area_ClustSize",
                    "Right_VentralDC_ClustSize",
                    "Right_vessel_ClustSize",
                    "Right_choroid_plexus_ClustSize",
                    "Fifth_Ventricle_ClustSize",
                    "WM_hypointensities_ClustSize",
                    "non_WM_hypointensities_ClustSize",
                    "Optic_Chiasm_ClustSize",
                    "CC_Posterior_ClustSize",
                    "CC_Mid_Posterior_ClustSize",
                    "CC_Central_ClustSize",
                    "CC_Mid_Anterior_ClustSize",
                    "CC_Anterior_ClustSize",
                    "ctx_lh_unknown_ClustSize",
                    "ctx_lh_bankssts_ClustSize",
                    "ctx_lh_caudalanteriorcingulate_ClustSize",
                    "ctx_lh_caudalmiddlefrontal_ClustSize",
                    "ctx_lh_cuneus_ClustSize",
                    "ctx_lh_entorhinal_ClustSize",
                    "ctx_lh_fusiform_ClustSize",
                    "ctx_lh_inferiorparietal_ClustSize",
                    "ctx_lh_inferiortemporal_ClustSize",
                    "ctx_lh_isthmuscingulate_ClustSize",
                    "ctx_lh_lateraloccipital_ClustSize",
                    "ctx_lh_lateralorbitofrontal_ClustSize",
                    "ctx_lh_lingual_ClustSize",
                    "ctx_lh_medialorbitofrontal_ClustSize",
                    "ctx_lh_middletemporal_ClustSize",
                    "ctx_lh_parahippocampal_ClustSize",
                    "ctx_lh_paracentral_ClustSize",
                    "ctx_lh_parsopercularis_ClustSize",
                    "ctx_lh_parsorbitalis_ClustSize",
                    "ctx_lh_parstriangularis_ClustSize",
                    "ctx_lh_pericalcarine_ClustSize",
                    "ctx_lh_postcentral_ClustSize",
                    "ctx_lh_posteriorcingulate_ClustSize",
                    "ctx_lh_precentral_ClustSize",
                    "ctx_lh_precuneus_ClustSize",
                    "ctx_lh_rostralanteriorcingulate_ClustSize",
                    "ctx_lh_rostralmiddlefrontal_ClustSize",
                    "ctx_lh_superiorfrontal_ClustSize",
                    "ctx_lh_superiorparietal_ClustSize",
                    "ctx_lh_superiortemporal_ClustSize",
                    "ctx_lh_supramarginal_ClustSize",
                    "ctx_lh_frontalpole_ClustSize",
                    "ctx_lh_temporalpole_ClustSize",
                    "ctx_lh_transversetemporal_ClustSize",
                    "ctx_lh_insula_ClustSize",
                    "ctx_rh_unknown_ClustSize",
                    "ctx_rh_bankssts_ClustSize",
                    "ctx_rh_caudalanteriorcingulate_ClustSize",
                    "ctx_rh_caudalmiddlefrontal_ClustSize",
                    "ctx_rh_cuneus_ClustSize",
                    "ctx_rh_entorhinal_ClustSize",
                    "ctx_rh_fusiform_ClustSize",
                    "ctx_rh_inferiorparietal_ClustSize",
                    "ctx_rh_inferiortemporal_ClustSize",
                    "ctx_rh_isthmuscingulate_ClustSize",
                    "ctx_rh_lateraloccipital_ClustSize",
                    "ctx_rh_lateralorbitofrontal_ClustSize",
                    "ctx_rh_lingual_ClustSize",
                    "ctx_rh_medialorbitofrontal_ClustSize",
                    "ctx_rh_middletemporal_ClustSize",
                    "ctx_rh_parahippocampal_ClustSize",
                    "ctx_rh_paracentral_ClustSize",
                    "ctx_rh_parsopercularis_ClustSize",
                    "ctx_rh_parsorbitalis_ClustSize",
                    "ctx_rh_parstriangularis_ClustSize",
                    "ctx_rh_pericalcarine_ClustSize",
                    "ctx_rh_postcentral_ClustSize",
                    "ctx_rh_posteriorcingulate_ClustSize",
                    "ctx_rh_precentral_ClustSize",
                    "ctx_rh_precuneus_ClustSize",
                    "ctx_rh_rostralanteriorcingulate_ClustSize",
                    "ctx_rh_rostralmiddlefrontal_ClustSize",
                    "ctx_rh_superiorfrontal_ClustSize",
                    "ctx_rh_superiorparietal_ClustSize",
                    "ctx_rh_superiortemporal_ClustSize",
                    "ctx_rh_supramarginal_ClustSize",
                    "ctx_rh_frontalpole_ClustSize",
                    "ctx_rh_temporalpole_ClustSize",
                    "ctx_rh_transversetemporal_ClustSize",
                    "ctx_rh_insula_ClustSize",
                ]
                self.qrep_dat[tracer] = self.qrep_dat[tracer][cols_in_order]

                # Save the output dataframe
                output_file = op.join(
                    self.paths["qreport"],
                    f"LEADS-PETCore-quarterly-report_{self.report_period}_{tracer.upper()}-ROI-means_{TODAY}.csv",
                )
                if save_output:
                    if overwrite or not op.isfile(output_file):
                        self.qrep_dat[tracer].to_csv(output_file, index=False)
                        print(f"Saved {output_file}")

            match tracer:
                case "fbb":
                    format_fbb(self)
                case "fdg":
                    format_fdg(self)
                case "ftp":
                    format_ftp(self)
                case _:
                    raise ValueError(f"Tracer {tracer} not recognized")

        # Prepare the final dataframes for each tracer
        self.qrep_dat = {}
        for tracer in self.tracers:
            # Make a copy of the PET index dataframe that we will use
            # to structure the quarterly report data
            self.qrep_dat[tracer] = self.pet_idx[tracer].copy()
            print("-" * 80)
            print(f"Building quarterly report for {tracer.upper()}")

            # Merge reference region scaling factors into the main dataframe
            self.qrep_dat[tracer] = self.qrep_dat[tracer].merge(
                self.ref_region_dat[tracer], on=["subject_id", "pet_date"], how="left"
            )
            if tracer == "fbb":
                # Merge amyloid visual read info into the main dataframe
                self.qrep_dat[tracer] = self.qrep_dat[tracer].merge(
                    self.screening, on="subject_id", how="left"
                )
                # Merge Centiloid data into the main dataframe
                self.qrep_dat[tracer] = self.qrep_dat[tracer].merge(
                    self.centiloid_dat, on=["subject_id", "pet_date"], how="left"
                )
            else:
                # Merge cohort assignment into the main dataframe
                self.qrep_dat[tracer] = self.qrep_dat[tracer].merge(
                    self.screening[["subject_id", "CohortAssgn"]],
                    on="subject_id",
                    how="left",
                )

            # Merge ROI mean SUVR and volume data into the main dataframe
            self.qrep_dat[tracer] = self.qrep_dat[tracer].merge(
                self.roi_dat[tracer], on=["subject_id", "pet_date"], how="left"
            )

            # Convert PET date to datetime
            self.qrep_dat[tracer]["pet_date"] = pd.to_datetime(
                self.qrep_dat[tracer]["pet_date"]
            )

            # Sort the data by subject_id and PET date
            self.qrep_dat[tracer] = (
                self.qrep_dat[tracer]
                .sort_values(["subject_id", "pet_date"])
                .reset_index(drop=True)
            )

            # Remove scans from the current quarter onward
            print(
                f"Checking for {tracer.upper()} scans after {self.stop_date.strftime('%Y-%m-%d')}"
            )
            scans_past_stop_date = self.qrep_dat[tracer].query(
                "(pet_date > @self.stop_date)"
            )
            if len(scans_past_stop_date) > 0:
                self.qrep_dat[tracer] = self.qrep_dat[tracer].query(
                    "(pet_date <= @self.stop_date)"
                )
                print(
                    "- Removing {:,} {} scans from {:,} subjects acquired after the stop date".format(
                        len(scans_past_stop_date),
                        tracer.upper(),
                        scans_past_stop_date["subject_id"].nunique(),
                    )
                )
                print(
                    "    "
                    + "\n    ".join(
                        scans_past_stop_date[["subject_id", "CohortAssgn", "pet_date"]]
                        .to_markdown(index=False, tablefmt="rst")
                        .split("\n")
                    )
                )
            else:
                print(f"- No {tracer.upper()} scans after the stop date")
            print()

            # Remove scans not in the eligibility list
            print("Checking eligibility list")
            ineligible = self.qrep_dat[tracer].query(
                "subject_id not in @self.eligible_subjects"
            )
            if len(ineligible) > 0:
                self.qrep_dat[tracer] = self.qrep_dat[tracer].query(
                    "subject_id in @self.eligible_subjects"
                )
                print(
                    "- Removing {} {} scans from {} subjects not in the eligibility list".format(
                        len(ineligible),
                        tracer.upper(),
                        ineligible["subject_id"].nunique(),
                    )
                )
                print(
                    "Ineligible subjects:",
                    ineligible["subject_id"].unique(),
                    sep="\n",
                )
            else:
                print("- No scans from ineligible subjects")
            print()

            # Convert PET date back to string
            self.qrep_dat[tracer]["pet_date"] = self.qrep_dat[tracer][
                "pet_date"
            ].dt.strftime("%Y-%m-%d")

            # Format the tracer-specific column data
            format_tracer(self, tracer)

            # Print final dataframe shape
            n_scans = len(self.qrep_dat[tracer])
            n_subjs = self.qrep_dat[tracer]["ID"].nunique()
            print(
                "",
                f"Final {tracer.upper()} quarterly report dataframe for {self.report_period}:",
                f"- Retained {n_scans:,} {tracer.upper()} scans from {n_subjs:,} subjects",
                f"- Dataframe shape: {self.qrep_dat[tracer].shape}",
                sep="\n",
                end="\n" * 2,
            )


def _parse_args():
    """Parse and return command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Save LEADS PET quarterly report files\n\n"
            + "Overview\n--------\n"
            + "This program aggregates ROI extraction files in LEADS processed PET scan\n"
            + "directories, creating 3 CSV files (one each for FBB, FTP, and FDG) that are\n"
            + "uploaded to LONI every quarter for other LEADS investigators to use.\n\n"
            + "Processed PET scans are filtered such that only scans that passed quality by\n"
            + "checks by both the University of Michigan (assessing PET reconstruction\n"
            + "quality) and UCSF (assessing MRI-based PET processing pipeline quality) are\n"
            + "retained in the final extraction files.\n\n"
            + "Column name and datatype formats for the extraction files is rigid, and any\n"
            + "deviations from the norm are likely to cause errors at the time we attempt\n"
            + "to upload these files to LONI. Thus, this script is designed to format the\n"
            + "extraction spreadsheets in exactly the manner expected by LONI."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        exit_on_error=False,
    )
    parser.add_argument(
        "-d",
        "--report_date",
        default=None,
        help=(
            "Target date for the quarterly report, formatted like YYYY-MM-DD.\n"
            "This is used to infer the report period and should be entered as a string\n"
            "like YYYY-MM-DD. You can also just directly enter the report period as a\n"
            "string like 'YYYY-QN', where N is an integer 1-4.\n\n"
            "Data will be included up to the end of the PREVIOUS quarter.\n"
            "For example, if you are compiling a report for 2025Q1, data will include\n"
            "all eligible PET scans through Dec. 31, 2024.\n\n"
            "By default, the quarterly report is generated for the current quarter,\n"
            "with data going up to the end of the previous quarter."
        ),
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
        help="Generate quarterly report dataframes but do not save them as CSVs",
    )
    parser.add_argument(
        "--no-overwrite",
        action="store_false",
        dest="overwrite",
        help="Do not overwrite existing quarterly report CSV files",
    )

    return parser.parse_args()


if __name__ == "__main__":
    start_time = time.time()

    # Get command line arguments.
    args = _parse_args()

    # Instantiate the QReport class and run the main function.
    qr = QReport(args.proj_dir, args.report_date)
    qr.compile(args.save, args.overwrite)

    # Report elapsed time.
    elapsed = time.time() - start_time
    minutes, seconds = divmod(elapsed, 60)
    print(f"Elapsed time: {int(minutes)} min, {int(seconds)} s", end="\n" * 2)

    # Exit successfully
    sys.exit(0)
