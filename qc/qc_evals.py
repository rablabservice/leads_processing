#!/usr/bin/env python

"""
Create and manage QC eval files.
"""


import argparse
import os
import shutil
import sys
import os.path as op

import numpy as np
import pandas as pd

utils_dir = op.join(op.dirname(__file__), "..", "utils")
if utils_dir not in sys.path:
    sys.path.append(utils_dir)
import utilities as uts

setup_dir = op.join(op.dirname(__file__), "..", "setup")
if setup_dir not in sys.path:
    sys.path.append(setup_dir)
import select_scans_to_process as sstp


def get_qc_eval_fields(scan_type):
    """Return the QC eval columns for a given scan type."""
    qc_fields = {
        "MRI-T1": [
            "subj",
            "scan_date",
            "rater",
            "manual_freesurfer_edits",
            "native_nu_rating",
            "aparc_rating",
            "spm_seg_ok",
            "affine_nu_ok",
            "warped_nu_ok",
            "notes",
        ],
        "FBB": [
            "subj",
            "scan_date",
            "rater",
            "native_pet_ok",
            "pet_to_mri_coreg_ok",
            "wcbl_mask_ok",
            "eroded_wm_and_brainstem_masks_ok",
            "affine_pet_ok",
            "warped_pet_ok",
            "notes",
        ],
        "FDG": [
            "subj",
            "scan_date",
            "rater",
            "native_pet_ok",
            "pet_to_mri_coreg_ok",
            "pons_mask_ok",
            "affine_pet_ok",
            "warped_pet_ok",
            "notes",
        ],
        "FTP": [
            "subj",
            "scan_date",
            "rater",
            "native_pet_ok",
            "pet_to_mri_coreg_ok",
            "infcblgm_mask_ok",
            "erodedwm_mask_ok",
            "affine_pet_ok",
            "warped_pet_ok",
            "notes",
        ],
    }
    return qc_fields[scan_type]


def create_qc_eval_file(scan_dir):
    """Create a blank QC eval file.

    Parameters
    ----------
    scan_dir : str
        Path to the processed scan directory.
    """
    # Make sure the scan directory exists
    scan_dir = op.abspath(scan_dir)
    if not op.isdir(scan_dir):
        raise FileNotFoundError(scan_dir)

    # Get scan info
    scan_tag = uts.get_scan_tag(scan_dir)
    subj, scan_type, scan_date = uts.parse_scan_tag(scan_tag)

    # Create the QC eval directory
    qc_dir = op.join(scan_dir, "qc")
    if not op.isdir(qc_dir):
        os.makedirs(qc_dir)

    # Initialize the QC eval dataframe
    fields = pd.Index(get_qc_eval_fields(scan_type), name="field")
    values = [subj, scan_date] + [np.nan] * (len(fields) - 2)
    qc_eval = pd.Series(index=fields, data=values, name="value").reset_index()
    # qc_eval = pd.DataFrame(
    #     data=[subj, scan_date] + [np.nan] * (len(columns) - 2), columns=columns
    # )

    # Save the dateframe as a CSV file
    timestamp = uts.now()
    qc_eval_file = op.join(qc_dir, f"{scan_tag}_qc-eval_{timestamp}.csv")
    qc_eval.to_csv(qc_eval_file, index=False)
    print(f"Saved {qc_eval_file}")

    return qc_eval_file


def get_qc_eval_file(scan_dir):
    """Return the most recent QC eval file for a single scan."""
    scan_dir = op.abspath(scan_dir)
    qc_eval_files = uts.glob_sort(op.join(scan_dir, "qc", "*_qc-eval_*.csv"))
    if len(qc_eval_files) == 0:
        return None
    qc_eval_file = qc_eval_files[-1]
    return qc_eval_file


def check_if_qc_complete(scan_dir):
    """Determine if a scan has been QC'd.

    Check the most recently created QC eval file, and return True
    only if all required columns have been correctly filled out.
    """
    # Get the most recent QC eval file
    qc_eval_file = get_qc_eval_file(scan_dir)
    if qc_eval_file is None:
        return False

    # Check each column in qc_eval_file against its completion
    # criteria. Return True only if all criteria are passed.
    qc_complete = True
    qc_eval = pd.read_csv(qc_eval_file).set_index("field")["value"]
    for key, val in qc_eval.items():
        if key == "subj":
            if not isinstance(val, str) or len(val) == 0:
                qc_complete = False
                break
        elif key == "scan_date":
            if not isinstance(val, str):
                qc_complete = False
                break
            try:
                year, month, day = map(int, val.split("-"))
                if any(
                    [
                        year < 2000,
                        year > pd.Timestamp("today").year,
                        month < 1,
                        month > 12,
                        day < 1,
                        day > 31,
                    ]
                ):
                    qc_complete = False
                    break
            except ValueError:
                qc_complete = False
                break
        elif key == "rater":
            if not isinstance(val, str) or len(val) == 0:
                qc_complete = False
                break
        elif key == "manual_freesurfer_edits":
            try:
                val = int(val)
                if not val in (0, 1):
                    qc_complete = False
                    break
            except ValueError:
                qc_complete = False
                break
        elif key.endswith("_rating"):
            try:
                val = int(val)
                if not val in (0, 1, 2):
                    qc_complete = False
                    break
            except ValueError:
                qc_complete = False
                break
        elif key.endswith("_ok"):
            try:
                val = int(val)
                if not val in (0, 1):
                    qc_complete = False
                    break
            except ValueError:
                qc_complete = False
                break
    return qc_complete


def get_processed_scan_dirs(scan_type, proc_dir):
    """Return a list of all processed directories for the scan type."""
    proc_dir = op.abspath(proc_dir)
    scan_dirs = []
    with os.scandir(proc_dir) as subj_dirs:
        for subj_dir in subj_dirs:
            # Get all directories in the subject directory
            subj_scan_dirs = [
                entry.path
                for entry in os.scandir(subj_dir)
                if entry.is_dir() and entry.name.startswith(scan_type)
            ]
            scan_dirs.extend(subj_scan_dirs)
    return scan_dirs


def merge_completed_qc_eval_files(
    proj_dir="/mnt/coredata/processing/leads",
    scan_types=("MRI-T1", "FBB", "FDG", "FTP"),
    save_output=True,
):
    """Merge all completed QC eval files in the project directory.

    Save a new QC eval spreadsheet for each scan type, by merging
    all completed QC eval files in the processed scan directories.
    """
    # Define paths to the needed directories
    proj_dir = op.abspath(proj_dir)
    proc_dir = op.abspath(op.join(proj_dir, "data", "processed"))
    qc_dir = op.abspath(op.join(proj_dir, "metadata", "qc"))
    qc_archive = op.join(qc_dir, "archive")
    os.makedirs(qc_archive, exist_ok=True)

    timestamp = uts.now()
    qc_evals = {}
    for scan_type in scan_types:
        # Get all processed scan directories
        scan_dirs = get_processed_scan_dirs(scan_type, proc_dir)

        # Get a list of all completed QC eval files
        completed_qc_eval_files = [
            get_qc_eval_file(scan_dir)
            for scan_dir in scan_dirs
            if check_if_qc_complete(scan_dir)
        ]

        # Load the QC eval files into a single dataframe
        qc_evals[scan_type] = pd.concat(
            [pd.read_csv(f).set_index("field").T for f in completed_qc_eval_files],
            ignore_index=True,
        )

        # Save the merged QC eval dataframe
        if save_output:
            # Move existing QC eval files to the archive
            old_files = uts.glob_sort(
                op.join(qc_dir, f"processed_{scan_type}_qc-evals_*.csv")
            )
            for f in old_files:
                shutil.move(f, qc_archive)

            # Save the new QC eval file
            new_file = op.join(
                qc_dir, f"processed_{scan_type}_qc-evals_{timestamp}.csv"
            )
            qc_evals[scan_type].to_csv(new_file, index=False)
            print(f"Saved {new_file}")

    return qc_evals


def find_scans_with_incomplete_qc(
    proj_dir="/mnt/coredata/processing/leads",
    scan_types=("MRI-T1", "FBB", "FDG", "FTP"),
):
    """Find all scans with incomplete QC eval."""
    # Define paths to the needed directories
    proj_dir = op.abspath(proj_dir)
    proc_dir = op.abspath(op.join(proj_dir, "data", "processed"))

    scan_dirs = {}
    processed_scans = {}
    unprocessed_scans = {}
    qc_complete_scans = {}
    qc_incomplete_scans = {}
    n_scans = {}
    n_processed_scans = {}
    n_unprocessed_scans = {}
    n_qc_complete_scans = {}
    n_qc_incomplete_scans = {}
    for scan_type in scan_types:
        # Get all processed scan directories
        scan_dirs[scan_type] = get_processed_scan_dirs(scan_type, proc_dir)

        # Find processed vs. unprocessed scans
        if scan_type == "MRI-T1":
            processed_scans[scan_type] = [
                d for d in scan_dirs[scan_type] if sstp.check_if_mri_processed(d)
            ]
        else:
            processed_scans[scan_type] = [
                d for d in scan_dirs[scan_type] if sstp.check_if_pet_processed(d)
            ]
        unprocessed_scans[scan_type] = [
            d for d in scan_dirs[scan_type] if d not in processed_scans[scan_type]
        ]

        # Find processed scans with complete vs. incomplete QC evals
        qc_complete_scans[scan_type] = [
            d for d in processed_scans[scan_type] if check_if_qc_complete(d)
        ]
        qc_incomplete_scans[scan_type] = [
            d
            for d in processed_scans[scan_type]
            if d not in qc_complete_scans[scan_type]
        ]

        # Sort the incomplete scans by scan date
        idx = np.argsort(
            [op.basename(d).split("_")[1] for d in qc_incomplete_scans[scan_type]]
        )
        qc_incomplete_scans[scan_type] = np.array(qc_incomplete_scans[scan_type])[
            idx
        ].tolist()

        # Get counts of scans in each category
        n_scans[scan_type] = len(scan_dirs[scan_type])
        n_processed_scans[scan_type] = len(processed_scans[scan_type])
        n_unprocessed_scans[scan_type] = len(unprocessed_scans[scan_type])
        n_qc_complete_scans[scan_type] = len(qc_complete_scans[scan_type])
        n_qc_incomplete_scans[scan_type] = len(qc_incomplete_scans[scan_type])

        # Print scan counts in each category, ending with a list
        # of scans with incomplete QC evals
        print(f"\n{scan_type}", "-" * len(scan_type), sep="\n")
        print(f"- {n_scans[scan_type]:,} scans in {proc_dir}")
        print(f"- {n_processed_scans[scan_type]:,} scans have been fully processed")
        print(f"- {n_qc_complete_scans[scan_type]:,} processed scans have been QC'd")

        if n_qc_incomplete_scans[scan_type] > 0:
            print(
                f"- {n_qc_incomplete_scans[scan_type]:,} processed scans with incomplete QC"
            )
            uts.print_list(qc_incomplete_scans[scan_type])
        else:
            print("")

    return qc_incomplete_scans


def _parse_args():
    """Parse and return command line arguments."""
    parser = argparse.ArgumentParser(
        description="QC evaluation tools!",
        formatter_class=argparse.RawTextHelpFormatter,
        exit_on_error=False,
    )
    parser.add_argument(
        "task",
        type=str,
        choices=[
            "create",
            "check",
            "merge",
            "incomplete",
        ],
        help=(
            "Which QC evaluation task to perform:\n"
            + "  'create'     : Create a new QC eval file for one or more processed scans\n"
            + "  'check'      : Check if QC eval is complete for one or more processed scans\n"
            + "  'merge'      : Merge all completed QC eval files in the database into a\n"
            + "                 single spreadsheet for each scan type\n"
            + "  'incomplete' : Find all scans with incomplete QC eval\n"
        ),
    )
    parser.add_argument(
        "-s",
        "--scan_dirs",
        type=str,
        nargs="*",
        help=(
            "One or more scan directories that you want to create new QC eval files for\n"
            + "(if task is 'create'), or check if QC evals have been completed (if task is 'check')"
        ),
    )
    parser.add_argument(
        "--proj_dir",
        type=str,
        default="/mnt/coredata/processing/leads",
        help="Path to the top-level project directory. Default: %(default)s",
    )

    # Parse the command line arguments
    if (len(sys.argv) == 1) or (sys.argv[1] in ("-h", "--help")):
        parser.print_help()
        sys.exit()
    return parser.parse_args()


if __name__ == "__main__":
    # Get command line arguments
    args = _parse_args()

    # Run the requested task
    if args.task == "create":
        for scan_dir in args.scan_dirs:
            scan_dir = op.abspath(scan_dir)
            if not op.isdir(scan_dir):
                raise FileNotFoundError(scan_dir)
            create_qc_eval_file(scan_dir)
    elif args.task == "check":
        for scan_dir in args.scan_dirs:
            scan_dir = op.abspath(scan_dir)
            if not op.isdir(scan_dir):
                raise FileNotFoundError(scan_dir)
            scan_tag = uts.get_scan_tag(scan_dir)
            if check_if_qc_complete(scan_dir):
                print(f"{scan_tag} QC is complete")
            else:
                print(f"{scan_tag} QC is incomplete")
    elif args.task == "merge":
        if not op.isdir(args.proj_dir):
            raise FileNotFoundError(args.proj_dir)
        merge_completed_qc_eval_files(args.proj_dir)
        find_scans_with_incomplete_qc(args.proj_dir)
    elif args.task == "incomplete":
        if not op.isdir(args.proj_dir):
            raise FileNotFoundError(args.proj_dir)
        find_scans_with_incomplete_qc(args.proj_dir)
