#!/usr/bin/env python

"""
Create processed MRI and PET directories and link PET scans to their associated MRIs
"""


import argparse
import os
import os.path as op
import shutil
import sys
from glob import glob

import pandas as pd


def setup_processed_scan_dirs(
    raw_mris=None,
    raw_pets=None,
    log_dir="/mnt/coredata/processing/leads/metadata/scans_to_process",
    overwrite=False,
):
    """Create processed MRI and PET directories and link PET scans to their associated MRIs

    Before this function is run, there must already be:
    1. A CSV file that lists the raw MRI scans to be processed
    2. A CSV file that lists the raw PET scans to be processed
    3. For each MRI to be processed, a raw MRI image with .nii extension
    4. For each PET scan to be processed, a raw PET image with .nii extension

    This function then does the following for each scan that is
    scheduled to be processed:
    1. Creates new processed MRI directory
    2. Symlinks from processed to raw MRI directory
    3. Creates new processed PET directory
    4. Symlinks from processed to raw PET directory
    5. Copies raw PET image to the processed PET directory
    6. Symlinks from processed PET to associated MRI directory

    Parameters
    ----------
    raw_mris : DataFrame or str, optional
        A pandas DataFrame or path to a CSV file with columns
        'subj', 'mri_date', 'mri_raw_niif', 'mri_proc_dir',
        'mri_image_id', and 'need_to_process'. If need_to_process==0 for
        a given MRI, it will be skipped. If this argument is None, the
        most recently saved Raw_MRI_Scan_Index_*.csv file in log_dir
        will be used
    raw_pets : DataFrame or str, optional
        A pandas DataFrame or path to a CSV file with columns
        'subj', 'tracer', 'pet_date', 'pet_raw_niif', 'pet_proc_dir',
        'mri_image_id', and 'need_to_process'. If need_to_process==0 for
        a given PET scan, it will be skipped. If this argument is None,
        the most recently saved Raw_PET_Scan_Index_*.csv file in log_dir
        will be used
    log_dir : str, optional
        Path to the directory containing the log files that list the
        raw MRI and PET scans to be processed
    overwrite : bool, optional
        If True, remove and recreate any existing processed MRI and PET
        directories that are scheduled to be processed. Note that this
        argument will not affect which scans will be processed, which is
        determined by the 'need_to_process' column in the raw_mris and
        raw_pets DataFrames

    Returns
    -------
    None
    """
    # Load the most recently saved raw_mris spreadsheet if not provided
    if raw_mris is None:
        raw_mrisf = glob_sort_mtime(op.join(log_dir, "Raw_MRI_Scan_Index_*.csv"))[0]
    elif isinstance(raw_mris, str):
        raw_mrisf = raw_mris
    print(f"- Reading {raw_mrisf}")
    raw_mris = pd.read_csv(raw_mrisf)

    # Loop over each MRI to be processed and setup the processed
    # directory structure
    for _, scan in raw_mris.query("(need_to_process==1)").iterrows():
        # Print the scan tag
        mri_tag = f"{scan['subj']}_MRI-T1_{scan['mri_date']}"
        print(f"  * {mri_tag}")

        # Make sure the raw MRI file exists
        if not op.isfile(scan["mri_raw_niif"]):
            print(f"    - {scan['mri_raw_niif']} does not exist, skipping scan")
            continue

        # Create the processed MRI directory
        if op.isdir(scan["mri_proc_dir"]):
            if overwrite:
                shutil.rmtree(scan["mri_proc_dir"])
                print(f"    $ rm -rf {scan['mri_proc_dir']}")
                os.makedirs(scan["mri_proc_dir"])
                print(f"    $ mkdir {scan['mri_proc_dir']}")
        else:
            os.makedirs(scan["mri_proc_dir"])
        print(f"    $ mkdir {scan['mri_proc_dir']}")

        # Create a symlink to the raw MRI directory
        link_src = op.dirname(scan["mri_raw_niif"])
        link_dst = op.join(scan["mri_proc_dir"], "raw")
        if not op.islink(link_dst):
            os.symlink(link_src, link_dst)
            print(f"    $ ln -s {link_src} {link_dst}")

    # Load the most recently saved raw_pets spreadsheet if not provided
    if raw_pets is None:
        raw_petsf = glob_sort_mtime(op.join(log_dir, "Raw_PET_Scan_Index_*.csv"))[0]
    elif isinstance(raw_pets, str):
        raw_petsf = raw_pets
    print(f"- Reading {raw_petsf}")
    raw_pets = pd.read_csv(raw_petsf)

    # Loop over each PET scan to be processed and setup the processed
    # directory structure
    for _, scan in raw_pets.query("(need_to_process==1)").iterrows():
        # Print the scan tag
        pet_tag = f"{scan['subj']}_{scan['tracer']}_{scan['pet_date']}"
        print(f"  * {pet_tag}")

        # Make sure the raw PET file exists
        if not op.isfile(scan["pet_raw_niif"]):
            print(f"    - {scan['pet_raw_niif']} does not exist, skipping scan")
            continue

        # Make sure the raw PET file ends in .nii
        if not scan["pet_raw_niif"].endswith(".nii"):
            print(f"    - {scan['pet_raw_niif']} does not end in .nii, skipping scan")
            continue

        # Make sure the processed MRI directory exists
        mri_image_id = scan["mri_image_id"]
        result = raw_mris.loc[raw_mris["mri_image_id"] == mri_image_id, "mri_proc_dir"]
        if result.empty:
            print(
                f"    - {scan['pet_raw_niif']} is missing an accompanying processed MRI directory, skipping scan"
            )
            continue
        else:
            mri_proc_dir = result.iloc[0]

        # Create the processed PET directory
        if op.isdir(scan["pet_proc_dir"]):
            if overwrite:
                shutil.rmtree(scan["pet_proc_dir"])
                print(f"    $ rm -rf {scan['pet_proc_dir']}")
                os.makedirs(scan["pet_proc_dir"])
                print(f"    $ mkdir {scan['pet_proc_dir']}")
        else:
            os.makedirs(scan["pet_proc_dir"])
            print(f"    $ mkdir {scan['pet_proc_dir']}")

        # Create a symlink to the raw PET directory
        link_src = op.dirname(scan["pet_raw_niif"])
        link_dst = op.join(scan["pet_proc_dir"], "raw")
        if not op.islink(link_dst):
            os.symlink(link_src, link_dst)
            print(f"    $ ln -s {link_src} {link_dst}")

        # Copy the raw PET file to the processed PET directory
        infile = scan["pet_raw_niif"]
        outfile = op.join(scan["pet_proc_dir"], f"{pet_tag}.nii")
        if not op.isfile(outfile):
            shutil.copy(infile, outfile)
            print(f"    $ cp {infile} {outfile}")

        # Create a symlink to the processed MRI directory
        link_src = mri_proc_dir
        link_dst = op.join(scan["pet_proc_dir"], "mri")
        if not op.islink(link_dst):
            os.symlink(link_src, link_dst)
            print(f"    $ ln -s {link_src} {link_dst}")

    print("")


def glob_sort_mtime(pattern):
    """Return files matching pattern in most recent modified order.

    Returns
    -------
    files : list of str
        List of files matching pattern, sorted by most recent modified
        (files[0] is the most recently modified).
    """
    files = sorted(glob(pattern), key=op.getmtime, reverse=True)
    return files


def _parse_args():
    """Parse and return command line arguments."""
    parser = argparse.ArgumentParser(
        description="Create processed MRI and PET directories and link PET scans to their associated MRIs",
        formatter_class=argparse.RawTextHelpFormatter,
        exit_on_error=False,
    )
    parser.add_argument(
        "--raw_mris",
        type=str,
        help=(
            "Path to the raw MRI scan CSV file that lists which MRIs will be processed\n"
            + "If None, the most recently saved Raw_MRI_Scan_Index_*.csv in log_dir will be used"
        ),
    )
    parser.add_argument(
        "--raw_pets",
        type=str,
        help=(
            "Path to the raw PET scan CSV file that lists which PET scans will be processed\n"
            + "If None, the most recently saved Raw_PET_Scan_Index_*.csv file in log_dir will be used"
        ),
    )
    parser.add_argument(
        "--log_dir",
        type=str,
        default="/mnt/coredata/processing/leads/metadata/scans_to_process",
        help="Path to the log directory (default: %(default)s)",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Overwrite existing files in processed",
    )

    # Parse the command line arguments
    return parser.parse_args()


if __name__ == "__main__":
    # Get command line arguments.
    args = _parse_args()

    # Call the main function
    setup_processed_scan_dirs(
        raw_mris=args.raw_mris,
        raw_pets=args.raw_pets,
        log_dir=args.log_dir,
        overwrite=args.overwrite,
    )

    # Exit successfully
    sys.exit(0)
