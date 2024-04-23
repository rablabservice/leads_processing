#!/usr/bin/env python
import argparse
import os
import os.path as op
import shutil
import sys
from glob import glob
import pandas as pd


def create_pet_proc_dirs(
    raw_scans=None,
    log_dir="/mnt/coredata/processing/leads/metadata/log",
    overwrite=False,
    verbose=True,
):
    """Create processed PET directories and link to associated MRIs.

    For each scan that needs to be processed, there must already be:
    1. A raw PET file in .nii format
    2. A processed MRI directory that will be linked to

    This function then does the following for each scan:
    1. Creates new processed PET directory
    2. Copies raw PET file to processed PET directory and renames it
    3. Creates symbolic link from processed PET directory to the
       associated MRI directory

    Parameters
    ----------
    raw_scans : DataFrame or str, optional
        A pandas DataFrame or path to a CSV file with columns
        'raw_petf', 'scan_tag', 'mri_dir', and 'proc_pet_dir', which
        hold paths to raw PET .nii files, scan tags
        ("<subj>_<tracer>_<scan_date>"), processed MRI directories that
        will be used to process each PET scan, and target directories
        for processed PET data that will be created by this function,
        respectively. If None, the most recently saved raw_scans CSV
        file in log_dir will be used.
    overwrite : bool, optional
        If True, overwrite existing processed PET directories if they
        exist. If False, skip scans with existing processed PET
        directories.
    verbose : bool, optional
        If True, output status messages during execution.

    Returns
    -------
    None
    """
    # Print the welcome message
    if verbose:
        title = f"\nCreating processed PET directories"
        print(title, "-" * len(title), sep="\n")

    # Load the most recently saved raw_scans spreadsheet if not provided
    if raw_scans is None:
        raw_scansf = glob_sort_mtime(op.join(log_dir, "raw_pet_scans_*.csv"))[0]
        if verbose:
            print(f"- Reading {raw_scansf}")
        raw_scans = pd.read_csv(raw_scansf)
    elif isinstance(raw_scans, str):
        if verbose:
            print(f"- Reading {raw_scans}")
        raw_scans = pd.read_csv(raw_scans)

    # Filter scans that need to be processed
    raw_scans = raw_scans.query("(need_to_process==True)").reset_index(drop=True)
    if verbose:
        print(f"- {raw_scans.shape[0]} scans to process")

    # Loop over each scan and do directory setup
    for idx, scan in raw_scans.iterrows():
        # Make sure the raw PET file exists
        if not op.isfile(scan["raw_petf"]):
            if verbose:
                print(
                    f"- Skipping {scan['scan_tag']} due to missing raw PET file: {scan['raw_petf']}"
                )
            continue
        elif not scan["raw_petf"].endswith(".nii"):
            if verbose:
                print(
                    f"- Skipping {scan['scan_tag']} as raw PET file does not end in .nii: {scan['raw_petf']}"
                )
            continue
        # Make sure the MRI directory exists
        if not op.isdir(scan["mri_dir"]):
            if verbose:
                print(
                    f"- Skipping {scan['scan_tag']} due to missing MRI directory: {scan['mri_dir']}"
                )
            continue
        # Remove existing processed PET directories if overwrite is True
        if op.isdir(scan["proc_pet_dir"]):
            if overwrite:
                if verbose:
                    print(
                        f"- Removing existing directory and its contents: {scan['proc_pet_dir']}"
                    )
                shutil.rmtree(scan["proc_pet_dir"])
            else:
                # if verbose:
                #     print(
                #         f"- Skipping {scan['scan_tag']} due to existing processed PET directory: {scan['proc_pet_dir']}"
                #     )
                continue

        # Create the processed PET directory
        os.makedirs(scan["proc_pet_dir"])

        # Copy the raw PET file to the processed PET directory
        infile = scan["raw_petf"]
        outfile = op.join(scan["proc_pet_dir"], f"{scan['scan_tag']}.nii")
        shutil.copy(infile, outfile)

        # Create a symlink to the processed MRI directory
        link_src = scan["mri_dir"]
        link_dst = op.join(scan["proc_pet_dir"], "mri")
        os.symlink(link_src, link_dst)

    if verbose:
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
        description="Create processed PET directories and link to associated MRIs",
        formatter_class=argparse.RawTextHelpFormatter,
        exit_on_error=False,
    )
    parser.add_argument(
        "--raw_scans",
        type=str,
        help=(
            "Path to a CSV file with columns 'raw_petf', 'scan_tag', 'mri_dir', and\n"
            + "'proc_pet_dir', which hold paths to raw PET .nii files, scan tags\n"
            + "('<subj>_<tracer>_<scan_date>'), processed MRI directories that will\n"
            + "be used to process each PET scan, and target directories for processed\n"
            + "PET data that will be created by this function, respectively. If None,\n"
            + "the most recently saved raw_scans CSV file in log_dir will be used."
        ),
    )
    parser.add_argument(
        "--log_dir",
        type=str,
        default="/mnt/coredata/processing/leads/metadata/log",
        help=(
            "Path to directory where raw_scans CSV files are saved. Default is\n"
            + "/mnt/coredata/processing/leads/metadata/log"
        ),
    )
    parser.add_argument(
        "-o", "--overwrite", action="store_true", help="Overwrite existing files in raw"
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Run without printing output"
    )

    # Parse the command line arguments
    return parser.parse_args()


if __name__ == "__main__":
    # Get command line arguments.
    args = _parse_args()

    # Format CL args
    verbose = not args.quiet

    # Run the function
    create_pet_proc_dirs(
        raw_scans=args.raw_scans,
        log_dir=args.log_dir,
        overwrite=args.overwrite,
        verbose=verbose,
    )

    sys.exit(0)
