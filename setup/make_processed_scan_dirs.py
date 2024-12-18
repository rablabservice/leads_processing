#!/usr/bin/env python

"""
Create processed MRI and PET directories and link PET scans to their associated MRIs
"""


import argparse
import os
import os.path as op
import shutil
import sys

import pandas as pd

utils_dir = op.join(op.dirname(__file__), "..", "utils")
if utils_dir not in sys.path:
    sys.path.append(utils_dir)
import utilities as uts


def make_processed_scan_dirs(
    scans_to_process_dir="/mnt/coredata/processing/leads/metadata/scans_to_process",
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
    5. Symlinks from processed PET to associated MRI directory

    Parameters
    ----------
    scans_to_process_dir : str, optional
        Path to the directory containing CSV files that list the MRI and
        PET scans in 'raw' that are scheduled to be processed
    overwrite : bool, optional
        If True, remove and recreate any existing processed MRI and PET
        directories that are scheduled to be processed. Note that this
        argument will not affect which scans will be processed, which is
        determined by the 'scheduled_for_processing' column in the
        raw_mris and raw_pets DataFrames

    Returns
    -------
    None
    """

    def print_tag(scan_tag, already_printed=False):
        """Print the scan tag if it hasn't already been printed"""
        if already_printed:
            return already_printed
        else:
            print(f"  * {scan_tag}")
            already_printed = True
            return already_printed

    # Load the most recently modified raw_MRI-T1_index* CSV file in
    # scans_to_process_dir
    try:
        raw_mrisf = uts.glob_sort_mtime(
            op.join(scans_to_process_dir, "raw_MRI-T1_index*.csv")
        )[0]
    except IndexError:
        print(f"ERROR: No raw_MRI-T1_index*.csv file found in {scans_to_process_dir}")
        sys.exit(1)
    raw_mris = pd.read_csv(raw_mrisf)

    # Load the most recently saved raw_PET_index* CSV file in
    # scans_to_process_dir
    try:
        raw_petsf = uts.glob_sort_mtime(
            op.join(scans_to_process_dir, "raw_PET_index*.csv")
        )[0]
    except IndexError:
        print(f"ERROR: No raw_PET_index_*.csv file found in {scans_to_process_dir}")
        sys.exit(1)
    raw_pets = pd.read_csv(raw_petsf)

    # Loop over each MRI to be processed and setup the processed
    # directory structure
    count_problems_mri = 0
    for _, scan in raw_mris.query("(scheduled_for_processing==1)").iterrows():
        # Print the scan tag
        already_printed = False
        mri_tag = f"{scan['subj']}_MRI-T1_{scan['mri_date']}"

        # Make sure the raw MRI file exists
        if not op.isfile(scan["mri_raw_niif"]):
            already_printed = print_tag(mri_tag, already_printed)
            print(f"!!  - {scan['mri_raw_niif']} DOES NOT EXIST; SKIPPING SCAN")
            count_problems_mri += 1
            continue

        # Create the processed MRI directory
        if op.isdir(scan["mri_proc_dir"]):
            if overwrite:
                if not already_printed:
                    print(f"  * {mri_tag}")
                    already_printed = True
                shutil.rmtree(scan["mri_proc_dir"])
                print(f"    $ rm -rf {scan['mri_proc_dir']}")
                os.makedirs(scan["mri_proc_dir"])
                print(f"    $ mkdir {scan['mri_proc_dir']}")
        else:
            already_printed = print_tag(mri_tag, already_printed)
            os.makedirs(scan["mri_proc_dir"])
            print(f"    $ mkdir {scan['mri_proc_dir']}")

        # Create a symlink to the raw MRI directory
        link_src = op.dirname(scan["mri_raw_niif"])
        link_dst = op.join(scan["mri_proc_dir"], "raw")
        if not op.islink(link_dst):
            already_printed = print_tag(mri_tag, already_printed)
            os.symlink(link_src, link_dst)
            print(f"    $ ln -s {link_src} {link_dst}")
        # If the symlink already exists, alert the user if it points to
        # an unexpected location
        else:
            if op.realpath(link_dst) != link_src:
                already_printed = print_tag(mri_tag, already_printed)
                print(
                    f"!!  - {link_dst} ALREADY EXISTS BUT POINTS TO {op.realpath(link_dst)}; EXPECTED LOCATION IS {link_src}"
                )
                count_problems_mri += 1

    # Loop over each PET scan to be processed and setup the processed
    # directory structure
    count_problems_pet = 0
    for _, scan in raw_pets.query("(scheduled_for_processing==1)").iterrows():
        # Print the scan tag
        already_printed = False
        pet_tag = f"{scan['subj']}_{scan['tracer']}_{scan['pet_date']}"

        # Make sure the raw PET file exists
        if not op.isfile(scan["pet_raw_niif"]):
            already_printed = print_tag(pet_tag, already_printed)
            print(f"!!  - {scan['pet_raw_niif']} DOES NOT EXIST; SKIPPING SCAN")
            count_problems_pet += 1
            continue

        # Make sure the raw PET file ends in .nii
        if not scan["pet_raw_niif"].endswith(".nii"):
            already_printed = print_tag(pet_tag, already_printed)
            print(f"!!  - {scan['pet_raw_niif']} DOES NOT END IN .nii; SKIPPING SCAN")
            count_problems_pet += 1
            continue

        # Make sure the processed MRI directory exists
        mri_image_id = scan["mri_image_id"]
        result = raw_mris.loc[raw_mris["mri_image_id"] == mri_image_id, "mri_proc_dir"]
        if result.empty:
            already_printed = print_tag(pet_tag, already_printed)
            print(
                f"!!  - {scan['pet_raw_niif']} IS MISSING AN ACCOMPANYING MRI; SKIPPING SCAN"
            )
            count_problems_pet += 1
            continue
        else:
            mri_proc_dir = result.iloc[0]

        # Create the processed PET directory
        if op.isdir(scan["pet_proc_dir"]):
            if overwrite:
                already_printed = print_tag(pet_tag, already_printed)
                shutil.rmtree(scan["pet_proc_dir"])
                print(f"    $ rm -rf {scan['pet_proc_dir']}")
                os.makedirs(scan["pet_proc_dir"])
                print(f"    $ mkdir {scan['pet_proc_dir']}")
        else:
            already_printed = print_tag(pet_tag, already_printed)
            os.makedirs(scan["pet_proc_dir"])
            print(f"    $ mkdir {scan['pet_proc_dir']}")

        # Create a symlink to the raw PET directory
        link_src = op.dirname(scan["pet_raw_niif"])
        link_dst = op.join(scan["pet_proc_dir"], "raw")
        if not op.islink(link_dst):
            already_printed = print_tag(pet_tag, already_printed)
            os.symlink(link_src, link_dst)
            print(f"    $ ln -s {link_src} {link_dst}")
        # If the symlink already exists, alert the user if it points to
        # an unexpected location
        else:
            if op.realpath(link_dst) != link_src:
                already_printed = print_tag(pet_tag, already_printed)
                print(
                    f"!!  - {link_dst} ALREADY EXISTS BUT POINTS TO {op.realpath(link_dst)}; EXPECTED LOCATION IS {link_src}"
                )
                count_problems_pet += 1

        # Create a symlink to the processed MRI directory
        link_src = mri_proc_dir
        link_dst = op.join(scan["pet_proc_dir"], "mri")
        if not op.islink(link_dst):
            already_printed = print_tag(pet_tag, already_printed)
            os.symlink(link_src, link_dst)
            print(f"    $ ln -s {link_src} {link_dst}")
        # If the symlink already exists, alert the user if it points to
        # an unexpected location
        else:
            if op.realpath(link_dst) != link_src:
                already_printed = print_tag(pet_tag, already_printed)
                print(
                    f"!!  - {link_dst} ALREADY EXISTS BUT POINTS TO {op.realpath(link_dst)}; EXPECTED LOCATION IS {link_src}"
                )
                count_problems_pet += 1

    # Notify the user of any problems detected
    print("")
    if count_problems_mri > 0:
        print(
            f"!! {count_problems_mri} PROBLEMS DETECTED WITH MRI SCANS; SEE ABOVE FOR DETAILS"
        )
    else:
        print(
            "  * No conflicts detected between existing data structure and latest MRI index"
        )

    if count_problems_pet > 0:
        print(
            f"!! PROBLEMS DETECTED WITH {count_problems_pet:,} PET SCANS; SEE ABOVE FOR DETAILS"
        )
    else:
        print(
            "  * No conflicts detected between existing data structure and latest PET index"
        )


def _parse_args():
    """Parse and return command line arguments."""
    parser = argparse.ArgumentParser(
        description="Create processed MRI and PET directories and link PET scans to their associated MRIs",
        formatter_class=argparse.RawTextHelpFormatter,
        exit_on_error=False,
    )
    parser.add_argument(
        "--scans_to_process_dir",
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
    make_processed_scan_dirs(
        scans_to_process_dir=args.scans_to_process_dir,
        overwrite=args.overwrite,
    )

    # Exit successfully
    sys.exit(0)
