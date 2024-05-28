#!/usr/bin/env python

"""
Move scans from newdata to raw, keeping file hierarchies intact
"""


import argparse
import os
import os.path as op
import shutil
import sys
from glob import glob


def move_newdata_to_raw(newdata_dir, raw_dir, overwrite=False, wipe_newdata=True):
    """Move scans from newdata to raw, keeping file hierarchies intact

    Parameters
    ----------
    newdata_dir : str
        The directory containing the new scan data. Must format like:
            <newdata_dir>/<subj>/<...>/<nifti_files>
    raw_dir : str
        The directory to move the new scan data to. Structure after
        moving will be:
            <raw_dir>/<subj>/<...>/<nifti_files>
    overwrite : bool
        If True, overwrite existing scan directories in raw_dir with
        directories from newdata_dir. If False, skip existing
        directories.
    wipe_newdata : bool
        If True, remove all files and folders from newdata_dir after
        moving everything eligible to be moved to raw_dir.

    Returns
    -------
    None
    """

    def do_cleanup():
        """Remove all files and folders from newdata."""
        print(f"  * Cleaning up {newdata_dir}")
        for file in os.listdir(newdata_dir):
            filepath = op.join(newdata_dir, file)
            if op.isdir(filepath):
                shutil.rmtree(filepath)
            else:
                os.remove(filepath)

    # Ensure the base directory paths are absolute and normalized
    newdata_dir = op.abspath(newdata_dir)
    raw_dir = op.abspath(raw_dir)

    # Find all niftis in newdata
    check_exts = (".nii", ".nii.gz")
    glob_files = []
    for ext in check_exts:
        glob_files.extend(glob(op.join(newdata_dir, f"**/*{ext}"), recursive=True))

    # If no niftis are found, print a message and return
    if len(glob_files) == 0:
        print(f"  * No niftis found in {newdata_dir}")
        if wipe_newdata:
            do_cleanup()
        print("")
        return

    # Find all unique nifti-containing directories in newdata
    source_dirs = set([op.dirname(f) for f in glob_files])
    print(f"  * Found {len(source_dirs)} directories with niftis")

    skip_dirs = ["LEADS", "ADNI", "IDEAS", "SCAN", "PAD"]
    count = 0
    for source_dir in source_dirs:
        # Create a matching file hierarchy in raw as in newdata
        scan_path = op.relpath(source_dir, newdata_dir)
        scan_path_parts = scan_path.split("/")
        if scan_path_parts[0].upper() in skip_dirs:
            scan_path = op.join(*scan_path_parts[1:])
        target_dir = op.join(raw_dir, scan_path)

        # Check if the target directory exists
        if op.exists(target_dir):
            # If overwrite is True, remove the existing directory
            if overwrite:
                print(f"  * Overwriting existing raw directory: {target_dir}")
                shutil.rmtree(target_dir)
            else:
                print(f"  * Skipping existing raw directory: {target_dir}")
                continue

        # Create the necessary directory structure, then copy source to
        # target
        os.makedirs(op.dirname(target_dir), exist_ok=True)
        shutil.move(source_dir, target_dir)
        print(f"  * Moved {source_dir} to {target_dir}")
        count += 1

    # Print how many scans we moved over
    print(f"  * Done; {count} scan directories moved to raw")

    # Clean up empty directories in newdata
    if wipe_newdata:
        do_cleanup()


def _parse_args():
    """Parse and return command line arguments."""
    parser = argparse.ArgumentParser(
        description="Move scans from newdata to raw, keeping file hierarchies intact",
        formatter_class=argparse.RawTextHelpFormatter,
        exit_on_error=False,
    )
    parser.add_argument(
        "--newdata",
        type=str,
        default="/mnt/coredata/processing/leads/data/newdata",
        help=(
            "Path to <newdata> directory (where you put unprocessed PET/MRI scan\n"
            + "directories that you want to organize in <raw>; can be named something\n"
            + "other than 'newdata'\n"
            + "(default: %(default)s)"
        ),
    )
    parser.add_argument(
        "--raw",
        type=str,
        default="/mnt/coredata/processing/leads/data/raw",
        help=(
            "Path to <raw> directory (contains all unprocessed PET/MRI scan directories\n"
            + "for a project; this is where newdata directories will be moved to; can be\n"
            + "named something other than 'raw'\n"
            + "(default: %(default)s)"
        ),
    )
    parser.add_argument(
        "-o", "--overwrite", action="store_true", help="Overwrite existing files in raw"
    )
    parser.add_argument(
        "--no-wipe-newdata",
        action="store_false",
        dest="wipe_newdata",
        help="Don't wipe 'newdata' after moving scans to 'raw'",
    )

    # Parse the command line arguments
    args = parser.parse_args()
    if (len(sys.argv) == 1) and not all(op.isdir(args.newdata), op.isdir(args.raw)):
        parser.print_help()
        sys.exit()
    return args


if __name__ == "__main__":
    # Get command line arguments.
    args = _parse_args()

    # Move newdata to raw
    move_newdata_to_raw(
        newdata_dir=args.newdata,
        raw_dir=args.raw,
        overwrite=args.overwrite,
        wipe_newdata=args.wipe_newdata,
    )

    # Exit successfully
    sys.exit(0)
