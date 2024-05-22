#!/usr/bin/env python

"""
Unzip all .zip files in a directory.
"""
import argparse
import os.path as op
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from glob import glob


def fast_unzip_dir(source_dir):
    """
    Unzips all .zip files in source_dir using multiprocessing if needed.

    Parameters
    ----------
    source_dir : str
        Path to the directory with .zip files that you wish to unzip

    Returns
    -------
    None
    """
    # Ensure the source directory exists
    if not op.isdir(source_dir):
        raise ValueError(f"{source_dir} does not exist or is not a directory")

    # Find all .zip files in source_dir
    zip_files = glob(op.join(source_dir, "*.zip"))

    # Unzip the files, parallelizing if there are multiple .zip files
    if len(zip_files) == 0:
        print(f"- No files to unzip in {source_dir}")
    elif len(zip_files) == 1:
        print(f"- Scanning {source_dir}, found 1 file to unzip...")
        unzip_file(zip_files[0], source_dir)
    else:
        print(f"- Scanning {source_dir}, found {len(zip_files)} files to unzip...")
        max_workers = min(16, len(zip_files))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_zip = {
                executor.submit(unzip_file, zipf, source_dir): zipf
                for zipf in zip_files
            }
            for future in as_completed(future_to_zip):
                zipf = future_to_zip[future]
                try:
                    future.result()
                except Exception as e:
                    print(f"Error unzipping {zipf}: {e}")


def unzip_file(zipf, target_dir=None):
    """Unzip a .zip file.

    Parameters
    ----------
    zipf : str
        Path to the .zip file to be unzipped

    Returns
    -------
    None
    """
    with zipfile.ZipFile(zipf) as zf:
        try:
            zf.extractall(target_dir)
            msg = f"  * Unzipped {zipf}"
            if target_dir is not None:
                msg += f" to {target_dir}"
            print(msg)
        except zipfile.BadZipFile:
            print(f"ERROR: {zipf} is not a valid zip file")


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Unzip all .zip files in a directory")
    parser.add_argument(
        "-d",
        "--source_dir",
        default=".",
        help="Directory containing the .zip files to be unzipped",
    )
    return parser.parse_args()


# If this script is run from the command line, parse arguments and run
if __name__ == "__main__":
    args = _parse_args()
    fast_unzip_dir(args.source_dir)
