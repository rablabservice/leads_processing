#!/usr/bin/env python

"""
Select MRI and PET scans to process and save the list to CSV files
"""


import argparse
import datetime
import os
import os.path as op
import re
import sys
from glob import glob

import numpy as np
import pandas as pd


PATHS = {}
SCAN_TYPES = None
TIMESTAMP = None


def set_globals():
    """Set module-level global variables"""
    global PATHS, SCAN_TYPES, TIMESTAMP

    root_dir = "/mnt/coredata"
    if not __file__.startswith(root_dir):
        raise ValueError(
            f"Could not parse the project path because {__file__} is not in {root_dir}"
        )

    # The project directory is assumed to be the 4th directory up from /
    # (e.g., /mnt/coredata/processing/leads)
    PATHS["proj"] = "/".join(__file__.split("/")[:5])

    # Set global paths and check that they exist
    PATHS["code"] = op.join(PATHS["proj"], "code")
    PATHS["config"] = op.join(PATHS["code"], "config")
    PATHS["metadata"] = op.join(PATHS["proj"], "metadata")
    PATHS["scans_to_process"] = op.join(PATHS["metadata"], "scans_to_process")
    PATHS["ssheets"] = op.join(PATHS["metadata"], "ssheets")
    PATHS["data"] = op.join(PATHS["proj"], "data")
    PATHS["newdata"] = op.join(PATHS["data"], "newdata")
    PATHS["raw"] = op.join(PATHS["data"], "raw")
    PATHS["processed"] = op.join(PATHS["data"], "processed")
    for k in PATHS:
        if not op.isdir(PATHS[k]):
            raise ValueError(f"Expected {PATHS[k]}, but this directory does not exist")

    # Define the SCAN_TYPES dict
    SCAN_TYPES = load_scan_typesf()

    # Set the timestamp
    TIMESTAMP = now()


def load_scan_typesf(scan_typesf=None):
    """Load scan types CSV and return a {name_in: name_out} dict."""
    if scan_typesf is None:
        scan_typesf = op.join(PATHS["config"], "scan_types_and_tracers.csv")
    scan_types = pd.read_csv(scan_typesf)
    scan_types["name_in"] = scan_types["name_in"].str.lower()
    scan_types = scan_types.drop_duplicates("name_in").dropna()
    scan_types = scan_types.set_index("name_in")["name_out"].to_dict()
    return scan_types


def now():
    """Return the current date and time down to seconds."""
    return datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")


set_globals()


def main(overwrite, process_unused_mris):
    """Save CSV files for MRI and PET scans in the raw directory.

    Indicate which scans need to be processed.

    Returns
    -------
    None
    """
    # Get a list of all directories containing .nii files
    print(f"  * Searching {PATHS['raw']} for all *.nii files")
    raw_niis = fast_recursive_glob_nii(PATHS["raw"])

    # Find the subject ID, scan type, acquisition date, and LONI image ID
    # for each nifti file in raw
    raw_niis = pd.DataFrame(raw_niis, columns=["raw_niif"])
    raw_niis.insert(0, "subj", raw_niis["raw_niif"].apply(get_subj))
    raw_niis.insert(1, "scan_type", raw_niis["raw_niif"].apply(get_scan_type))
    raw_niis.insert(2, "scan_date", raw_niis["raw_niif"].apply(get_scan_date))
    raw_niis.insert(3, "image_id", raw_niis["raw_niif"].apply(get_image_id))

    # Convert the date column to datetime
    raw_niis["scan_date"] = pd.to_datetime(raw_niis["scan_date"])

    # Separate MRI and PET scans into separate dataframes
    raw_mris = raw_niis.query("(scan_type=='MRI-T1')").reset_index(drop=True)
    raw_pets = raw_niis.query("(scan_type!='MRI-T1')").reset_index(drop=True)

    # Rename columns
    raw_mris = raw_mris.rename(
        columns={
            "scan_date": "mri_date",
            "image_id": "mri_image_id",
            "raw_niif": "mri_raw_niif",
        }
    )
    raw_pets = raw_pets.rename(
        columns={
            "scan_type": "tracer",
            "scan_date": "pet_date",
            "image_id": "pet_image_id",
            "raw_niif": "pet_raw_niif",
        }
    )

    # Get rid of scan_type column in the MRI dataframe (all scans are T1)
    raw_mris = raw_mris.drop(columns=["scan_type"])

    # Sort MRIs by subject and scan date
    raw_mris = raw_mris.sort_values(["subj", "mri_date"]).reset_index(drop=True)

    # Sort PET scans by subject, tracer, and scan date
    raw_pets = raw_pets.sort_values(["subj", "tracer", "pet_date"]).reset_index(
        drop=True
    )

    # Get PET resolution from filename
    raw_pets.insert(
        raw_pets.columns.tolist().index("pet_date") + 1,
        "pet_res",
        raw_pets["pet_raw_niif"].apply(get_pet_resolution),
    )

    # Report how many MRI and PET scans are in the raw directory
    print(
        "  * Found {:,} niftis in {:,} subdirectories".format(
            len(raw_niis), len([op.dirname(f) for f in raw_niis["raw_niif"]])
        )
    )
    print(f"    - {len(raw_mris):,} T1 MRIs")
    print(f"    - {len(raw_pets):,} PET scans")
    for scan_type in raw_pets["tracer"].unique():
        print(
            "      * {} {}".format(
                len(raw_pets.loc[raw_pets["tracer"] == scan_type]), scan_type
            )
        )

    # Calculate the number of MRIs for each subject and the time between
    # them
    raw_mris = add_mri_date_columns(raw_mris)

    # Calculate the number of PET scans for each subject and tracer, and
    # the time between them
    raw_pets = add_pet_date_columns(raw_pets)

    # Match each PET scan to its closest MRI
    raw_pets = find_closest_mri_to_pet(raw_pets, raw_mris)
    print("  * Matching PET scans to their closest T1 MRIs")

    # Figure out which MRIs are actually used for PET processing
    raw_mris["mri_used_for_pet_proc"] = raw_mris["mri_image_id"].apply(
        lambda x: 1 if np.isin(x, raw_pets["mri_image_id"]) else 0
    )
    idx = raw_mris.loc[raw_mris["mri_used_for_pet_proc"] == 0].index.tolist()
    if process_unused_mris:
        msg = "      to any PET scan but will be processed in any case"
    else:
        msg = (
            "      to any PET scan and will not be added to the list of MRIs to process"
        )
    print(
        "  * Auditing MRI scans",
        "    - {:,}/{:,} MRIs in {} are not the closest MRI".format(
            len(idx),
            len(raw_mris),
            PATHS["raw"],
        ),
        msg,
        sep="\n",
    )

    # Flag PET scans with issues that preclude processing
    raw_pets = audit_pet(raw_pets)
    print("  * Auditing PET scans")
    print(
        "    - Flagged {:,} scans with issues to be resolved before processing".format(
            len(raw_pets.loc[raw_pets["flag"] == 1])
        )
    )

    # Add path to the processed MRI directory
    ii = raw_mris.columns.tolist().index("mri_raw_niif")
    raw_mris["mri_proc_dir"] = raw_mris.apply(
        lambda x: get_mri_proc_dir(x["subj"], x["mri_date"]), axis=1
    )

    # Add path to the processed PET directory
    ii = raw_pets.columns.tolist().index("pet_raw_niif")
    raw_pets.insert(
        ii + 1,
        "pet_proc_dir",
        raw_pets.apply(
            lambda x: get_pet_proc_dir(x["subj"], x["tracer"], x["pet_date"]), axis=1
        ),
    )

    # Determine which MRIs have been processed
    raw_mris["freesurfer_run"] = raw_mris["mri_proc_dir"].apply(check_if_freesurfer_run)
    print(
        "  * {:,}/{:,} MRIs have been processed through Freesurfer".format(
            len(raw_mris.loc[raw_mris["freesurfer_run"] == 1]),
            len(raw_mris),
        )
    )

    raw_mris["mri_processed"] = raw_mris["mri_proc_dir"].apply(check_if_mri_processed)
    print(
        "  * {:,}/{:,} MRIs have been fully processed".format(
            len(raw_mris.loc[raw_mris["mri_processed"] == 1]),
            len(raw_mris),
        )
    )

    # Determine which PET scans have been processed
    raw_pets.insert(
        ii + 1, "pet_processed", raw_pets["pet_proc_dir"].apply(check_if_pet_processed)
    )
    print(
        "  * {:,}/{:,} PET scans have already been processed".format(
            len(raw_pets.loc[raw_pets["pet_processed"] == 1]),
            len(raw_pets),
        )
    )

    # Determine which MRI scans need to be processed
    raw_mris["need_to_process"] = raw_mris.apply(
        lambda x: get_mri_to_process(
            x["mri_used_for_pet_proc"],
            x["mri_processed"],
            overwrite,
            process_unused_mris,
        ),
        axis=1,
    )
    print(
        "  * {:,}/{:,} MRIs are scheduled for processing".format(
            len(raw_mris.loc[raw_mris["need_to_process"] == 1]),
            len(raw_mris),
        )
    )

    # Determine which PET scans need to be processed
    raw_pets["need_to_process"] = raw_pets.apply(
        lambda x: get_pet_to_process(x["flag"], x["pet_processed"], overwrite), axis=1
    )
    print(
        "  * {:,}/{:,} PET scans are scheduled for processing".format(
            len(raw_pets.loc[raw_pets["need_to_process"] == 1]),
            len(raw_pets),
        )
    )

    # Convert date columns to string
    fmt = "%Y-%m-%d"
    raw_mris["mri_date"] = raw_mris["mri_date"].dt.strftime(fmt)
    for col in ["pet_date", "mri_date"]:
        raw_pets[col] = raw_pets[col].dt.strftime(fmt)

    # Save the raw MRI scans dataframe to a CSV file
    save_mri_scan_index(raw_mris)

    # Save the raw PET scans dataframe to a CSV file
    save_pet_scan_index(raw_pets)


def fast_recursive_glob_nii(path):
    """Return a list of all files in path that end in .nii"""

    def _path_recurse(path):
        """Recursively traverse a directory and its subdirectories"""
        with os.scandir(path) as entries:
            for entry in entries:
                if entry.is_dir():
                    _path_recurse(entry.path)
                elif entry.is_file() and entry.name.endswith(".nii"):
                    nii_files.append(entry.path)

    nii_files = []
    _path_recurse(path)
    return nii_files


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


def get_subj(filepath, raw_dir=None):
    """Return the subject ID from filepath to the recon'd nifti.

    Parameters
    ----------
    filepath : str
        The filepath to the reconstructed nifti.

    Returns
    -------
    subj : str
        The subject ID parsed from the input file basename.
    """
    if raw_dir is None:
        raw_dir = PATHS["raw"]
    try:
        subj = filepath.replace(raw_dir + "/", "").split("/")[0]
        if len(subj) > 0:
            return subj
        else:
            return np.nan
    except IndexError:
        return np.nan


def get_scan_type(filepath):
    """Search filepath and return the scan type"""

    def find_matches(text, scan_types, scan_type_keys):
        """Helper function to find matches"""
        matches = [key for key in scan_type_keys if key in text]
        unique_values = list(set([scan_types[key] for key in matches]))
        return unique_values

    # Load the SCAN_TYPES dict and get lowercase keys
    scan_type_keys = list(SCAN_TYPES)

    # Convert filepath to lowercase and get the basename
    filepath = filepath.lower()
    basename = op.basename(filepath)

    # First attempt to match in the basename
    match_values = find_matches(basename, SCAN_TYPES, scan_type_keys)
    if len(match_values) == 1:
        return match_values[0]
    elif len(match_values) > 1:
        return "FILENAME MATCHED MULTIPLE PET TRACERS OR MRI MODALITIES"

    # Check for specific substring in basename
    if "coreg,_avg,_std_img_and_vox_siz,_uniform" in basename:
        return "FDG"

    # Attempt to match in the directory path after "raw"
    dirname = op.dirname(filepath)
    raw_base = "/{}/".format(op.basename(PATHS["raw"]))
    raw_idx = dirname.find(raw_base)
    if raw_idx == -1:
        return "FAILED TO IDENTIFY PET TRACER OR MRI MODALITY FROM FILENAME"

    # Only consider part of the path after "raw"
    search_path = dirname[raw_idx + len(raw_base) :]
    match_values = find_matches(search_path, SCAN_TYPES, scan_type_keys)
    if len(match_values) == 1:
        return match_values[0]
    elif len(match_values) > 1:
        return "FILENAME MATCHED MULTIPLE PET TRACERS OR MRI MODALITIES"
    else:
        if "coreg,_avg,_std_img_and_vox_siz,_uniform" in search_path:
            return "FDG"
        else:
            return "FAILED TO IDENTIFY PET TRACER OR MRI MODALITY FROM FILENAME"


def get_scan_date(filepath):
    """Search filepath and return the scan date"""

    def test_datestr(str_to_search, date_start="2000-01-01", date_stop=None):
        """Return date if start of string is in YYYYMMDD format

        The date must be >= date_start and <= today's date.
        """
        # Remove hypens and underscores
        str_to_search = str_to_search.replace("-", "").replace("_", "")

        # Check if the first 8 characters are numeric
        if str_to_search is None or not re.match(r"^\d{8}", str_to_search):
            return False

        # Extract the date part of the string
        datestr = "{}-{}-{}".format(
            str_to_search[:4], str_to_search[4:6], str_to_search[6:8]
        )

        # Check if the date is valid and within the specified range
        date_fmt = "%Y-%m-%d"
        try:
            date_test = datetime.datetime.strptime(datestr, date_fmt).date()
            date_start = datetime.datetime.strptime(date_start, date_fmt).date()
            if date_stop is None:
                date_stop = datetime.date.today()
            else:
                date_stop = datetime.datetime.strptime(date_stop, date_fmt).date()
            if date_start <= date_test <= date_stop:
                return datestr
            else:
                return False
        except ValueError:
            return False

    basename = op.basename(filepath)
    parts = basename.split("_")
    for part in parts:
        datestr = test_datestr(part)
        if datestr:
            return datestr

    # If no date was found in the basename, try the dirname but limit
    # search to everything forward from the raw/ directory
    dirname = op.dirname(filepath)
    raw_base = "/{}/".format(op.basename(PATHS["raw"]))
    raw_idx = dirname.find(raw_base)
    if raw_idx == -1:
        return np.nan
    search_path = dirname[raw_idx + len(raw_base) :]
    parts = search_path.split("/")
    for part in parts:
        datestr = test_datestr(part)
        if datestr:
            return datestr

    return np.nan


def datetime_to_datestr(dt):
    """Convert a datetime object to a date string."""
    try:
        return dt.strftime("%Y-%m-%d")
    except ValueError:
        return np.nan


def get_image_id(filepath):
    """Search basename and return the LONI Image ID"""

    def test_image_id(str_to_search, pattern=r"^I[0-9]+$"):
        """Return str_to_search if it starts with 'I' and continues with numbers"""
        if str_to_search is None:
            return False
        if bool(re.match(pattern, str_to_search)):
            return str_to_search
        else:
            return False

    basename = op.basename(filepath)
    parts = basename.split("_")
    for part in parts:
        image_id = test_image_id(part)
        if image_id:
            return image_id
    return np.nan


def get_pet_resolution(filepath):
    """Search basename and return the PET resolution"""

    def test_resolution(str_to_search, pattern=r"uniform_([0-9])mm_res"):
        """Return scan resolution if the expected pattern is found"""
        if str_to_search is None:
            return False
        match = re.search(pattern, str_to_search)
        if match:
            return int(match.group(1))
        return False

    res = test_resolution(op.basename(filepath).lower())
    if res:
        return res
    else:
        return np.nan


def add_mri_date_columns(mri_scans):
    """Add info on time between PET and MRI and adjacent PET scans"""
    # Copy the input dataframe
    mri_scans_cp = mri_scans.copy()

    # Add columns for PET scan number and total number of scans per tracer
    ii = mri_scans_cp.columns.tolist().index("mri_date")
    grp = mri_scans_cp.groupby("subj")
    mri_scans_cp.insert(ii + 1, "mri_scan_number", grp.cumcount() + 1)
    mri_scans_cp.insert(ii + 2, "n_mri_scans", grp["mri_date"].transform("count"))

    # Add columns for days from each PET scan to baseline and days between
    # consecutive PET scans per tracer
    baseline_mri_dates = grp["mri_date"].min()
    mri_scans_cp.insert(
        ii + 3,
        "days_from_baseline_mri",
        mri_scans_cp.apply(
            lambda x: (x["mri_date"] - baseline_mri_dates[x["subj"]]).days,
            axis=1,
        ),
    )
    mri_scans_cp.insert(
        ii + 4, "days_from_last_mri", grp["mri_date"].diff().dt.days.fillna(0)
    )

    return mri_scans_cp


def add_pet_date_columns(pet_scans):
    """Add info on time between PET and MRI and adjacent PET scans"""
    # Copy the input dataframe
    pet_scans_cp = pet_scans.copy()

    # Add columns for PET scan number and total number of scans per tracer
    ii = pet_scans_cp.columns.tolist().index("pet_date")
    grp = pet_scans_cp.groupby(["subj", "tracer"])
    pet_scans_cp.insert(ii + 1, "pet_scan_number", grp.cumcount() + 1)
    pet_scans_cp.insert(ii + 2, "n_pet_scans", grp["pet_date"].transform("count"))

    # Add columns for days from each PET scan to baseline and days between
    # consecutive PET scans per tracer
    baseline_pet_dates = grp["pet_date"].min()
    pet_scans_cp.insert(
        ii + 3,
        "days_from_baseline_pet",
        pet_scans_cp.apply(
            lambda x: (
                x["pet_date"] - baseline_pet_dates[(x["subj"], x["tracer"])]
            ).days,
            axis=1,
        ),
    )
    pet_scans_cp.insert(
        ii + 4, "days_from_last_pet", grp["pet_date"].diff().dt.days.fillna(0)
    )

    return pet_scans_cp


def find_closest_mri_to_pet(pet_scans, mri_scans):
    """Return a merged dataframe with the closest MRI to each PET scan"""
    # Check that the necessary columns are present
    assert np.all(
        np.isin(["subj", "tracer", "pet_date", "pet_image_id"], pet_scans.columns)
    )
    assert np.all(
        np.isin(
            ["subj", "mri_date", "mri_scan_number", "n_mri_scans"], mri_scans.columns
        )
    )

    # Copy the input dataframes
    pet_scans_cp = pet_scans.copy()
    mri_scans_cp = mri_scans.copy()

    # Merge the dataframes
    merged = pet_scans_cp.merge(mri_scans_cp, on="subj", how="left")

    # Compute the date difference between each PET and MRI scan
    ii = merged.columns.tolist().index("pet_date")
    merged.insert(
        ii + 1,
        "pet_to_mri_days",
        (merged["pet_date"] - merged["mri_date"]).abs().dt.days,
    )

    # Sort by |pet_date - mri_date| and drop duplicates
    merged_min = merged.sort_values("pet_to_mri_days").drop_duplicates("pet_image_id")

    # Resort the dataframe and reset index before returning
    merged_min = merged_min.sort_values(["subj", "tracer", "pet_date"]).reset_index(
        drop=True
    )

    return merged_min


def audit_pet(pet_scans):
    """Audit each PET scan and flag scans with potential issues"""
    pet_scans_cp = pet_scans.copy()
    pet_scans_cp["flag"] = 0
    pet_scans_cp["notes"] = ""

    # Missing tracer
    idx = pet_scans_cp.loc[
        pet_scans_cp["tracer"].str.startswith("FAILED TO IDENTIFY")
    ].index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[
        idx, "notes"
    ] += "Unknown tracer, check raw PET filename against scan_types_and_tracers.csv\n"

    # Multiple tracers
    idx = pet_scans_cp.loc[
        pet_scans_cp["tracer"].str.startswith("FILENAME MATCHED MULTIPLE")
    ].index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[
        idx, "notes"
    ] += (
        "Matched >1 tracer, check raw PET filename against scan_types_and_tracers.csv\n"
    )

    # Missing PET date
    idx = pet_scans_cp.loc[pd.isna(pet_scans_cp["pet_date"])].index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[idx, "notes"] += "Missing PET date\n"

    # PET resolution unclear or not at 6mm
    idx = pet_scans_cp.loc[pet_scans_cp["pet_res"] != 6].index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[idx, "notes"] += "PET resolution is unclear or not at 6mm\n"

    # Missing MRI
    idx = pet_scans_cp.loc[pd.isna(pet_scans_cp).any(axis=1)].index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[idx, "notes"] += f"No available MRI in {PATHS['raw']}\n"

    # PET and MRI scan dates > 365 days apart
    idx = pet_scans_cp.query("(pet_to_mri_days>365)").index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[idx, "notes"] += "Closest MRI is more than 1 year from PET date\n"

    # Same MRI used to process multiple PET scans for the same tracer
    idx = pet_scans_cp.loc[
        pet_scans_cp.groupby(["subj", "tracer"])["mri_image_id"].transform(
            lambda x: (
                True if (np.max(np.unique(x, return_counts=True)[1]) > 1) else False
            )
        )
    ].index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[
        idx, "notes"
    ] += "Same MRI would be used to process multiple timepoints for the same PET tracer\n"

    return pet_scans_cp


def get_mri_proc_dir(subj, mri_date):
    """Return the processed MRI directory for each MRI scan"""
    return op.join(
        PATHS["processed"], subj, "MRI-T1_{}".format(datetime_to_datestr(mri_date))
    )


def get_pet_proc_dir(subj, tracer, pet_date):
    """Return the processed PET directory for each PET scan"""
    return op.join(
        PATHS["processed"],
        subj,
        "{}_{}".format(tracer, datetime_to_datestr(pet_date)),
    )


def check_if_freesurfer_run(mri_proc_dir):
    """Return True if the MRI scan has been processed by Freesurfer"""
    freesurfer_dir = op.join(mri_proc_dir, "freesurfer")
    if not op.isdir(freesurfer_dir):
        return 0
    # Search for the Freesurfer output directory
    nuf = op.join(freesurfer_dir, "mri", "nu.mgz")
    aparcf = op.join(freesurfer_dir, "mri", "aparc+aseg.mgz")
    bstemf = op.join(freesurfer_dir, "mri", "brainstemSsLabels.v12.FSvoxelSpace.mgz")
    if op.isfile(nuf) and op.isfile(aparcf) and op.isfile(bstemf):
        return 1
    else:
        return 0


def check_if_mri_processed(mri_proc_dir):
    """Return True if the PET scan has been processed"""
    if not op.isdir(mri_proc_dir):
        return 0
    # Search for the warped and affine transformed nu.nii,
    # respectively, as these are among the last files created in the MRI
    # processing pipeline. Also make sure we find at least one mask.
    mfiles = glob(op.join(mri_proc_dir, "*mask-*.nii"))
    wfiles = glob(op.join(mri_proc_dir, "w*_nu.nii"))
    afiles = glob(op.join(mri_proc_dir, "a*_nu.nii"))
    if mfiles and wfiles and afiles:
        return 1
    else:
        return 0


def check_if_pet_processed(pet_proc_dir):
    """Return True if the PET scan has been processed"""
    if not op.isdir(pet_proc_dir):
        return 0
    # Search for the warped and affine transformed PET SUVRs,
    # respectively, as these are among the last files created in the PET
    # processing pipeline
    wfiles = glob(op.join(pet_proc_dir, "w*suvr*.nii"))
    afiles = glob(op.join(pet_proc_dir, "a*suvr*.nii"))
    if wfiles and afiles:
        return 1
    else:
        return 0


def get_mri_to_process(used_for_pet, processed, overwrite, process_unused_mris):
    """Return 1 if the MRI scan should be processed, otherwise 0"""
    if process_unused_mris or used_for_pet:
        if overwrite or not processed:
            return 1
    return 0


def get_pet_to_process(flagged, processed, overwrite):
    """Return 1 if the PET scan should be processed, otherwise 0"""
    if not flagged:
        if overwrite or not processed:
            return 1
    return 0


def save_mri_scan_index(mri_scans):
    """Save the MRI scan index to a CSV file"""
    # Move any existing files to archive
    archive_dir = op.join(PATHS["scans_to_process"], "archive")
    files = glob(op.join(PATHS["scans_to_process"], "Raw_MRI_Scan_Index_*.csv"))
    if files:
        if not op.isdir(archive_dir):
            os.makedirs(archive_dir)
        for f in files:
            os.rename(f, op.join(archive_dir, op.basename(f)))

    # Save the mri_scans dataframe
    outf = op.join(PATHS["scans_to_process"], f"Raw_MRI_Scan_Index_{TIMESTAMP}.csv")
    mri_scans.to_csv(outf, index=False)
    print(f"  * Saved raw MRI scan index to {outf}")


def save_pet_scan_index(pet_scans):
    """Save the PET scan index to a CSV file"""
    # Move any existing files to archive
    archive_dir = op.join(PATHS["scans_to_process"], "archive")
    files = glob(op.join(PATHS["scans_to_process"], "Raw_PET_Scan_Index_*.csv"))
    if files:
        if not op.isdir(archive_dir):
            os.makedirs(archive_dir)
        for f in files:
            os.rename(f, op.join(archive_dir, op.basename(f)))

    # Save the pet_scans dataframe
    outf = op.join(PATHS["scans_to_process"], f"Raw_PET_Scan_Index_{TIMESTAMP}.csv")
    pet_scans.to_csv(outf, index=False)
    print(f"  * Saved raw PET scan index to {outf}")


def _parse_args():
    """Parse and return command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Find MRI and PET scans that need to be processed by comparing data\n"
            + "in raw vs. processed directories"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        exit_on_error=False,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help=(
            "Process all PET scans and their closest MRIs in the raw directory,\n"
            + "whether or not they have already been processed"
        ),
    )
    parser.add_argument(
        "-p",
        "--process_unused_mris",
        action="store_true",
        help=(
            "Process all MRIs scans in the raw directory, even if they are not the\n"
            + "closest MRI to any PET scan currently in the raw directory"
        ),
    )

    # Parse the command line arguments
    return parser.parse_args()


if __name__ == "__main__":
    # Get command line arguments.
    args = _parse_args()

    # Call the main function
    main(overwrite=args.overwrite, process_unused_mris=args.process_unused_mris)

    # Exit successfully
    sys.exit(0)
