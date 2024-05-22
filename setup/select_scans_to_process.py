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


# Define global variables
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

    def is_dst():
        today = datetime.datetime.now()
        year = today.year
        dst_start = datetime.datetime(year, 3, 8, 2, 0)  # Second Sunday in March
        dst_end = datetime.datetime(year, 11, 1, 2, 0)  # First Sunday in November
        dst_start += datetime.timedelta(
            days=(6 - dst_start.weekday())
        )  # Find the second Sunday
        dst_end += datetime.timedelta(
            days=(6 - dst_end.weekday())
        )  # Find the first Sunday
        return dst_start <= today < dst_end

    # Get the current UTC time
    utc_time = datetime.datetime.now(datetime.timezone.utc)

    # Define the offset for Los Angeles time
    la_offset = -8 if not is_dst() else -7

    # Create a timezone-aware datetime object for Los Angeles
    la_time = utc_time + datetime.timedelta(hours=la_offset)

    # Format the time as specified
    formatted_time = la_time.strftime("%Y-%m-%d-%H-%M-%S")

    return formatted_time


set_globals()


def main(
    process_all_mris=False, overwrite=False, max_days_from_present=100, save_csv=True
):
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
    print("  * Matching each PET scan to its closest MRI")
    raw_pets = find_closest_mri_to_pet(raw_pets, raw_mris)

    # Figure out which MRIs are actually used for PET processing
    raw_mris = get_orphan_mris(raw_mris, raw_pets, max_days_from_present)
    print(
        f"  * Auditing MRIs in {PATHS['raw']}",
        "    - {:,} MRIs are orphans (>{} days old and not linked to any PET scan).".format(
            raw_mris["mri_is_orphan"].sum(), max_days_from_present
        ),
        sep="\n",
    )
    if process_all_mris:
        msg = "      These scans will still be scheduled for processing at the user's request"
    else:
        msg = "      These scans will not be scheduled for processing (process anyway with flag --all)"
    print(msg)

    # Flag PET scans with issues that preclude processing
    raw_pets = audit_pet(raw_pets)
    print(f"  * Auditing PET scans in {PATHS['raw']}")
    print(
        "    - Flagged {:,} scans with issues to be resolved before processing (see 'flag_notes')".format(
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
    raw_mris["freesurfer_complete"] = raw_mris["mri_proc_dir"].apply(
        check_if_freesurfer_run
    )
    raw_mris["mri_processing_complete"] = raw_mris["mri_proc_dir"].apply(
        check_if_mri_processed
    )

    # Schedule MRIs for processing
    baseline_mri_processed = (
        raw_mris.loc[raw_mris["mri_scan_number"] == 1]
        .set_index("subj")["mri_processing_complete"]
        .to_dict()
    )
    raw_mris["scheduled_for_processing"] = raw_mris.apply(
        lambda x: get_mri_to_process(
            x["mri_scan_number"],
            x["mri_is_orphan"],
            x["mri_processing_complete"],
            baseline_mri_processed[x["subj"]],
            overwrite,
            process_all_mris,
        ),
        axis=1,
    )

    # Summarize how many MRIs have been or still need to be processed
    raw_mris_req = raw_mris.query("(mri_scan_number==1) or (mri_is_orphan==0)")
    raw_mris_opt = raw_mris.query("(mri_scan_number>1) and (mri_is_orphan==1)")
    print("\n  * MRI scan review:")
    print(
        "    - {:,}/{:,} required MRIs have been processed through Freesurfer".format(
            len(raw_mris_req.loc[raw_mris_req["freesurfer_complete"] == 1]),
            len(raw_mris_req),
        )
    )
    print(
        "    - {:,}/{:,} required MRIs have been fully processed".format(
            len(raw_mris_req.loc[raw_mris_req["mri_processing_complete"] == 1]),
            len(raw_mris_req),
        )
    )
    print(
        "    - {:,}/{:,} optional MRIs have been processed through Freesurfer".format(
            len(raw_mris_opt.loc[raw_mris_opt["freesurfer_complete"] == 1]),
            len(raw_mris_opt),
        )
    )
    print(
        "    - {:,}/{:,} optional MRIs have been fully processed".format(
            len(raw_mris_opt.loc[raw_mris_opt["mri_processing_complete"] == 1]),
            len(raw_mris_opt),
        )
    )
    _n = len(
        raw_mris_req.query(
            "(mri_processing_complete == 0) & (scheduled_for_processing == 0)"
        )
    )
    if _n:
        print(
            f"    - {_n:,} required, follow-up MRIs cannot be processed until baseline MRI is processed"
        )

    # Determine which PET scans have been processed
    ii = raw_pets.columns.tolist().index("mri_raw_niif")
    mri_processed = raw_mris.set_index("mri_image_id")[
        "mri_processing_complete"
    ].to_dict()
    raw_pets.insert(
        ii + 1,
        "mri_processing_complete",
        raw_pets["mri_image_id"].apply(lambda x: mri_processed.get(x, np.nan)),
    )
    raw_pets.insert(
        ii + 2,
        "pet_processing_complete",
        raw_pets["pet_proc_dir"].apply(check_if_pet_processed),
    )

    # Schedule PET scans for processing
    raw_pets["scheduled_for_processing"] = raw_pets.apply(
        lambda x: get_pet_to_process(
            x["pet_processing_complete"],
            x["mri_processing_complete"],
            x["flag"],
            overwrite,
        ),
        axis=1,
    )

    # Summarize how many PET scans have been or still need to be processed
    print("\n  * PET scan review:")
    print(
        "    - {:,}/{:,} PET scans have been fully processed".format(
            len(raw_pets.loc[raw_pets["pet_processing_complete"] == 1]),
            len(raw_pets),
        )
    )
    qry = (
        "(pet_processing_complete == 0) & (flag == 0) & (scheduled_for_processing == 0)"
    )
    _n = len(raw_pets.query(qry))
    if _n:
        print(
            f"    - {_n:,} PET scans cannot be processed until their corresponding MRI is processed"
        )

    # Save the raw scan dataframes to CSV files
    if save_csv:
        save_raw_mri_index(raw_mris)
        save_raw_pet_index(raw_pets)

    # Report how many MRI and PET scans are scheduled for processing
    print("")
    print("  +.." + ("=" * 44) + "..+")
    print("  |" + (" " * 47) + " |")
    print(
        "  |  {:>5,} MRIs are scheduled for processing       |".format(
            len(raw_mris.loc[raw_mris["scheduled_for_processing"] == 1])
        )
    )
    print(
        "  |  {:>5,} PET scans are scheduled for processing  |".format(
            len(raw_pets.loc[raw_pets["scheduled_for_processing"] == 1])
        )
    )
    print("  |" + (" " * 47) + " |")
    print("  +.." + ("=" * 44) + "..+")
    print("")

    return raw_mris, raw_pets


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
        "days_mri_to_pet",
        (merged["pet_date"] - merged["mri_date"]).dt.days,
    )
    merged.insert(
        ii + 2,
        "abs_days_mri_to_pet",
        merged["days_mri_to_pet"].abs(),
    )

    # Sort by |pet_date - mri_date| and drop duplicates
    merged_min = merged.sort_values("abs_days_mri_to_pet").drop_duplicates(
        ["subj", "tracer", "pet_date"]
    )

    # Resort the dataframe and reset index before returning
    merged_min = merged_min.sort_values(["subj", "tracer", "pet_date"]).reset_index(
        drop=True
    )

    return merged_min


def get_orphan_mris(mri_scans, pet_scans, max_days_from_present=100):
    """Flag MRIs that are not used for PET processing"""
    ii = mri_scans.columns.tolist().index("n_mri_scans")
    mri_scans.insert(
        ii + 1,
        "mri_used_for_pet_proc",
        mri_scans["mri_image_id"].apply(
            lambda x: 1 if np.isin(x, pet_scans["mri_image_id"]) else 0
        ),
    )
    today = pd.Timestamp("today")
    mri_scans.insert(
        ii + 2,
        "mri_is_orphan",
        mri_scans.apply(
            lambda x: (
                1
                if (x["mri_used_for_pet_proc"] == 0)
                and ((today - x["mri_date"]).days > max_days_from_present)
                else 0
            ),
            axis=1,
        ),
    )
    return mri_scans


def audit_pet(
    pet_scans,
    audit_pet_res=True,
    audit_pet_to_mri_days=True,
    audit_repeat_mri=True,
    presumed_pet_res=6,
    max_pet_to_mri_days=365,
):
    """Audit each PET scan and flag scans with potential issues

    Parameters
    ----------
    pet_scans : pd.DataFrame
        A DataFrame with columns 'subj', 'tracer', 'pet_date', 'pet_res',
        'mri_image_id', 'flag', and 'flag_notes'
    audit_pet_res : bool, optional
        If True, audit the PET resolution
    audit_pet_to_mri_days : bool, optional
        If True, audit the days between PET and MRI scans
    audit_repeat_mri : bool, optional
        If True, audit for repeated MRIs used for multiple PET scans
    presumed_pet_res : int, optional
        The presumed PET resolution in mm
    max_pet_to_mri_days : int, optional
        The maximum allowed number of days between PET and MRI scans

    Returns
    -------
    pet_scans_cp : pd.DataFrame
        A copy of the input DataFrame with added 'flag' and 'flag_notes'
        columns. Scans with flag==1 have potential issues that need to
        be resolved before processing. These issues are noted in the
        'flag_notes' column. Otherwise flag==0 and flag_notes==""
    """
    pet_scans_cp = pet_scans.copy()
    pet_scans_cp["flag"] = 0
    pet_scans_cp["flag_notes"] = ""

    # Missing tracer
    idx = pet_scans_cp.loc[
        pet_scans_cp["tracer"].str.startswith("FAILED TO IDENTIFY")
    ].index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[
        idx, "flag_notes"
    ] += "Unknown tracer--check raw PET filename against scan_types_and_tracers.csv; "

    # Multiple tracers
    idx = pet_scans_cp.loc[
        pet_scans_cp["tracer"].str.startswith("FILENAME MATCHED MULTIPLE")
    ].index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[
        idx, "flag_notes"
    ] += (
        "Matched >1 tracer--check raw PET filename against scan_types_and_tracers.csv; "
    )

    # Missing PET date
    idx = pet_scans_cp.loc[pd.isna(pet_scans_cp["pet_date"])].index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[idx, "flag_notes"] += "Missing PET date; "

    # PET resolution unclear or not at 6mm
    if audit_pet_res:
        idx = pet_scans_cp.loc[
            pet_scans_cp["pet_res"] != presumed_pet_res
        ].index.tolist()
        pet_scans_cp.loc[idx, "flag"] = 1
        pet_scans_cp.loc[
            idx, "flag_notes"
        ] += f"PET resolution could not be parsed or is not at {presumed_pet_res}mm; "

    # Missing MRI
    idx = pet_scans_cp.loc[pd.isna(pet_scans_cp).any(axis=1)].index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[idx, "flag_notes"] += f"No available MRI in {PATHS['raw']}; "

    # PET and MRI scan dates too far apart
    if audit_pet_to_mri_days:
        idx = pet_scans_cp.query(
            f"abs_days_mri_to_pet > {max_pet_to_mri_days}"
        ).index.tolist()
        pet_scans_cp.loc[idx, "flag"] = 1
        pet_scans_cp.loc[
            idx, "flag_notes"
        ] += f"Closest MRI is more than {max_pet_to_mri_days} days from PET date; "

    # Same MRI used to process multiple PET scans for the same tracer
    if audit_repeat_mri:
        idx = pet_scans_cp.loc[
            pet_scans_cp.groupby(["subj", "tracer"])["mri_image_id"].transform(
                lambda x: (
                    True if (np.max(np.unique(x, return_counts=True)[1]) > 1) else False
                )
            )
        ].index.tolist()
        pet_scans_cp.loc[idx, "flag"] = 1
        pet_scans_cp.loc[
            idx, "flag_notes"
        ] += "Same MRI would be used to process multiple timepoints for the same PET tracer; "

    pet_scans_cp["flag_notes"] = pet_scans_cp["flag_notes"].apply(
        lambda x: x[:-2] if (len(x) >= 2) else x
    )
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


def get_mri_to_process(
    mri_scan_number,
    mri_is_orphan,
    already_processed,
    baseline_processed,
    overwrite,
    process_all_mris,
):
    """Return 1 if the MRI scan should be processed, otherwise 0"""
    # Don't process MRIs that have already been processed unless
    # overwrite is True
    if already_processed and not overwrite:
        return 0

    # Always process baseline MRIs as later MRIs are coreg'd to them
    if mri_scan_number == 1:
        return 1

    # Only process follow-up MRIs if the baseline MRI has been processed
    if baseline_processed:
        # Don't process orphan MRIs unless process_all_mris is True
        if process_all_mris or not mri_is_orphan:
            return 1
        else:
            return 0
    else:
        return 0


def get_pet_to_process(pet_processed, mri_processed, pet_flagged, overwrite):
    """Return 1 if the PET scan should be processed, otherwise 0"""
    # Don't process PET scans that have already been processed unless
    # overwrite is True
    if pet_processed and not overwrite:
        return 0

    # Don't process PET scans that have been flagged with issues
    if pet_flagged:
        return 0

    # Only process PET scans for which the corresponding MRI has already
    # been processed
    if mri_processed:
        return 1
    else:
        return 0


def save_raw_mri_index(mri_scans):
    """Save the MRI scan index to a CSV file"""
    # Move any existing files to archive
    archive_dir = op.join(PATHS["scans_to_process"], "archive")
    files = glob(op.join(PATHS["scans_to_process"], "raw_MRI_index_*.csv"))
    if files:
        if not op.isdir(archive_dir):
            os.makedirs(archive_dir)
        for f in files:
            os.rename(f, op.join(archive_dir, op.basename(f)))

    # Resort the dataframe so we process scans that have already been
    # run through FreeSurfer first
    mri_scans = mri_scans.sort_values(
        [
            "freesurfer_complete",
            "subj",
            "mri_date",
        ],
        ascending=[False, True, True],
    ).reset_index(drop=True)

    # Save the mri_scans dataframe
    outf = op.join(PATHS["scans_to_process"], f"raw_MRI_index_{TIMESTAMP}.csv")

    mri_scans.to_csv(outf, index=False)
    print(f"  * Saved raw MRI scan index to {outf}")


def save_raw_pet_index(pet_scans):
    """Save the PET scan index to a CSV file"""
    # Move any existing files to archive
    archive_dir = op.join(PATHS["scans_to_process"], "archive")
    files = glob(op.join(PATHS["scans_to_process"], "raw_PET_index_*.csv"))
    if files:
        if not op.isdir(archive_dir):
            os.makedirs(archive_dir)
        for f in files:
            os.rename(f, op.join(archive_dir, op.basename(f)))

    # Save the pet_scans dataframe
    outf = op.join(PATHS["scans_to_process"], f"raw_PET_index_{TIMESTAMP}.csv")
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
        "-a",
        "--all",
        action="store_true",
        help=(
            "Process all MRIs scans in the raw directory, even if they are not the\n"
            + "closest MRI to any PET scan currently in the raw directory"
        ),
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
        "-m",
        "--max_days",
        type=int,
        default=100,
        help=(
            "The maximum number of days from present that an MRI scan can have been acquired\n"
            + "without being linked to a PET scan before it is considered an orphan and is not\n"
            + "scheduled for processing when --all is not used"
        ),
    )

    # Parse the command line arguments
    return parser.parse_args()


if __name__ == "__main__":
    # Get command line arguments.
    args = _parse_args()

    # Call the main function
    main(
        process_all_mris=args.all,
        overwrite=args.overwrite,
        max_days_from_present=args.max_days,
    )

    # Exit successfully
    sys.exit(0)
