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

utils_dir = op.join(op.dirname(__file__), "..", "utils")
if utils_dir not in sys.path:
    sys.path.append(utils_dir)
import utilities as uts

# Define globals
AMYLOID_TRACERS = ["FBB", "FBP", "FLUTE", "NAV", "PIB"]


def main(
    proj_dir,
    orphan_mri_tolerance=182,
    schedule_all_mris=False,
    overwrite=False,
    audit_pet_res=True,
    expected_pet_res=6,
    audit_pet_to_mri_days=True,
    max_pet_to_mri=365,
    audit_repeat_mri=False,
    check_brainstem=True,
    save_csv=True,
):
    """Save CSV files for MRI and PET scans in the raw directory.

    Indicate which scans need to be processed.

    Returns
    -------
    None
    """
    # Define paths to the needed directories
    raw_dir = op.abspath(op.join(proj_dir, "data", "raw"))
    raw_base = op.basename(raw_dir)
    proc_dir = op.abspath(op.join(proj_dir, "data", "processed"))
    config_dir = op.abspath(op.join(proj_dir, "code", "config"))
    scans_to_process_dir = op.abspath(op.join(proj_dir, "metadata", "scans_to_process"))
    for d in [raw_dir, proc_dir, config_dir, scans_to_process_dir]:
        if not op.isdir(d):
            raise ValueError(f"{d} does not exist")

    # Load the scan types CSV files
    scan_typesf = op.join(config_dir, "scan_types_and_tracers.csv")
    scan_types = load_scan_typesf(scan_typesf)

    # Get a list of all directories containing .nii files
    print(f"  * Searching {raw_dir} for all *.nii files")
    raw_niis = fast_recursive_glob_nii(raw_dir)

    # Find the subject ID, scan type, acquisition date, and LONI image ID
    # for each nifti file in raw
    raw_niis = pd.DataFrame(raw_niis, columns=["raw_niif"])
    raw_niis.insert(
        0, "subj", raw_niis["raw_niif"].apply(lambda x: get_subj(x, raw_dir))
    )
    raw_niis.insert(
        1,
        "scan_type",
        raw_niis["raw_niif"].apply(lambda x: get_scan_type(x, scan_types, raw_base)),
    )
    raw_niis.insert(
        2, "scan_date", raw_niis["raw_niif"].apply(lambda x: get_scan_date(x, raw_base))
    )
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
    raw_mris = get_orphan_mris(raw_mris, raw_pets, orphan_mri_tolerance)
    n_orphans = raw_mris["mri_is_orphan"].sum()
    print("  * Auditing MRIs in {raw_dir}")
    if n_orphans == 0:
        print("    - All MRIs are linked to a PET scan.")
    else:
        print(
            f"    - {n_orphans:,} MRIs are orphans (>{orphan_mri_tolerance} days old and not linked to any PET scan)."
        )
    if schedule_all_mris:
        msg = "      These scans will still be scheduled for processing at the user's request"
    else:
        msg = "      These scans will not be scheduled for processing (process anyway with flag --all)"
    print(msg)

    # Flag PET scans with issues that preclude processing
    raw_pets = audit_pet(
        raw_pets,
        audit_pet_res=audit_pet_res,
        expected_pet_res=expected_pet_res,
        audit_pet_to_mri_days=audit_pet_to_mri_days,
        max_pet_to_mri=max_pet_to_mri,
        audit_repeat_mri=audit_repeat_mri,
    )
    print(f"  * Auditing PET scans in {raw_dir}")
    print(
        "    - Flagged {:,} scans with issues to be resolved before processing (see 'flag_notes')".format(
            len(raw_pets.loc[raw_pets["flag"] == 1])
        )
    )

    # Add path to the processed MRI directory
    ii = raw_mris.columns.tolist().index("mri_raw_niif")
    raw_mris["mri_proc_dir"] = raw_mris.apply(
        lambda x: get_mri_proc_dir(x["subj"], x["mri_date"], proc_dir), axis=1
    )

    # Add path to the processed PET directory
    ii = raw_pets.columns.tolist().index("pet_raw_niif")
    raw_pets.insert(
        ii + 1,
        "pet_proc_dir",
        raw_pets.apply(
            lambda x: get_pet_proc_dir(x["subj"], x["tracer"], x["pet_date"], proc_dir),
            axis=1,
        ),
    )

    # Determine which MRIs have been processed
    raw_mris["freesurfer_complete"] = raw_mris["mri_proc_dir"].apply(
        lambda x: check_if_freesurfer_run(x, check_brainstem)
    )
    raw_mris["mri_seg_complete"] = raw_mris["mri_proc_dir"].apply(
        check_if_mri_seg_complete
    )
    raw_mris["mri_processing_complete"] = raw_mris["mri_proc_dir"].apply(
        check_if_mri_processed
    )

    # Schedule MRIs for processing
    baseline_mri_seg_complete = (
        raw_mris.loc[raw_mris["mri_scan_number"] == 1]
        .set_index("subj")["mri_seg_complete"]
        .to_dict()
    )
    raw_mris["scheduled_for_processing"] = raw_mris.apply(
        lambda x: get_mri_to_process(
            x["mri_scan_number"],
            x["mri_is_orphan"],
            x["mri_processing_complete"],
            baseline_mri_seg_complete[x["subj"]],
            overwrite,
            schedule_all_mris,
        ),
        axis=1,
    )

    # Summarize how many MRIs have been or still need to be processed
    raw_mris_req = raw_mris.query("(mri_scan_number==1) or (mri_is_orphan==0)")
    raw_mris_opt = raw_mris.query("(mri_scan_number>1) and (mri_is_orphan==1)")
    print("\n  MRI scan review\n  ---------------")
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
        "    - {:,}/{:,} orphan MRIs have been processed through Freesurfer".format(
            len(raw_mris_opt.loc[raw_mris_opt["freesurfer_complete"] == 1]),
            len(raw_mris_opt),
        )
    )
    print(
        "    - {:,}/{:,} orphan MRIs have been fully processed".format(
            len(raw_mris_opt.loc[raw_mris_opt["mri_processing_complete"] == 1]),
            len(raw_mris_opt),
        )
    )
    qry = "(mri_processing_complete == 0) & (scheduled_for_processing == 0)"
    _n = len(raw_mris_req.query(qry))
    if _n:
        print(
            f"    - {_n:,} required, follow-up MRIs cannot be processed until baseline MRI is processed:"
        )
        uts.print_list(
            raw_mris_req.query(qry)
            .apply(
                lambda x: f"{x['subj']}_MRI-T1_{datetime_to_datestr(x['mri_date'])}",
                axis=1,
            )
            .tolist()
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
    print("\n  PET scan review\n  ---------------")
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
            f"    - {_n:,} PET scans cannot be processed until their corresponding MRI is processed:"
        )
        uts.print_list(
            raw_pets.query(qry)["pet_proc_dir"].apply(uts.get_scan_tag).tolist()
        )

    # Save the raw scan dataframes to CSV files
    if save_csv:
        timestamp = uts.now()
        save_raw_mri_index(raw_mris, scans_to_process_dir, timestamp)
        save_raw_pet_index(raw_pets, scans_to_process_dir, timestamp)
        print("")

    # Report how many MRI and PET scans are scheduled for processing
    n_mri_scheduled = len(raw_mris.loc[raw_mris["scheduled_for_processing"] == 1])
    n_mri_scheduled_full = len(
        raw_mris.loc[
            (raw_mris["scheduled_for_processing"] == 1)
            & (raw_mris["freesurfer_complete"] == 0)
        ]
    )
    n_mri_scheduled_partial = len(
        raw_mris.loc[
            (raw_mris["scheduled_for_processing"] == 1)
            & (raw_mris["freesurfer_complete"] == 1)
        ]
    )
    n_pet_scheduled = len(raw_pets.loc[raw_pets["scheduled_for_processing"] == 1])
    if n_mri_scheduled == 0 and n_pet_scheduled == 0:
        sp = 1
    else:
        sp = int(np.log10(max(n_mri_scheduled, n_pet_scheduled))) + 1
    sp += int(sp / 3)  # Space to add a comma every 3 digits
    pad = 42 + sp
    if not save_csv:
        pad += 25 + sp
    print("")
    print("+.." + ("=" * pad) + "..+")
    print("|" + (" " * (pad + 3)) + " |")
    msg = {
        "mri": f"|  {n_mri_scheduled:>{sp},} MRI scans are scheduled for processing",
        "pet": f"|  {n_pet_scheduled:>{sp},} PET scans are scheduled for processing",
    }
    if not save_csv:
        for scan_type in msg:
            msg[scan_type] = msg[scan_type].replace(
                "are scheduled for processing",
                "can be scheduled for processing (rerun with save_csv=True)",
            )
    print(msg["mri"])
    print(msg["pet"])
    print("|" + (" " * (pad + 3)) + " |")
    print("+.." + ("=" * pad) + "..+")
    print("")

    # Report the individual scans to be processed
    if n_mri_scheduled > 0:
        print("  All MRIs scheduled to be processed:")
        if n_mri_scheduled_full > 0:
            mris_to_process_full = (
                raw_mris.query(
                    "(scheduled_for_processing==1) & (freesurfer_complete==0)"
                )["mri_proc_dir"]
                .apply(uts.get_scan_tag)
                .tolist()
            )
            print("\n  FreeSurfer + SPM processing", "  " + ("-" * 27), sep="\n")
            uts.print_list(mris_to_process_full)
        if n_mri_scheduled_partial > 0:
            mris_to_process_partial = (
                raw_mris.query(
                    "(scheduled_for_processing==1) & (freesurfer_complete==1)"
                )["mri_proc_dir"]
                .apply(uts.get_scan_tag)
                .tolist()
            )
            print(
                "\n  Just post-FreeSurfer SPM processing", "  " + ("-" * 35), sep="\n"
            )
            uts.print_list(mris_to_process_partial)

    if n_pet_scheduled > 0:
        tracers_to_process = raw_pets.loc[
            raw_pets["scheduled_for_processing"] == 1, "tracer"
        ].unique()
        print("\n  All PET scans scheduled to be processed:")
        for tracer in tracers_to_process:
            pets_to_process = (
                raw_pets.loc[
                    (raw_pets["scheduled_for_processing"] == 1)
                    & (raw_pets["tracer"] == tracer),
                    "pet_proc_dir",
                ]
                .apply(uts.get_scan_tag)
                .tolist()
            )
            print(f"\n  {tracer}", "  " + ("-" * len(tracer)), sep="\n")
            uts.print_list(pets_to_process)

    return raw_mris, raw_pets


def load_scan_typesf(scan_typesf):
    """Load scan types CSV and return a {name_in: name_out} dict."""
    scan_types = pd.read_csv(scan_typesf)
    scan_types["name_in"] = scan_types["name_in"].str.lower()
    scan_types = scan_types.drop_duplicates("name_in").dropna()
    scan_types = scan_types.set_index("name_in")["name_out"].to_dict()
    return scan_types


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


def get_subj(filepath, raw_dir):
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
    try:
        subj = filepath.replace(raw_dir + "/", "").split("/")[0]
        if len(subj) > 0:
            return subj
        else:
            return np.nan
    except IndexError:
        return np.nan


def get_scan_type(filepath, scan_types, raw_base="raw"):
    """Search filepath and return the scan type"""

    def find_matches(text, scan_types, scan_type_keys):
        """Helper function to find matches"""
        matches = [key for key in scan_type_keys if key in text]
        unique_values = list(set([scan_types[key] for key in matches]))
        return unique_values

    # Load the SCAN_TYPES dict and get lowercase keys
    scan_type_keys = list(scan_types)

    # Convert filepath to lowercase and get the basename
    filepath = filepath.lower()
    basename = op.basename(filepath)

    # First attempt to match in the basename
    match_values = find_matches(basename, scan_types, scan_type_keys)
    if len(match_values) == 1:
        return match_values[0]
    elif len(match_values) > 1:
        return "FILENAME MATCHED MULTIPLE PET TRACERS OR MRI MODALITIES"

    # Check for specific substring in basename
    if "coreg,_avg,_std_img_and_vox_siz,_uniform" in basename:
        return "FDG"

    # Attempt to match in the directory path after "raw"
    dirname = op.dirname(filepath)
    raw_base = f"/{raw_base}/"
    raw_idx = dirname.find(raw_base)
    if raw_idx == -1:
        return "FAILED TO IDENTIFY PET TRACER OR MRI MODALITY FROM FILENAME"

    # Only consider part of the path after "raw"
    search_path = dirname[raw_idx + len(raw_base) :]
    match_values = find_matches(search_path, scan_types, scan_type_keys)
    if len(match_values) == 1:
        return match_values[0]
    elif len(match_values) > 1:
        return "FILENAME MATCHED MULTIPLE PET TRACERS OR MRI MODALITIES"
    else:
        if "coreg,_avg,_std_img_and_vox_siz,_uniform" in search_path:
            return "FDG"
        else:
            return "FAILED TO IDENTIFY PET TRACER OR MRI MODALITY FROM FILENAME"


def get_scan_date(filepath, raw_base="raw"):
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
    raw_base = f"/{raw_base}/"
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


def get_orphan_mris(mri_scans, pet_scans, orphan_mri_tolerance=182):
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
                and ((today - x["mri_date"]).days > orphan_mri_tolerance)
                else 0
            ),
            axis=1,
        ),
    )
    return mri_scans


def audit_pet(
    pet_scans,
    audit_pet_res=True,
    expected_pet_res=6,
    audit_pet_to_mri_days=True,
    max_pet_to_mri=365,
    audit_repeat_mri=False,
):
    """Audit each PET scan and flag scans with potential issues

    Parameters
    ----------
    pet_scans : pd.DataFrame
        A DataFrame with columns 'subj', 'tracer', 'pet_date', 'pet_res',
        'mri_image_id', 'flag', and 'flag_notes'
    audit_pet_res : bool, optional
        If True, audit the PET resolution
    expected_pet_res : int, optional
        The expected PET resolution in mm
    audit_pet_to_mri_days : bool, optional
        If True, audit the days between PET and MRI scans
    max_pet_to_mri : int, optional
        The maximum allowed number of days between PET and MRI scans
    audit_repeat_mri : bool, optional
        If True, audit for repeated MRIs used for multiple PET scans

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
            pet_scans_cp["pet_res"] != expected_pet_res
        ].index.tolist()
        pet_scans_cp.loc[idx, "flag"] = 1
        pet_scans_cp.loc[
            idx, "flag_notes"
        ] += f"PET resolution could not be parsed or is not at {expected_pet_res}mm; "

    # Missing MRI
    idx = pet_scans_cp.loc[pd.isna(pet_scans_cp).any(axis=1)].index.tolist()
    pet_scans_cp.loc[idx, "flag"] = 1
    pet_scans_cp.loc[idx, "flag_notes"] += "Missing MRI; "

    # PET and MRI scan dates too far apart
    if audit_pet_to_mri_days:
        idx = pet_scans_cp.query(
            f"abs_days_mri_to_pet > {max_pet_to_mri}"
        ).index.tolist()
        pet_scans_cp.loc[idx, "flag"] = 1
        pet_scans_cp.loc[
            idx, "flag_notes"
        ] += f"Closest MRI is more than {max_pet_to_mri} days from PET date; "

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


def get_mri_proc_dir(subj, mri_date, proc_dir):
    """Return the processed MRI directory for each MRI scan"""
    return op.join(proc_dir, subj, "MRI-T1_{}".format(datetime_to_datestr(mri_date)))


def get_pet_proc_dir(subj, tracer, pet_date, proc_dir):
    """Return the processed PET directory for each PET scan"""
    return op.join(
        proc_dir,
        subj,
        "{}_{}".format(tracer, datetime_to_datestr(pet_date)),
    )


def check_if_freesurfer_run(mri_proc_dir, check_brainstem=True):
    """Return True if the MRI has been processed by Freesurfer"""
    freesurfer_dir = op.join(mri_proc_dir, "freesurfer")
    if not op.isdir(freesurfer_dir):
        return 0
    # Search for the Freesurfer output directory
    nuf = op.join(freesurfer_dir, "mri", "nu.mgz")
    aparcf = op.join(freesurfer_dir, "mri", "aparc+aseg.mgz")
    check_files = [nuf, aparcf]
    if check_brainstem:
        bstemf = op.join(
            freesurfer_dir, "mri", "brainstemSsLabels.v12.FSvoxelSpace.mgz"
        )
        check_files.append(bstemf)
    if all([op.isfile(f) for f in check_files]):
        return 1
    else:
        return 0


def check_if_mri_seg_complete(mri_proc_dir):
    """Return True if the MRI has been segmented"""
    if not op.isdir(mri_proc_dir):
        return 0

    # Search for the processed MRI files
    proc_files = {
        "seg_nu": glob(op.join(mri_proc_dir, "c*nu.nii")),
    }

    # Check if all file types are present
    if all([len(proc_files[key]) > 0 for key in proc_files]):
        return 1
    else:
        return 0


def check_if_mri_processed(mri_proc_dir):
    """Return True if the MRI has been fully processed"""
    if not op.isdir(mri_proc_dir):
        return 0

    # Search for the processed MRI files
    proc_files = {
        "affine_nu": glob(op.join(mri_proc_dir, "a*nu.nii")),
        "mask_brainstem": glob(op.join(mri_proc_dir, "*mask-brainstem.nii")),
        "mask_eroded_subcortwm": glob(
            op.join(mri_proc_dir, "*mask-eroded-subcortwm.nii")
        ),
        "mask_infcblgm": glob(op.join(mri_proc_dir, "*mask-infcblgm.nii")),
        "mask_pons": glob(op.join(mri_proc_dir, "*mask-pons.nii")),
        "mask_wcbl": glob(op.join(mri_proc_dir, "*mask-wcbl.nii")),
        "warp_nu": glob(op.join(mri_proc_dir, "w*nu.nii")),
    }

    # Check if all file types are present
    if all([len(proc_files[key]) > 0 for key in proc_files]):
        return 1
    else:
        return 0


def check_if_pet_processed(pet_proc_dir):
    """Return True if the PET scan has been fully processed"""
    if not op.isdir(pet_proc_dir):
        return 0

    # Get the scan info
    pet_tag = uts.get_scan_tag(pet_proc_dir)
    _, tracer, _ = uts.parse_scan_tag(pet_tag)

    # Get the list of reference regions used for each tracer
    code_dir = op.dirname(op.dirname(__file__))
    ref_regions = pd.read_csv(op.join(code_dir, "config", "ref_regions.csv"))
    ref_regions_by_tracer = {
        tracer: ref_regions.loc[ref_regions["tracer"] == tracer, "ref_region"].tolist()
        for tracer in ref_regions["tracer"].unique()
    }

    # List the files whose existence we will check for
    proc_files = {
        "native_pet_pet": op.join(pet_proc_dir, f"{pet_tag}.nii"),
        "native_mri_pet": op.join(pet_proc_dir, f"r{pet_tag}.nii"),
        "ref_region_means": op.join(pet_proc_dir, f"{pet_tag}_ref-region-means.csv"),
    }
    if tracer in AMYLOID_TRACERS:
        proc_files["cortical_summary_values"] = op.join(
            pet_proc_dir, f"{pet_tag}_amyloid-cortical-summary.csv"
        )
    for rr in ref_regions_by_tracer[tracer]:
        proc_files[f"native_mri_suvr_{rr}"] = op.join(
            pet_proc_dir, f"r{pet_tag}_suvr-{rr}.nii"
        )
        proc_files[f"warped_mni_suvr_{rr}"] = op.join(
            pet_proc_dir, f"wr{pet_tag}_suvr-{rr}.nii"
        )
        proc_files[f"affine_mni_suvr_{rr}"] = op.join(
            pet_proc_dir, f"ar{pet_tag}_suvr-{rr}.nii"
        )
        proc_files[f"native_mri_suvr_extractions_{rr}"] = op.join(
            pet_proc_dir, f"r{pet_tag}_suvr-{rr}_roi-extractions.csv"
        )

    # Check if all files are present
    if all([op.isfile(f) for f in proc_files.values()]):
        return 1
    else:
        return 0


def get_mri_to_process(
    mri_scan_number,
    mri_is_orphan,
    already_processed,
    baseline_processed,
    overwrite,
    schedule_all_mris,
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
        # Don't process orphan MRIs unless schedule_all_mris is True
        if schedule_all_mris or not mri_is_orphan:
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


def save_raw_mri_index(mri_scans, scans_to_process_dir, timestamp):
    """Save the MRI scan index to a CSV file"""
    # Move any existing files to archive
    archive_dir = op.join(scans_to_process_dir, "archive")
    files = glob(op.join(scans_to_process_dir, "raw_MRI-T1_index*.csv"))
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
    outf = op.join(scans_to_process_dir, f"raw_MRI-T1_index_{timestamp}.csv")

    mri_scans.to_csv(outf, index=False)
    print(f"  * Saved raw MRI scan index to {outf}")


def save_raw_pet_index(pet_scans, scans_to_process_dir, timestamp):
    """Save the PET scan index to a CSV file"""
    # Move any existing files to archive
    archive_dir = op.join(scans_to_process_dir, "archive")
    files = glob(op.join(scans_to_process_dir, "raw_PET_index*.csv"))
    if files:
        if not op.isdir(archive_dir):
            os.makedirs(archive_dir)
        for f in files:
            os.rename(f, op.join(archive_dir, op.basename(f)))

    # Save the pet_scans dataframe
    outf = op.join(scans_to_process_dir, f"raw_PET_index_{timestamp}.csv")
    pet_scans.to_csv(outf, index=False)
    print(f"  * Saved raw PET scan index to {outf}")


def _parse_args():
    """Parse and return command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "High-level program to schedule MRI and PET scans for later processing\n\n"
            + "Overview\n--------\n"
            + "This program is designed for an MRI-based PET processing pipeline in which\n"
            + "scans are downloaded from LONI (with each PET scan being a preprocessed\n"
            + "'Step 4' image in common resolution).\n\n"
            + "Before running this program, a zip file containing 1+ scans must be added by\n"
            + "the user to the 'newdata' directory, unzipped, converted from DICOM to NIfTI\n"
            + "format, and moved to the 'raw' directory (see 'setup_leads_processing.m' to\n"
            + "see how these last three steps are automated).\n\n"
            + "Presumably, the user is normally engaging with this program through the\n"
            + "run_leads_mri_based_pipeline.m script that calls it. That said, it is\n"
            + "possible to run this program as a standalone from the Linux command line to\n"
            + "have a more granular level of control.\n\n"
            + "The following steps are completed:\n"
            + "  1.  All scan directories in 'raw' are identified by a recursive search for\n"
            + "      *.nii files\n"
            + "  2.  Filenames are parsed to resolve:\n"
            + "      -  Subject ID\n"
            + "      -  Scan type (MRI or PET; and for PET, which tracer)\n"
            + "      -  Scan acquisition date\n"
            + "      -  LONI Image ID\n"
            + "      -  For PET scans, the spatial resolution\n"
            + "  3.  Each PET scan is matched to the closest MRI scan\n"
            + "  4.  Orphan MRIs (those not recently acquired and not closest to any PET\n"
            + "      scan) are identified, and under default settings will not be scheduled\n"
            + "      for processing\n"
            + "  5.  PET scans are audited for potential issues that would prevent processing.\n"
            + "      These include:\n"
            + "      - Inability to parse any of the scan information above\n"
            + "      - PET not at the expected resolution (6mm by default)\n"
            + "      - PET and MRI scans too far apart (>365 days, by default)\n"
            + "      - Same MRI would be used for multiple PET scans of the same tracer (e.g.\n"
            + "        two FTP timepoints)\n"
            + "  6.  A search is conducted to identify which MRIs in 'raw' have been\n"
            + "      processed by FreeSurfer, and which have completed post-FreeSurfer,\n"
            + "      SPM-based processing, respectively\n"
            + "  7.  A search is conducted to identify which PET scans in 'raw' have been\n"
            + "      fully processed\n"
            + "  8.  MRIs are scheduled for processing if they have not been fully processed,\n"
            + "      are not orphans, and (for follow-up MRIs) their corresponding baseline\n"
            + "      MRI has been processed\n"
            + "  9.  PET scans are scheduled for processing if they have not been fully\n"
            + "      processed, are not flagged with an issue, and their corresponding MRI\n"
            + "      has been fully processed\n"
            + "  10. The user is informed of how many MRI and PET scans are in 'raw', and how\n"
            + "      many of these scans have completed processing, have been flagged with an\n"
            + "      issue, and have been scheduled for processing\n"
            + "  11. Two CSV files containing scan-level information are saved to\n"
            + "      'scans_to_process': one for MRIs and one for PET scans. These files are\n"
            + "      read by downstream scripts in the processing pipeline. Manually editing\n"
            + "      these files will affect which scans the scripts attempt to process,\n"
            + "      though could also introduce errors, so edit the CSVs with caution and at\n"
            + "      your own risk.\n\n"
            + "The remainder of this documentation lists the command line arguments that can\n"
            + "be passed to this program to change some of its default behavior."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        exit_on_error=False,
    )
    parser.add_argument(
        "-p",
        "--proj-dir",
        required=True,
        help=(
            "Full path to the top-level project directory. Must contain\n"
            + "'data/raw' and 'data/processed' subdirectories, along with\n"
            + "'metadata/scans_to_process' where output CSVs will be saved"
        ),
    )
    parser.add_argument(
        "--orphan-tol",
        type=int,
        default=182,
        dest="orphan_mri_tolerance",
        help=(
            "Defines the maximum number of days from present that an MRI\n"
            + "can have been acquired without being linked to a PET scan\n"
            + "before it is considered an orphan and is not scheduled for\n"
            + "processing when --all is not specified"
        ),
    )
    parser.add_argument(
        "-a",
        "--all",
        dest="schedule_all_mris",
        action="store_true",
        help="Schedule all MRIs in 'raw' for processing, even orphans",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help=(
            "Schedule all eligible PET and MRI scans in 'raw' to be\n"
            + "processed, regardless of whether they have already been\n"
            + "processed. Orphan MRIs are still not scheduled unless --all\n"
            + "is specified, and PET scans that are flagged with an issue\n"
            + "are not scheduled"
        ),
    )
    parser.add_argument(
        "--no-audit-pet-res",
        action="store_false",
        dest="audit_pet_res",
        help=(
            "Don't audit the PET resolution. By default, PET scans are\n"
            + "flagged if their resolution is not 6mm. This flag disables\n"
            + "that check"
        ),
    )
    parser.add_argument(
        "--expected-pet-res",
        type=int,
        default=6,
        dest="expected_pet_res",
        help=(
            "The expected PET resolution in mm. By default, PET scans are\n"
            + "flagged if their resolution is not 6mm"
        ),
    )
    parser.add_argument(
        "--no-audit-pet-to-mri-days",
        action="store_false",
        dest="audit_pet_to_mri_days",
        help=(
            "Don't audit the days between PET and MRI scans. By default,\n"
            + "PET scans are flagged if the closest MRI is more than 365\n"
            + "days away. This flag disables that check"
        ),
    )
    parser.add_argument(
        "--max-pet-to-mri",
        type=int,
        default=365,
        help=(
            "The maximum allowed number of days between a PET scan and\n"
            + "the closest MRI scan. By default, PET scans are flagged if\n"
            + "the closest MRI is more than 365 days away"
        ),
    )
    parser.add_argument(
        "--audit-repeat-mri",
        action="store_true",
        dest="audit_repeat_mri",
        help=(
            "Audit for repeated MRIs used for multiple PET scans.\n"
            + "When this argument is passed, PET scans are flagged if the same\n"
            + "MRI would need to be used to process more than one PET scan\n"
            + "from the same subject and tracer"
        ),
    )
    parser.add_argument(
        "--no-check-brainstem",
        action="store_false",
        dest="check_brainstem",
        help=(
            "Don't check for the presence of the brainstemSsLabels.v12\n"
            + "file in the Freesurfer output directory. By default, this\n"
            + "file is checked and if missing, the MRI is not scheduled\n"
            + "for processing"
        ),
    )
    parser.add_argument(
        "--no-save-csv",
        action="store_false",
        dest="save_csv",
        help=(
            "Parse 'raw' and view the summary information for this\n"
            + "program, but don't save new CSV files (i.e., don't schedule\n"
            + "new scans for processing yet, but see how many scans would\n"
            + "be scheduled if the program is rerun without --view-only"
        ),
    )

    # Parse the command line arguments
    return parser.parse_args()


if __name__ == "__main__":
    # Get command line arguments.
    args = _parse_args()

    # Call the main function
    main(
        proj_dir=args.proj_dir,
        orphan_mri_tolerance=args.orphan_mri_tolerance,
        schedule_all_mris=args.schedule_all_mris,
        overwrite=args.overwrite,
        audit_pet_res=args.audit_pet_res,
        expected_pet_res=args.expected_pet_res,
        audit_pet_to_mri_days=args.audit_pet_to_mri_days,
        max_pet_to_mri=args.max_pet_to_mri,
        audit_repeat_mri=args.audit_repeat_mri,
        check_brainstem=args.check_brainstem,
        save_csv=args.save_csv,
    )

    # Exit successfully
    sys.exit(0)
