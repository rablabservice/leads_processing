#!/usr/bin/env python

import os
import os.path as op
import datetime
from glob import glob
import numpy as np
import pandas as pd

# Define directory paths
PATHS = {
    "proj": "/mnt/coredata/processing/leads",
}
PATHS["data"] = op.join(PATHS["proj"], "data")
PATHS["freesurfer"] = op.join(PATHS["data"], "freesurfer")
PATHS["metadata"] = op.join(PATHS["proj"], "metadata")
PATHS["raw"] = op.join(PATHS["data"], "raw")
PATHS["processed"] = op.join(PATHS["data"], "processed")


# Define functions
def scrape_raw(raw_dir, verbose=True):
    """Scrape raw directory for all nifti files and parse them.

    Returns a pandas DataFrame with columns:
    - subj: subject ID
    - scan_date: scan date (YYYY-MM-DD)
    - scan_type: scan type (MRI modality or PET tracer)
    - raw_petf: full path to the nifti file in raw
    """
    def _get_subj(filepath):
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
            subj = filepath.replace(raw_dir + "/", "").split("/")[1]
            if len(subj) > 0:
                return subj
            else:
                return np.nan
        except IndexError:
            if verbose:
                print(f"WARNING: Could not parse subject ID from {filepath}")
            return np.nan

    def _get_scan_date(filepath):
        """Return the scan date from filepath to the recon'd nifti.

        Iterates over filepath directories from right to left until it finds
        a filename or directory whose first 10 characters matches the date
        format YYYY-MM-DD.

        Returns np.nan if no scan date is found, otherwise a string like
        'YYYY-MM-DD'.
        """
        for d in filepath.split(op.sep)[::-1]:
            try:
                acqdate = check_dt_fmt(d[:10], raise_error=True)
                return acqdate
            except ValueError:
                pass
        if verbose:
            print(f"WARNING: Could not parse scan date from {filepath}")
        return np.nan

    def _get_scan_type(filepath, scan_type_map_file=None):
        """Parse the filepath and return the scan type."""
        if scan_type_map_file is None:
            scan_type_map_file = op.join(
                PATHS["metadata"], "ssheets", "scan_types_and_tracers.csv"
            )
        scan_type_map = (
            pd.read_csv(scan_type_map_file).set_index("name_in")["name_out"].to_dict()
        )
        basename = op.basename(filepath).lower()
        for k, v in scan_type_map.items():
            if k in basename:
                return v
        if "fdg" in filepath.lower():
            return "FDG"
        if verbose:
            print(f"WARNING: Could not parse scan type from {filepath}")
        return np.nan

    # Find all nifti PET files in raw
    search_raw_dirs = [d for d in glob(op.join(PATHS["raw"], "*")) if (op.isdir(d) and not (d.split(op.sep)[-1] == "mri"))]
    raw_scans = []
    for d in search_raw_dirs:
        raw_scans += glob(op.join(d, "**", "*.nii"), recursive=True)

    # Parse scan info for each nifti
    output = []
    for scanf in raw_scans:
        subj = _get_subj(scanf, raw_dir)
        scan_date = _get_scan_date(scanf)
        scan_type = _get_scan_type(scanf)
        output.append([subj, scan_date, scan_type, scanf])

    cols = ["subj", "scan_date", "scan_type", "raw_petf"]
    output = pd.DataFrame(output, columns=cols)
    return output


def find_closest_mri(
    subj, scan_date, freesurfer_dir, limit_days=365, strict_limit=False, verbose=True
):
    """Return closest MRI date, days from PET, and Freesurfer path.

    Parameters
    ----------
    subj : str
        The subject ID.
    scan_date : str
        Scan date (YYYY-MM-DD) to match the closest MRI scan to.
    freesurfer_dir : str
        Path to the top-level freesurfer directory containing individual
        processed MRI directories like <subj>_<scan_date>.
    limit_days : int
        A warning is raised if no MRI scan is found within limit days
        and np.nan is returned.
    strict_limit : bool
        If True, np.nan is returned if no MRI scan is found within limit days.
    verbose : bool
        If True, print warnings if no MRI scan is found within limit days.
    """
    proc_mris = glob(op.join(freesurfer_dir, f"{subj}_*"))
    if len(proc_mris) == 0:
        if verbose:
            print(
                (f"WARNING: {subj} scan on {scan_date} has no processed MRI scans in " +
                 f"{freesurfer_dir}")
            )
        return np.nan, np.nan, np.nan

    proc_mri_dates = [check_dt_fmt(op.basename(p).split("_")[1]) for p in proc_mris]

    days_to_scan = []
    for d in proc_mri_dates:
        days_to_scan.append(
            date_diff(datestr_to_datetime(d), datestr_to_datetime(scan_date), abs=True)
        )
    closest_mri = proc_mris[np.argmin(days_to_scan)]
    closest_mri_date = proc_mri_dates[np.argmin(days_to_scan)]
    min_days = min(days_to_scan)

    if min_days > limit_days:
        if verbose:
            print(
                (f"WARNING: {subj} scan on {scan_date} has no matching MRI within " +
                 f"{limit_days} days. Closest MRI is {min_days} days away")
            )
        if strict_limit:
            return np.nan, np.nan, np.nan

    return closest_mri_date, min_days, closest_mri


def date_diff(date1, date2, abs=False):
    """Return date2 - date1 in days."""
    try:
        diff = (date2 - date1).days
        if abs:
            return np.abs(diff)
        else:
            return diff
    except TypeError:
        return np.nan


def check_dt_fmt(datestr, raise_error=False):
    """Return datestr if formatted like YYYY-MM-DD.

    If raise_error is True, raise a ValueError if datestr is not
    formatted like YYYY-MM-DD. Otherwise return np.nan.
    """
    try:
        datestr_to_datetime(datestr)
        return datestr
    except ValueError:
        if raise_error:
            raise ValueError(f"{datestr} is not formatted like YYYY-MM-DD")
        else:
            return np.nan


def datestr_to_datetime(datestr):
    """Convert a date string to a datetime object."""
    return datetime.datetime.strptime(datestr, "%Y-%m-%d")


def datetime_to_datestr(dt):
    """Convert a datetime object to a date string."""
    return dt.strftime("%Y-%m-%d")


def now():
    """Return the current date and time down to seconds."""
    return datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")


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


def main():
    # Scrape the raw directory for all PET niftis
    raw_scans = scrape_raw(PATHS["raw"])

    # Search Freesurfer directories for the closest MRI to each PET scan
    raw_scans["mri_date"], raw_scans["days_to_mri"], raw_scans["fs_dir"] = zip(
        *raw_scans.apply(
            lambda x: find_closest_mri(x["subj"], x["scan_date"], PATHS["freesurfer"]),
            axis=1,
        )
    )

    # Add paths to processed MRI-T1 directories that match the selected
    # Freesurfer directory for each PET scan
    # Add processed MRI directories to raw_scans
    link_to_old_leads_mris = True
    mri_dirs = []
    for idx, scan in raw_scans.iterrows():
        if scan["fs_dir"] is np.nan:
            mri_dirs.append(np.nan)
            continue

        # Get the name of the processed MRI directory
        subj = scan["subj"]
        mri_date = scan["fs_dir"].split("_")[1]
        proc_mri_dir = op.join(PATHS["processed"], subj, f"MRI-T1_{mri_date}")

        # Symlink to processed MRI directories in the old LEADS database
        if link_to_old_leads_mris:
            proc_dir_old = "/mnt/coredata/Projects/LEADS/data_f7p1/processed"
            for tp in range(1, 6):
                proc_mri_dir_old = op.join(
                    proc_dir_old, subj, f"Timepoint{tp}", f"MRI_T1_{mri_date}"
                )
                if op.isdir(proc_mri_dir_old):
                    os.makedirs(op.dirname(proc_mri_dir), exist_ok=True)
                    if not op.exists(proc_mri_dir):
                        os.symlink(proc_mri_dir_old, proc_mri_dir)
                    break
                if tp == 4:
                    proc_mri_dir_old = np.nan

        # Store a reference to the processed MRI directory, if it exists
        if op.exists(proc_mri_dir):
            mri_dirs.append(proc_mri_dir)
        else:
            mri_dirs.append(np.nan)

    raw_scans["mri_dir"] = mri_dirs

    # Drop raw_scans rows with missing data (e.g. if an MRI is not
    # available or has not been processed). Then sort raw_scan rows.
    raw_scans = (
        raw_scans.dropna()
        .sort_values(["subj", "scan_type", "scan_date"])
        .reset_index(drop=True)
    )

    # Convert days_to_mri back to int
    raw_scans["days_to_mri"] = raw_scans["days_to_mri"].astype(int)

    # Track the visit number of each PET scan, with visit 1 being the
    # earliest scan date for a given subject and scan type, visit 2
    # being the next earliest scan date, and so on.
    raw_scans["visit"] = raw_scans.groupby(["subj", "scan_type"]).cumcount() + 1
    cols = raw_scans.columns.tolist()
    cols.insert(cols.index("scan_type") + 1, cols.pop(cols.index("visit")))
    raw_scans = raw_scans[cols]

    # Add a diagnosis column that describes each subject's cohort
    # assignment (EOAD, EOnonAD, or CN)
    dxf = op.join(PATHS["metadata"], "ssheets", "LEADS_Internal_PET-Screening.xlsx")
    if op.isfile(dxf):
        dx = pd.read_excel(dxf)
        dx_map = {"ID": "subj", "Cohort": "dx"}
        dx = dx.rename(columns=dx_map)[["subj", "dx"]]
        raw_scans = dx.merge(raw_scans, on="subj", how="right")

    subj_regf = op.join(
        PATHS["metadata"], "ssheets", "Participant Registration_vertical.csv"
    )
    if op.isfile(subj_regf):
        subj_reg = pd.read_csv(subj_regf)
        subj_reg_map = {
            "subject.label": "subj", "dd_revision_field.translated_value": "dx"
        }
        subj_reg = subj_reg.rename(columns=subj_reg_map)[["subj", "dx"]]
    cn_subjs = subj_reg.query("(dx=='Cognitively Normal Participant')")["subj"].tolist()
    raw_scans.loc[pd.isna(raw_scans["dx"]), "dx"] = raw_scans.loc[
        pd.isna(raw_scans["dx"]), "subj"
    ].apply(lambda x: "CN" if np.isin(x, cn_subjs) else np.nan)

    # Find which PET scans are ready and needing to be processed
    proc_pet_dirs = []
    need_to_process = []
    for idx, scan in raw_scans.iterrows():
        if scan["mri_dir"] is np.nan:
            proc_pet_dirs.append(np.nan)
            need_to_process.append(False)
            continue

        proc_pet_dir = op.join(
            PATHS["processed"], scan["subj"], f"{scan['scan_type']}_{scan['scan_date']}"
        )
        proc_pet_dirs.append(proc_pet_dir)
        need_to_process.append(not op.exists(proc_pet_dir))

    raw_scans["proc_pet_dir"] = proc_pet_dirs
    raw_scans["need_to_process"] = need_to_process

    # Save the raw_scans dataframe to a CSV file
    outf = op.join(PATHS["metadata"], "log", f"raw_pet_scans_{now()}.csv")
    raw_scans.to_csv(outf, index=False)
    print(f"Found {len(raw_scans):,} PET scans in {PATHS["raw"]}:")
    for tracer in raw_scans["scan_type"].unique():
        n = len(raw_scans.query(f"scan_type=='{tracer}'"))
        print(f"  {n:>3} {tracer}")
    print(f"Saved raw_scans to {outf}")
