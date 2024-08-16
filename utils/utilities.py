#!/usr/bin/env python

"""
A collection of functions that are used across multiple subdirectories
within the processing pipeline.
"""

import datetime
import os.path as op
from glob import glob


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


def glob_sort(pattern):
    """Return files matching pattern in alphanumeric order.

    Returns
    -------
    files : list of str
        List of files matching pattern, sorted alphanumerically.
    """
    files = sorted(glob(pattern))
    return files


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


def get_scan_tag(scan_dir):
    """Return the scan tag from the processed scan directory."""
    scan_dir = op.abspath(scan_dir)
    if op.isfile(scan_dir):
        scan_dir = op.dirname(scan_dir)
    subj = op.basename(op.dirname(scan_dir))
    scan_type = op.basename(scan_dir).split("_")[0]
    scan_date = op.basename(scan_dir).split("_")[-1]
    scan_tag = f"{subj}_{scan_type}_{scan_date}"
    return scan_tag


def parse_scan_tag(scan_tag):
    """Return subject ID, scan type, and scan date from the scan tag."""
    subj, scan_type, scan_date = scan_tag.split("_")
    return subj, scan_type, scan_date


def print_list(lst, max_line=115, sep="  "):
    """Print a list of strings to the console"""
    current_line = ""
    for entry in lst:
        new_entry = entry
        new_line = f"{current_line}{sep}{new_entry}"
        if len(new_line) > max_line:
            print(current_line)
            current_line = f"{sep}{new_entry}"
        else:
            current_line = new_line
    if len(current_line) > 0:
        print(current_line)
    print()
