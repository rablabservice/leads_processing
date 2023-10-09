#!/usr/bin/env python

"""
$ update_reads.py main_ssheet reads_topdir

"""

import sys
import os.path as op
import pandas as pd


def update_reads_df(main_ssheet, reads_topdir, verbose=True):
    """Update visual reads ssheet with info from individual read files.

    Parameters
    ----------
    main_ssheet : str
        Path to the main spreadsheet with visual reads for all subjects,
        which this code will update.
    reads_topdir : str
        Path to the parent directory in which individual read CSVs are
        contained. The expected file hierarchy is strict:
        [reads_topdir]/[subject]/results/reads/[subject]_read.csv
    verbose : bool, optional
        Whether to print updates to the terminal. Default is True.
    """
    # Check that the user input paths exist.
    assert op.isfile(main_ssheet)
    assert op.isdir(reads_topdir)

    # Note which columns of the visual reads spreadsheet should be
    # checked for missing values, and which columns should be updated
    # with values from the individual visual read files.
    check_cols = ["date_of_read", "read_result"]
    update_cols = [
        "date_of_read",
        "metaroi_suvr",
        "read_result",
        "read_notes",
        "read_quant_disagreement",
    ]

    # Get a dict of visual readers and their initials.
    readers = _get_readers()

    # Open the main visual reads spreadsheet.
    reads_df = pd.read_excel(main_ssheet)
    if verbose:
        print("Opened {}".format(main_ssheet))

    # Find scans with missing values in `check_cols` columns. These are
    # the scans that we will attempt to update based on info in the
    # individual visual read files.
    df_incomplete = reads_df[reads_df[check_cols].isnull().any(axis=1)]
    if verbose:
        print(
            "Found {}/{} scans with missing reads in {}".format(
                len(df_incomplete), len(reads_df), op.basename(main_ssheet)
            )
        )

    # Iterate over scans in df_incomplete.
    if verbose:
        print("\n", "-" * 30)
    count = 0
    for idx, row in df_incomplete.iterrows():
        # Get filepath to the current scan's visual read file.
        reads_csvf = op.join(
            reads_topdir, row["id"], "results", "read", "{}_read.csv".format(row["id"])
        )

        # Check that the visual read file exists.
        if not op.isfile(reads_csvf):
            if verbose:
                print("{} - FILE NOT FOUND >> {}".format(row["id"], reads_csvf))
            continue

        # Try to parse the visual read file.
        try:
            reads_csv = pd.read_csv(reads_csvf).iloc[0]
        except UnicodeDecodeError:
            print("Can't read {} - UnicodeDecodeError".format(reads_csvf))
            continue

        # Make all column names lowercase.
        reads_csv.index = reads_csv.index.str.lower()

        # Check if the assigned reader has completed the read (we just
        # check if they added their name to the `reader` column in
        # the visual read file).
        if pd.isna(reads_csv["reader"]):
            if verbose:
                print(
                    "{} - {} has not completed the read".format(
                        row["id"], row["reader"]
                    )
                )
            continue
        # Check if any of the required columns have not been filled out.
        if pd.isna(reads_csv).sum() > 0:
            if verbose:
                print(
                    "{} - Read partially complete but missing the following fields:".format(
                        row["id"]
                    ),
                    *reads_csv[pd.isna(reads_csv)].index.tolist(),
                    sep="\n" + " " * (len(row["id"]) + 3)
                )
        # Check if the assigned reader is who actually completed the
        # read.
        if (
            reads_csv["reader"].strip().lower()
            not in readers[row["reader"].strip().lower()]
        ):
            if verbose:
                print(
                    "{} - Different readers listed in {} and {}:".format(
                        row["id"], op.basename(main_ssheet), op.basename(reads_csvf)
                    ),
                    "{} != {}".format(row["reader"], reads_csv["reader"]),
                    sep="\n" + " " * (len(row["id"]) + 3),
                )

        # For all columns in `update_cols`, update empty values in the
        # main visual reads spreadsheet with info from the corresponding
        # field of the visual read CSV. Does not overwrite any existing
        # values in the main visual reads spreadsheet.
        _count = 0
        for col in update_cols:
            col = col.strip().lower()
            if col not in reads_csv:
                print("{} - {} not found in {}".format(row["id"], col, reads_csvf))
                continue
            if pd.isna(reads_csv[col]):
                continue
            if not pd.isna(reads_df.loc[idx, col]):
                if verbose:
                    print(
                        "{} - {} value already exists in {}; will not overwrite".format(
                            row["id"], col, op.basename(main_ssheet)
                        )
                    )
                continue
            else:
                reads_df.loc[idx, col] = reads_csv[col]
                _count += 1
        if _count > 0:
            count += 1

    # Save the updated visual reads spreadsheet.
    if count > 0:
        reads_df.to_excel(main_ssheet, index=False)
        if verbose:
            print("-" * 30, "\n")
            print(
                "Saved updated info for {} scans in {}".format(
                    count, op.basename(main_ssheet)
                )
            )
    else:
        if verbose:
            print("-" * 30, "\n")
            print("Nothing to update!")

    return None


def _get_readers():
    """Get a list of readers from the visual reads spreadsheet."""
    readers = {
        "charles": ["charles", "ccw", "cw"],
        "courtney": ["courtney", "clh", "ch"],
        "david": ["david", "dsm", "dm", "dns"],
        "ehud": ["ehud", "ez", "e-z"],
        "gil": ["gil", "gdr", "gr"],
        "hanna": ["hanna", "hc"],
    }
    return readers


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            __doc__,
            update_reads_df.__doc__,
            sep="\n",
        )
        exit()

    # Parse arguments.
    if op.isabs(sys.argv[1]):
        main_ssheet = sys.argv[1]
    else:
        main_ssheet = op.join(op.getcwd(), sys.argv[1])
    if op.isabs(sys.argv[2]):
        reads_topdir = sys.argv[2]
    else:
        reads_topdir = op.join(op.getcwd(), sys.argv[2])

    # Update the reads_df.
    update_reads_df(main_ssheet, reads_topdir)
