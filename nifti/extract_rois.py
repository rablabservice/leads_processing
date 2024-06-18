#!/usr/bin/env python

"""
$ extract_rois.py -i [pet1.nii pet2.nii ...] -m [mask1.nii mask2.nii ...]
OR
$ extract_rois.py -i [pet1.nii pet2.nii ...] -a [aparc+aseg1.nii aparc+aseg2.nii ...]
"""

import sys
import os.path as op
import importlib.resources
import argparse
from collections import OrderedDict as od
from inspect import isroutine
import re
from time import time

import nibabel as nib
import numpy as np
import pandas as pd
import general.nifti.nifti_ops as nops


class Timer(object):
    """I say how long things take to run."""

    def __init__(self, msg=None):
        """Start the global timer."""
        self.reset()
        if msg is None:
            self.msg = "Time elapsed: "
        else:
            self.msg = msg

    def __str__(self):
        """Print how long the global timer has been running."""
        elapsed = self.check()
        hours = int(elapsed / 3600)
        minutes = int((elapsed % 3600) / 60)
        seconds = elapsed % 60
        if hours > 0:
            msg = self.msg + "{}h, {}m, {:.3f}s".format(hours, minutes, seconds)
        elif minutes > 0:
            msg = self.msg + "{}m, {:.3f}s".format(minutes, seconds)
        else:
            msg = self.msg + f"{elapsed}s"
        return msg

    def check(self, reset=False):
        """Report the global runtime."""
        runtime = time() - self.start
        if reset:
            self.reset()
        return runtime

    def loop(self, key=None, verbose=True):
        """Report the loop runtime and reset the loop timer."""
        if not hasattr(self, "loops"):
            self.loops = od([])
        if not hasattr(self, "last_loop_start"):
            self.last_loop_start = self.start
        if key is None:
            key = "loop {}".format(len(self.loops) + 1)

        loop_runtime = time() - self.last_loop_start
        self.loops[key] = loop_runtime
        self.last_loop_start = time()
        if verbose:
            print("{}: {:.1f}s".format(key, self.loops[key]))

    def reset(self):
        """Reset the global timer."""
        self.start = time()


class TextFormatter(argparse.RawTextHelpFormatter):
    """Custom formatter for argparse help text."""

    # use defined argument order to display usage
    def _format_usage(self, usage, actions, groups, prefix):
        if prefix is None:
            prefix = "usage: "

        # if usage is specified, use that
        if usage is not None:
            usage = usage % dict(prog=self._prog)

        # if no optionals or positionals are available, usage is just prog
        elif usage is None and not actions:
            usage = "%(prog)s" % dict(prog=self._prog)
        elif usage is None:
            prog = "%(prog)s" % dict(prog=self._prog)
            # build full usage string
            action_usage = self._format_actions_usage(actions, groups)  # NEW
            usage = " ".join([s for s in [prog, action_usage] if s])
            # omit the long line wrapping code
        # prefix with 'usage:'
        return "%s%s\n\n" % (prefix, usage)


def _fmt_long_str(x, maxlen=50):
    """Truncate long strings of semicolon-separated values."""
    if len(x) <= maxlen:
        return x
    elif len(x) > maxlen:
        stop = x[maxlen + 4 :].find(";")
        if stop == -1:
            return x
        else:
            # find the last ';'
            start_last = x.rfind(";") + 1
            return x[: stop + maxlen] + "..." + x[start_last:]


def load_rois(roi_file):
    """Load dictionary of ROI names to lists of 1+ int labels."""
    rois = pd.read_csv(roi_file)
    rois = od(zip(rois.iloc[:, 0], rois.iloc[:, 1]))
    try:
        rois = {
            k: list(np.unique([int(x) for x in v.split(";")])) for k, v in rois.items()
        }
    except AttributeError:
        pass
    rois = pd.Series(rois)
    rois = pd.Series(rois)
    return rois


def roi_desc(dat, rois, subrois=None, aggf=np.mean, conv_nan=0):
    """Apply `aggf` over `dat` values within each ROI mask.

    Parameters
    ----------
    dat :
        Filepath string, nifti image, or array-like object.
    rois : str, list[str], or dict-like {str: obj}
        Map each ROI name to its filepath string(s), nifti image, or
        array.
    subrois : dict of {str: int or list}
        Map each sub-ROI within the main ROI mask to a value or list of
        mask values that comprise it. The classic example is of an
        aparc+aseg file containing multiple regions with different
        labels. Note: subrois cannot be passed if len(rois) > 1.
    aggf : function, list of functions, or dict of functions
        Function or functions to apply over `dat` values within each
        ROI.
    conv_nan : bool, number, or NoneType object
        Convert NaNs in `dat` to `conv_nan`. No conversion is applied if
        `conv_nan` is np.nan, None, or False.

    Returns
    -------
    output : DataFrame
        `aggf` output for each agg function, for each ROI. Index is the
        ROI names, columns are the function names. The last column is
        ROI volume (number of voxels in the mask).
    """
    if (not isinstance(rois, str)) and (len(rois) > 1) and (subrois is not None):
        raise ValueError("Cannot define multiple rois and subrois")

    # Load the data array.
    dat = load_nii(dat, dat_only=True, flatten=True, conv_nan=conv_nan)

    # Format the ROIs to be dict-like.
    if isinstance(rois, str):
        rois = [rois]

    if isinstance(rois, (list, tuple)):
        rois_dict = od([])
        for roi in rois:
            splits = splitt(roi, ["_", "."])
            for ii, string in enumerate(splits):
                if string.startswith("mask-"):
                    roi_name = "-".join(string.split("-")[1:])
                    rois_dict[roi_name] = roi
                    break
                elif ii == len(splits) - 1:
                    rois_dict[".".join(op.basename(roi).split(".")[:-1])] = roi
        rois = rois_dict
    elif hasattr(rois, "keys"):
        pass
    else:
        raise ValueError("rois must be str, list, tuple, or dict-like")

    # Format the aggregation functions to be dict-like.
    if isroutine(aggf):
        aggf = od({aggf.__name__: aggf})
    elif not isinstance(aggf, dict):
        aggf = od({func.__name__: func for func in aggf})

    # Prepare the output DataFrame.
    if subrois is not None:
        output_idx = list(subrois.keys())
    else:
        output_idx = list(rois.keys())
    output_cols = list(aggf.keys()) + ["voxels"]
    output = pd.DataFrame(index=output_idx, columns=output_cols)
    output = output.rename_axis("roi")

    # Loop over the ROIs and sub-ROIs.
    for roi, roi_mask in rois.items():
        if subrois is not None:
            mask = load_nii(roi_mask, dat_only=True, flatten=True, binarize=False)
            assert dat.shape == mask.shape
            for subroi, subroi_vals in subrois.items():
                mask_idx = np.where(np.isin(mask, subroi_vals))
                for func_name, func in aggf.items():
                    output.at[subroi, func_name] = func(dat[mask_idx])
                output.at[subroi, "voxels"] = mask_idx[0].size
        else:
            mask = load_nii(roi_mask, dat_only=True, flatten=True, binarize=True)
            assert dat.shape == mask.shape
            mask_idx = np.where(mask)
            for func_name, func in aggf.items():
                output.at[roi, func_name] = func(dat[mask_idx])
            output.at[roi, "voxels"] = mask_idx[0].size

    return output


def load_nii(
    infile,
    dtype=np.float32,
    squeeze=True,
    flatten=False,
    conv_nan=0,
    binarize=False,
    int_rounding="nearest",
):
    """Load a NIfTI file and return the NIfTI image and data array.

    Returns (img, dat), with dat being an instance of img.dataobj loaded
    from disk. You can modify or delete dat and get a new version from
    disk: ```dat = np.asanyarray(img.dataobj)```

    Parameters
    ----------
    infile : str
        The nifti file to load.
    dtype : data-type
        Determines the data type of the data array returned.
    flatten : bool
        If true, `dat` is returned as a flattened copy of the
        `img`.dataobj array. Otherwise `dat`.shape == `img`.shape.
    conv_nan : bool, number, or NoneType object
        Convert NaNs to `conv_nan`. No conversion is applied if
        `conv_nan` is np.nan, None, or False.
    binarize : bool
        If true, `dat` values > 0 become 1 and all other values are 0.
        `dat` type is recast to np.uint8.
    int_rounding : str
        Determines how the data array is recast if `binarize` is false
        and `dtype` is an integer.
        `nearest` : round to the nearest integer
        `floor` : round down
        `ceil` : round up

    Returns
    -------
    img : Nifti1Image
    dat : ndarray or ndarray subclass
    """
    # Get the right file extension.
    infile = find_gzip(infile)

    # Load the NIfTI image and data array.
    img = nib.load(infile)
    dat = np.asanyarray(img.dataobj)

    # Format the data array.
    dat = _format_array(
        dat,
        dtype=dtype,
        squeeze=squeeze,
        flatten=flatten,
        conv_nan=conv_nan,
        binarize=binarize,
        int_rounding=int_rounding,
    )

    return img, dat


def find_gzip(infile, raise_error=False, return_infile=False):
    """Find the existing file, gzipped or gunzipped.

    Return the infile if it exists, otherwise return the gzip-toggled
    version of the infile if it exists, otherwise return None or raise
    a FileNotFoundError.

    Parameters
    ----------
    infile : str
        The input file string.
    raise_error : bool
        If true, a FileNotFoundError is raised if the outfile does not
        exist.
    return_infile : bool
        If true, the infile is returned if the outfile does not exist.
        Otherwise None is returned if the outfile does not exist. This
        argument is ignored if raise_error is true.
    """
    if op.isfile(infile):
        outfile = infile
        return outfile
    elif op.isfile(toggle_gzip(infile)):
        outfile = toggle_gzip(infile)
        return outfile
    else:
        if raise_error:
            raise FileNotFoundError(
                "File not found: {}[.gz]".format(infile.replace(".gz", ""))
            )
        elif return_infile:
            return infile
        else:
            return None


def toggle_gzip(infile):
    """Return the gzip-toggled filepath.

    Parameters
    ----------
    infile : str
        The input file string.

    Returns
    -------
    outfile : str
        The output file string, which is the input file string minus
        the ".gz" extension if it exists in infile, or the input file
        string plus the ".gz" extension if it does not exist in infile.
    """
    if infile.endswith(".gz"):
        outfile = infile[:-3]
    else:
        outfile = infile + ".gz"
    return outfile


def _format_array(
    dat,
    dtype=np.float32,
    squeeze=True,
    flatten=False,
    conv_nan=0,
    binarize=False,
    int_rounding="nearest",
):
    """Format an array.

    Formatting options:
    - Flattening
    - NaN handling
    - Data type conversion

    Parameters
    ----------
    dtype : data-type
        Determines the data type returned.
    flatten : bool
        Return `dat` as a flattened copy of the input array.
    conv_nan : bool, number, or NoneType object
        Convert NaNs to `conv_nan`. No conversion is applied if
        `conv_nan` is np.nan, None, or False.
    binarize : bool
        If true, `dat` values > 0 become 1 and all other values are 0.
        `dat` type is recast to np.uint8.
    int_rounding : str
        Determines how the data array is recast if `binarize` is false
        and `dtype` is an integer.
        `nearest` : round to the nearest integer
        `floor` : round down
        `ceil` : round up

    Returns
    -------
    dat : ndarray or ndarray subclass
    """
    # Flatten the array.
    if flatten:
        dat = dat.ravel()

    # Squeeze the array.
    elif squeeze:
        dat = np.squeeze(dat)

    # Convert NaNs.
    if not np.any((conv_nan is None, conv_nan is False, conv_nan is np.nan)):
        dat[np.invert(np.isfinite(dat))] = conv_nan

    # Recast the data type.
    if binarize or (dtype is bool):
        idx = dat > 0
        dat[idx] = 1
        dat[~idx] = 0
        if dtype is bool:
            dat = dat.astype(bool)
        else:
            dat = dat.astype(np.uint8)
    elif "int" in str(dtype):
        if int_rounding == "nearest":
            dat = np.rint(dat)
        elif int_rounding == "floor":
            dat = np.floor(dat)
        elif int_rounding == "ceil":
            dat = np.ceil(dat)
        else:
            raise ValueError("int_rounding='{}' not valid".format(int_rounding))
        dat = dat.astype(dtype)
    else:
        dat = dat.astype(dtype)

    return dat


def splitt(string, delimiters):
    """Split a string into a list of substrings by 1+ delimiters."""
    if isinstance(delimiters, str):
        delimiters = [delimiters]

    pattern = "[" + re.escape("".join(delimiters)) + "]"
    return [substring for substring in re.split(pattern, string) if len(substring) > 0]


def _parse_args():
    """Parse and return command line arguments."""
    parser = argparse.ArgumentParser(
        description="""
$ extract_rois.py -i [pet1.nii pet2.nii ...] -m [mask1.nii mask2.nii ...]
OR
$ extract_rois.py -i [pet1.nii pet2.nii ...] -a [aparc+aseg1.nii aparc+aseg2.nii ...] -f [roi_file.csv]

But see options for further customization.

Extract image mean and ROI volume info from 1+ nifti images (-i|--images) and 1+ regions
of interest, which can be provided in the form of binary masks (-m|--masks) or
aparc+aseg style parcellation files (-a|--aparcs).

For parcellations, a 2-column CSV file (-f|--roi_file) is used that contains ROI names
and corresponding integer labels (or, for aggregate ROIs, a list of semicolon-separated
integer labels). It is possible to extract values from a subset of ROIs by specifying
them by name (-r|--aparc_rois). Otherwise, all ROIs in the CSV file are extracted by
default. --list_rois can be used to print the ROI names and labels from the CSV file
without doing any actual extractions.

Output is printed to the console by default, but can also be saved to a CSV file
(-o|--outputf). It is also possible to suppress printing the output dataframe to the
console (-q|--quiet).
        """,
        formatter_class=TextFormatter,
        exit_on_error=False,
    )
    parser.add_argument(
        "-i",
        "--images",
        type=str,
        nargs="+",
        help="Paths to 1+ images to extract values from",
    )
    parser.add_argument(
        "-m",
        "--masks",
        type=str,
        nargs="*",
        help="Paths to 1+ masks to apply over each image",
    )
    parser.add_argument(
        "-a",
        "--aparcs",
        type=str,
        nargs="*",
        help=(
            "Paths to 1+ parcellation files with ROI labels. Should match\n"
            + "the number of images, as images[i] is paired with aparcs[i]\n"
            + "for i = 1...len(images)"
        ),
    )
    parser.add_argument(
        "-r",
        "--aparc_rois",
        type=str,
        nargs="*",
        help="Names of ROIs to extract values from",
    )
    parser.add_argument(
        "-f",
        "--roi_file",
        type=str,
        help=(
            "Path to the 2-column CSV file with ROI names and int or\n"
            "semicolon-separated labels"
        ),
    )
    parser.add_argument(
        "-l",
        "--list_rois",
        action="store_true",
        help="List ROI names and labels from the ROIs CSV file",
    )
    parser.add_argument(
        "-o",
        "--outputf",
        type=str,
        help="Output CSV filepath. If not specified, output is printed but not saved",
    )
    parser.add_argument(
        "-s",
        "--shape",
        type=str,
        default="long",
        choices=["long", "wide", "l", "w"],
        help=(
            "Shape of the output dataframe. '-o wide' pivots the 'roi' column into\n"
            "multiple columns (one for each region)\n"
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Don't print the output dataframe to the terminal",
    )
    # Print help if no arguments are given.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args()
        return args


if __name__ == "__main__":
    timer = Timer()

    # Get command line arguments.
    args = _parse_args()

    # Print ROI names and labels in a nicely-formatted table.
    if args.list_rois:
        all_rois = pd.read_csv(args.roi_file)
        all_rois["n_labels"] = all_rois.iloc[:, 1].apply(lambda x: len(x.split(";")))
        all_rois.iloc[:, 1] = all_rois.iloc[:, 1].apply(_fmt_long_str)
        print(all_rois.to_markdown(index=False, tablefmt="rst"))
        print(args.roi_file, end="\n" * 2)
        sys.exit(0)

    # Load the ROI dictionary
    all_rois = load_rois(args.roi_file)

    # Check that at least one of masks or aparcs is specified
    if args.masks is None and args.aparcs is None:
        print(
            "ERROR: At least one of --masks (-m) or --aparcs (-a) must be specified\n"
        )
        sys.exit(1)

    # Check that the number of images and parcellations match
    if args.aparcs is not None:
        if (len(args.aparcs) > 1) and (len(args.images) != len(args.aparcs)):
            print(
                "ERROR: Number of images and parcellation files must match,\n"
                + "or there must be only one parcellation file specified.\n"
                + "Found {} images and {} parcellation files".format(
                    len(args.images), len(args.aparcs)
                )
            )
            sys.exit(1)

    # Extract ROI values from masks
    output = []
    if args.masks is not None:
        for img in args.images:
            _output = nops.roi_desc(dat=img, rois=args.masks)
            _output = _output.reset_index()
            _output.insert(0, "image_file", img)
            _output.insert(1, "roi_file", args.masks)
            output.append(_output)

    # Extract ROI values from parcellations
    if args.aparcs is not None:
        # Get ROIs
        keep_rois = {}
        if args.aparc_rois is None:
            keep_rois = all_rois
        else:
            for roi in args.aparc_rois:
                if roi not in all_rois.keys():
                    print(f"WARNING: {roi} missing from {args.roi_file}")
                else:
                    keep_rois[roi] = all_rois[roi]
        # Broadcast inputs if needed
        if (len(args.images) > 1) and (len(args.aparcs) == 1):
            args.aparcs = args.aparcs * len(args.images)
        # Extract ROI values
        for img, aparc in zip(args.images, args.aparcs):
            _output = nops.roi_desc(dat=img, rois=aparc, subrois=keep_rois)
            _output = _output.reset_index()
            _output.insert(0, "image_file", img)
            _output.insert(1, "roi_file", aparc)
            output.append(_output)

    output = pd.concat(output).reset_index(drop=True)
    output = output.rename(columns={"voxels": "voxel_count"})

    # Pivot the output dataframe
    if args.shape in ["wide", "w"]:
        output = output.pivot(
            index=["image_file", "roi_file"],
            columns="roi",
            values=["mean", "voxel_count"],
        )
        output.columns = ["_".join(col[::-1]).strip() for col in output.columns.values]
        output = output.reset_index()

    # Save output.
    if args.outputf is not None:
        output.to_csv(args.outputf, index=False)
        print(f"\nSaved output to {op.abspath(args.outputf)}")

    # Print output.
    if not args.quiet:
        output["image_file"] = output["image_file"].apply(op.basename)
        output["roi_file"] = output["roi_file"].apply(op.basename)
        for col in output.columns:
            if "mean" in col:
                output[col] = output[col].astype(float)
            elif "voxel_count" in col:
                output[col] = output[col].astype(float)
        output.columns = output.columns.str.replace("_", "\n")
        print(
            output.to_markdown(
                index=False,
                tablefmt="rst",
                floatfmt=".4f",
                intfmt=",",
            )
        )

    print(timer)
    sys.exit(0)
