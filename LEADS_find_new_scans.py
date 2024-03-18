#!/usr/bin/env python

import sys
import os
import os.path as op

PATHS = {
    "proj": "/mnt/coredata/processing/leads",
}
PATHS["data"] = op.join(PATHS["proj"], "data")
PATHS["freesurfer"] = op.join(PATHS["data"], "freesurfer")
PATHS["metadata"] = op.join(PATHS["proj"], "metadata")
PATHS["processed"] = op.join(PATHS["data"], "processed")
PATHS["raw"] = op.join(PATHS["data"], "raw")
