import sys
import os.path as op

# Get the absolute path to the top-level directory (where this __init__.py is located)
this_pkg = op.abspath(op.dirname(__file__))

# Add the top-level directory to sys.path
if this_pkg not in sys.path:
    sys.path.append(this_pkg)
