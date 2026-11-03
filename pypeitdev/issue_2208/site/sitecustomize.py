"""
Applies pca_patch to any Python process started with this directory on
PYTHONPATH, when I2208_PCA_FIX or I2208_PCA_CAPTURE is set.  Issue #2208.
"""
import os
import sys

if os.environ.get('I2208_PCA_FIX') or os.environ.get('I2208_PCA_CAPTURE'):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import pca_patch
    pca_patch.install_from_env()
