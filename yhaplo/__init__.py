"""Yhaplo."""

import warnings
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    warnings.warn(f"{__package__} version unavailable. Is it installed?", stacklevel=1)
    __version__ = "Unknown"
