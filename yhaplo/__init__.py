"""Yhaplo."""

import warnings
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__package__)
except PackageNotFoundError:
    warnings.warn(
        f"{__package__} version unavailable. Is it installed?",
        Warning,
        stacklevel=1,
    )
    __version__ = "Unknown"
