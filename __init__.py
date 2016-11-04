# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import pkg_resources
import warnings

try:
    __version__ = pkg_resources.get_distribution(__name__).version
except pkg_resources.DistributionNotFound as e:  # pragma: no cover
    warnings.warn("Version not available for package %s (is it installed?)" % __package__,
                  Warning)
    __version__ = None
