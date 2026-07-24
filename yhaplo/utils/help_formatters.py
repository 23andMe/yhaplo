"""Argparse help formatters."""

import argparse


class CapitalizedHelpFormatter(argparse.RawTextHelpFormatter):
    """RawTextHelpFormatter with a capitalized usage prefix and section headings."""

    def add_usage(self, usage, actions, groups, prefix=None):
        """Capitalize the usage prefix."""

        if prefix is None:
            prefix = "Usage: "

        super().add_usage(usage, actions, groups, prefix)

    def start_section(self, heading):
        """Capitalize the first letter of a section heading."""

        if heading is not None:
            heading = heading[0].upper() + heading[1:]

        super().start_section(heading)
