"""Call haplogroups.

For details, run:
    yhaplo --help

"""

import logging

from yhaplo.api.call_haplogroups import call_haplogroups
from yhaplo.api.command_line_args import get_command_line_args

root_logger = logging.getLogger()


def main() -> None:
    """Configure logging and call haplogroups."""

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    command_line_args = get_command_line_args()
    call_haplogroups(
        command_line_args=command_line_args,
        root_logger=root_logger,
    )


if __name__ == "__main__":
    call_haplogroups()
