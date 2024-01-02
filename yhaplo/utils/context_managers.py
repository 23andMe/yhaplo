"""Utility context managers."""

import logging
import time
from collections.abc import Iterator
from contextlib import contextmanager

logger = logging.getLogger(__name__)


@contextmanager
def logging_disabled(level: int = logging.INFO) -> Iterator[None]:
    """Disable log messages at specified level and below.

    Parameters
    ----------
    level : int
        Logging level to disable.

    Returns
    -------
    Context manager for disabling logging.

    """
    previous_disable_level = logging.root.manager.disable  # noqa
    logging.disable(level)

    try:
        yield
    finally:
        logging.disable(previous_disable_level)


@contextmanager
def timer(label: str = "Time") -> Iterator[None]:
    """Time an operation.

    Parameters
    ----------
    label : str
        Label for what is being timed, for log message.

    Returns
    -------
    Context manager for timing a block of code.

    """
    start_time = time.time()
    yield
    elapsed_time = time.time() - start_time
    logger.info(f"{label}: {elapsed_time:.3f} seconds\n")
