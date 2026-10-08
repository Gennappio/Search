"""One logger, to the console and to the run's log file."""

from __future__ import annotations

import logging
from pathlib import Path

FORMAT = "[%(asctime)s] %(message)s"


def get_logger(log_file: Path | None = None, level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger("search")
    logger.setLevel(level)
    logger.propagate = False
    if not any(
        isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
        for h in logger.handlers
    ):
        stream = logging.StreamHandler()
        stream.setFormatter(logging.Formatter(FORMAT, datefmt="%Y-%m-%d %H:%M:%S"))
        logger.addHandler(stream)
    if log_file is not None:
        wanted = str(Path(log_file))
        if not any(
            isinstance(h, logging.FileHandler) and h.baseFilename == wanted for h in logger.handlers
        ):
            handler = logging.FileHandler(log_file)
            handler.setFormatter(logging.Formatter(FORMAT, datefmt="%Y-%m-%d %H:%M:%S"))
            logger.addHandler(handler)
    return logger


def close_file_handlers() -> None:
    logger = logging.getLogger("search")
    for handler in list(logger.handlers):
        if isinstance(handler, logging.FileHandler):
            handler.close()
            logger.removeHandler(handler)
