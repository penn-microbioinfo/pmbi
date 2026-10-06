import logging
from typing import Union

DEFAULT_FMT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

# Logger for this module
logger = logging.getLogger(__name__)
h = logging.StreamHandler()
h.setFormatter(logging.Formatter(fmt=DEFAULT_FMT))
logger.addHandler(h)
logger.setLevel(level=logging.INFO)

def list_loggers():
    return [logging.getLogger(name) for name in logging.root.manager.loggerDict]

def handler_with_formatter(handler: logging.Handler, formatter: logging.Formatter):
        handler.setFormatter(formatter)
        return handler

def logger_exists(name):
    if name in [l.name for l in list_loggers()]:
        return True
    else:
        return False

class LoggerConstructor:
    def __init__(
        self,
        name,
        handlers: list[logging.Handler],
        level="INFO",
    ):

        if logger_exists(name):
            preexisting = logging.getLogger(name)
            logger.warning(f"Overwriting preexisting logger with name: {preexisting.name}")
            if preexisting.hasHandlers():
                logger.warning(f"Clearing handlers for logger: {preexisting.name}")
                preexisting.handlers.clear()
            

        self.logger = None
        self.name = name
        self.level = getattr(logging, level)

        # Handlers
        self.handlers = handlers

    def get(self) -> logging.Logger:
        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(self.level)

        for handler in self.handlers:
            self.logger.addHandler(handler)

        return self.logger


########################################
# Some preset logger configurations

def stream_file_logger(name, filename, level="INFO"):
    constructor = LoggerConstructor(
        name=name,
        handlers=[
            handler_with_formatter(handler=logging.StreamHandler(), formatter=logging.Formatter(fmt=DEFAULT_FMT)),
            handler_with_formatter(handler=logging.FileHandler(filename=filename), formatter=logging.Formatter(fmt=DEFAULT_FMT))
        ],
        level=level)

    return constructor.get()

def streamLogger(name, level="INFO"):
    constructor = LoggerConstructor(
        name=name,
        handlers=[
            handler_with_formatter(handler=logging.StreamHandler(), formatter=logging.Formatter(fmt=DEFAULT_FMT)),
        ],
        level=level)

    return constructor.get()


