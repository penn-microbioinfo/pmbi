import logging

def list_loggers():
    return [logging.getLogger(name) for name in logging.root.manager.loggerDict]

def logger_exists(name):
    if name in [l.name for l in list_loggers()]:
        return True
    else:
        return False

def streamLogger(name, level="INFO"):
    if logger_exists(name):
        return logging.getLogger(name)

    else:
        logger = logging.getLogger(name)
        logger.setLevel(getattr(logging, level))
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")) 
        logger.addHandler(handler)

        return logger

