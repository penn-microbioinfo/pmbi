from pmbi.config._config import Config
from pmbi.config.item import Item, RegexItem, PathItem

class STARsolo_config(Config):
    def __init__(self):
        super().__init__() 
        self._items = {
            RegexItem("filename_patterns.sample", capture_groups=1),
            RegexItem("filename_patterns.modality", capture_groups=1),
            RegexItem("filename_patterns.include", capture_groups=1, optional=True),
            RegexItem("filename_patterns.exclude", capture_groups=1, optional=True),
            Item("run.wd"),
            Item("run.flavor"),
            Item("run.barcode_whitelist"),
            Item("run.n_jobs", type_=int, optional=True),
        }

