import pandas as pd

class ProcessUnit:
    def __init__(
        self,
        table: pd.DataFrame,
        sample_name_column: str = "sample_name",
        file_path_column: str = "path",
    ):
        self._clsname = self.__class__.__name__
        self.table = table
        self.files = self.table[file_path_column].to_list()
        self.sample_name_column = sample_name_column
        self.file_path_column = file_path_column

        self.sample = self._validate_single_sample()

    @property
    def sample_id(self):
        return self.sample

    def _validate_single_sample(self) -> str:
        """Ensure unit contains exactly one sample"""
        unique_samples = self.table[self.sample_name_column].unique()
        if len(unique_samples) != 1:
            raise ValueError(
                f"ProcessUnit must contain exactly one sample. Found: {unique_samples}"
            )
        return unique_samples[0]

