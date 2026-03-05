import numpy as np
import pandas as pd

global ChromiumNextGEMSingleCell3primeHTv3_1__multiplet_rates
global ChromiumNextGEMSingleCell5primeHTv2__multiplet_rates

ChromiumNextGEMSingleCell3primeHTv3_1__multiplet_rates = pd.DataFrame(
    {
        "multiplet_rate": np.array([0.8, 1.6, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 8.0])
        / 100,
        "n_cells_loaded": np.array(
            [3130, 6320, 9550, 12800, 16100, 19500, 22900, 26300, 29800, 33300]
        ),
        "n_cells_recovered": np.array(
            [2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000]
        ),
    }
)

ChromiumNextGEMSingleCell5primeHTv2__multiplet_rates = dict(ChromiumNextGEMSingleCell3primeHTv3_1__multiplet_rates)
# %%

