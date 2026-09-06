import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error

def RMSE(
    y_data1: pd.Series,
    y_data2: pd.Series,
) -> float:
    y_data1_align, y_data2_align = y_data1.align(y_data2, join = "inner")
    RMSE = np.sqrt(mean_squared_error(y_data1_align, y_data2_align))
    return RMSE

