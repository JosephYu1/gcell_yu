import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d


def min_max(x):
    """normalize an array to 0-1"""
    if isinstance(x, np.ndarray):
        return (x - x.min()) / (x.max() - x.min())
    elif isinstance(x, list):
        return [(a - min(x)) / (max(x) - min(x)) for a in x]
    elif isinstance(x, pd.DataFrame):
        arr_copy = x.copy()
        arr_copy["plddt"] = (arr_copy["plddt"] - arr_copy["plddt"].min()) / (
            arr_copy["plddt"].max() - arr_copy["plddt"].min()
        )
        return arr_copy


def normalize(x):
    new_x = x - x.min()
    return new_x / new_x.max()


def smooth(x, window_size=10):
    result = gaussian_filter1d(x, sigma=window_size / 2)
    return normalize(result)


def get_3d_avg(x, pairwise_interaction):
    x = x * pairwise_interaction
    x = x.sum(1) / ((x > 0).sum(1) + 0.01)
    return x / x.max()


def square_grad(f):
    return normalize(np.gradient(f) ** 2)


def extract_wt_from_mut(str):
    return str[0:1]


def extract_alt_from_mut(str):
    return str[-1:]


def extract_pos_from_mut(str):
    return int(str[1:-1])
