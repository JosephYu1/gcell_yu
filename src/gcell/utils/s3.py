from __future__ import annotations

from pathlib import Path

import numpy as np
import s3fs
import zarr
from scipy.sparse import load_npz


def open_file_with_s3(file_path, mode="r", s3_file_sys=None):
    if s3_file_sys:
        return s3_file_sys.open(file_path, mode=mode)
    else:
        return Path(file_path).open(mode=mode)


def path_exists_with_s3(file_path, s3_file_sys=None):
    if s3_file_sys:
        return s3_file_sys.exists(file_path)
    else:
        return Path(file_path).exists()


def glob_with_s3(file_path, s3_file_sys=None):
    if s3_file_sys:
        return s3_file_sys.glob(file_path)
    else:
        return list(Path().glob(file_path))


def load_zarr_with_s3(file_path, mode="r", s3_file_sys=None):
    if s3_file_sys:
        return zarr.open(s3fs.S3Map(file_path, s3=s3_file_sys), mode=mode)
    else:
        return zarr.open(str(Path(file_path)), mode=mode)


def load_np_with_s3(file_path, s3_file_sys=None):
    if s3_file_sys:
        return np.load(s3_file_sys.open(file_path))
    else:
        return np.load(Path(file_path))


def load_npz_with_s3(file_path, s3_file_sys=None):
    if s3_file_sys:
        return load_npz(s3_file_sys.open(file_path))
    else:
        return load_npz(Path(file_path))
