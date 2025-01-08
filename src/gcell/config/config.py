from dataclasses import field
from typing import TypeVar

import hydra
from hydra.core.global_hydra import GlobalHydra
from omegaconf import OmegaConf

T = TypeVar("T")


def default_cfg():
    return field(default_factory=T)


def load_config(config_name, config_path="./"):
    """
    Load a configuration file using Hydra.

    Parameters
    ----------
    config_name : str
        The name of the configuration file to load.
    config_path : str, optional
        The path to the configuration file relative to where config.py is located.


    Note
    -----
    Hydra only supports relative paths to the parent of the caller, that is, this config.py file.

    Returns
    -------
    cfg : dict
        The configuration dictionary.
    """
    GlobalHydra.instance().clear()
    hydra.initialize(config_path=config_path, version_base="1.3")
    cfg = hydra.compose(config_name=config_name)
    return cfg


def pretty_print_config(cfg):
    """
    Pretty print the configuration dictionary.
    """
    print(OmegaConf.to_yaml(cfg))
