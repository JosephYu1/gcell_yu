from __future__ import annotations

from dataclasses import field
from typing import TypeVar

import hydra
from hydra.core.global_hydra import GlobalHydra
from omegaconf import OmegaConf

T = TypeVar("T")


def default_cfg():
    return field(default_factory=T)


def load_config(config_name, config_path="./"):
    # Initialize Hydra to load the configuration
    GlobalHydra.instance().clear()
    hydra.initialize(config_path=config_path, version_base="1.3")
    cfg = hydra.compose(config_name=config_name)
    return cfg


def pretty_print_config(cfg):
    print(OmegaConf.to_yaml(cfg))
