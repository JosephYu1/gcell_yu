from dataclasses import dataclass, field
from typing import TypeVar

import hydra
from hydra.core.global_hydra import GlobalHydra
from hydra_zen import MISSING
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


@dataclass
class AnnotationConfig:
    annotation_dir: str = field(default="")


@dataclass
class CellTypeConfig:
    data_dir: str = field(default="")
    interpret_dir: str = field(default="")
    assets_dir: str = field(default="")
    input: bool = field(default=False)
    jacob: bool = field(default=False)
    embed: bool = field(default=False)
    num_cls: int = field(default=2)
    motif_dir: str = field(default="")


@dataclass
class DemoConfig:
    celltype: CellTypeConfig = field(default_factory=CellTypeConfig)
    s3_uri: str = field(default="")
    s3_file_sys = MISSING
