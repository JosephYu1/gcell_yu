import os
from pathlib import Path

import pooch

# Default settings
DEFAULT_SETTINGS = {
    "annotation_dir": str(Path.home() / ".gcell_data" / "annotations"),
    "genome_dir": str(Path.home() / ".gcell_data" / "genomes"),
    "cache_dir": str(Path.home() / ".gcell_data" / "cache"),
}

# Current settings (can be modified at runtime)
SETTINGS = DEFAULT_SETTINGS.copy()


def setup_dirs():
    """Create necessary directories if they don't exist"""
    for dir_path in SETTINGS.values():
        Path(dir_path).mkdir(parents=True, exist_ok=True)


def get_setting(key: str) -> str:
    """Get a setting value, checking environment variables first"""
    env_key = f"GCELL_{key.upper()}"
    return os.environ.get(env_key, SETTINGS[key])


def update_settings(settings_dict: dict) -> None:
    """Update settings with new values"""
    SETTINGS.update(settings_dict)
    setup_dirs()


# Create a pooch instance for managing downloads
POOCH = pooch.create(
    path=get_setting("cache_dir"),
    base_url="https://example.com/",  # Replace with actual base URL
    registry={},  # Will be populated by individual modules
    retry_if_failed=3,
)

# Initialize directories
setup_dirs()


def download_with_pooch(fname, url, target_dir="cache_dir"):
    POOCH.path = Path(get_setting(target_dir))
    POOCH.registry[fname] = None
    POOCH.urls[fname] = url
    return POOCH.fetch(fname)
