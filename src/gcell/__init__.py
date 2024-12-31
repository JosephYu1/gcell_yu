from importlib.util import find_spec

from . import _settings

GRADIO_AVAILABLE = find_spec("gradio") is not None

__version__ = "0.1.0"

# Expose settings functions
get_setting = _settings.get_setting
update_settings = _settings.update_settings
