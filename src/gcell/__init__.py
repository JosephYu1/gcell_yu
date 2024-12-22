from __future__ import annotations

from importlib.util import find_spec

GRADIO_AVAILABLE = find_spec("gradio") is not None

__version__ = "0.1.0"
