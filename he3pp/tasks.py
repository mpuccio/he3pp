"""Compatibility facade for task entry points.

Task implementations are split across focused modules to keep this interface stable.
"""

from .tasks_checkpoint import checkpoint
from .tasks_data import analyse_data
from .tasks_mc import analyse_mc
from .tasks_signal import signal
from .tasks_systematics import systematics


__all__ = [
    "analyse_data",
    "analyse_mc",
    "signal",
    "systematics",
    "checkpoint",
]
