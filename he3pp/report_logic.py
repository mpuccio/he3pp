from typing import Any


def status_from_metrics(available: bool, metrics: dict[str, Any], alpha: float) -> tuple[str, str]:
    if not available:
        return "MISSING", "missing"
    if not metrics:
        return "UNK", "unknown"
    return ("OK", "ok") if float(metrics.get("p_value", -1.0)) >= float(alpha) else ("KO", "ko")
