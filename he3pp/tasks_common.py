from typing import Any

import ROOT

from .settings import RuntimeConfig


def default_tpc_model_name(runtime_config: RuntimeConfig) -> str:
    return "ExpGaus" if "ExpGaus" in runtime_config.tpc_function_names else str(runtime_config.tpc_function_names[0])


def weighted_eff_name(base: str) -> str:
    # Weighted efficiencies intentionally use "Weff*" naming.
    return f"W{base}"


def collect_rresult_ptrs(obj: Any) -> list[Any]:
    out: list[Any] = []
    if isinstance(obj, dict):
        for value in obj.values():
            out.extend(collect_rresult_ptrs(value))
        return out
    if isinstance(obj, (list, tuple)):
        for value in obj:
            out.extend(collect_rresult_ptrs(value))
        return out
    if hasattr(obj, "GetValue") and hasattr(obj, "GetPtr"):
        out.append(obj)
    return out


def run_graphs(actions: list[Any]) -> None:
    if not actions:
        return
    run_graphs_impl = getattr(getattr(ROOT, "RDF", None), "RunGraphs", None)
    if run_graphs_impl:
        run_graphs_impl(actions)
        return
    # Fallback for ROOT builds without RunGraphs.
    for action in actions:
        action.GetValue()


def find_object_by_class(root_dir: Any, class_name: str) -> Any | None:
    """Depth-first search for the first object with the requested ROOT class."""
    if not root_dir or not hasattr(root_dir, "GetListOfKeys"):
        return None
    for key in root_dir.GetListOfKeys():
        obj = key.ReadObj()
        if obj and hasattr(obj, "ClassName") and obj.ClassName() == class_name:
            return obj
        if obj and hasattr(obj, "InheritsFrom") and obj.InheritsFrom("TDirectory"):
            found = find_object_by_class(obj, class_name)
            if found:
                return found
    return None
