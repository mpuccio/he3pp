#!/usr/bin/env python3
import argparse
from dataclasses import dataclass
import os

ROOT = None


@dataclass
class Result:
    ref_keys: int
    cand_keys: int
    common: int
    missing: int
    extra: int
    content_diffs: int
    error_diffs: int
    max_content_diff: float
    max_error_diff: float


def collect(tdir, prefix=""):
    out = {}
    for key in tdir.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()
        path = f"{prefix}/{name}" if prefix else name
        if obj.InheritsFrom("TDirectory"):
            out.update(collect(obj, path))
        else:
            # Keep only histograms and detach them from ROOT files to avoid
            # use-after-free issues in PyROOT.
            if obj.ClassName().startswith("TH"):
                clone = obj.Clone(f"{obj.GetName()}__cmp")
                clone.SetDirectory(0)
                out[path] = clone
    return out


def normalize_key(path: str, normalize_weff: bool) -> str:
    if not normalize_weff:
        return path
    # Canonicalize both naming schemes to Weff* so comparisons are symmetric.
    return (
        path.replace("/effWTPCA", "/WeffTPCA")
        .replace("/effWTPCM", "/WeffTPCM")
        .replace("/effWTOFA", "/WeffTOFA")
        .replace("/effWTOFM", "/WeffTOFM")
    )


def compare(
    ref_path: str,
    cand_path: str,
    content_tol: float,
    error_tol: float,
    normalize_weff: bool,
    ignore_errors: bool,
    allow_extra_keys: list[str],
    allow_extra_prefixes: list[str],
) -> Result:
    global ROOT
    if ROOT is None:
        import ROOT as _ROOT
        ROOT = _ROOT

    ref = collect(ROOT.TFile(os.path.expandvars(os.path.expanduser(ref_path))))
    cand_raw = collect(ROOT.TFile(os.path.expandvars(os.path.expanduser(cand_path))))
    cand = {normalize_key(k, normalize_weff): v for k, v in cand_raw.items()}

    common = sorted(set(ref) & set(cand))
    missing = sorted(set(ref) - set(cand))
    extra_raw = sorted(set(cand) - set(ref))
    extra = [
        k
        for k in extra_raw
        if k not in allow_extra_keys and not any(k.startswith(prefix) for prefix in allow_extra_prefixes)
    ]

    content_diffs = 0
    error_diffs = 0
    max_content_diff = 0.0
    max_error_diff = 0.0

    for k in common:
        h1 = ref[k]
        h2 = cand[k]
        d = h1.GetDimension()
        nx = h1.GetNbinsX() + 2
        ny = h1.GetNbinsY() + 2 if d > 1 else 1
        nz = h1.GetNbinsZ() + 2 if d > 2 else 1
        local_max_content = 0.0
        local_max_error = 0.0
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    b = h1.GetBin(ix, iy, iz) if d > 1 else h1.GetBin(ix)
                    local_max_content = max(local_max_content, abs(h1.GetBinContent(b) - h2.GetBinContent(b)))
                    local_max_error = max(local_max_error, abs(h1.GetBinError(b) - h2.GetBinError(b)))
        max_content_diff = max(max_content_diff, local_max_content)
        max_error_diff = max(max_error_diff, local_max_error)
        if local_max_content > content_tol:
            content_diffs += 1
        if (not ignore_errors) and local_max_error > error_tol:
            error_diffs += 1

    return Result(
        ref_keys=len(ref),
        cand_keys=len(cand),
        common=len(common),
        missing=len(missing),
        extra=len(extra),
        content_diffs=content_diffs,
        error_diffs=error_diffs,
        max_content_diff=max_content_diff,
        max_error_diff=max_error_diff,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare ROOT outputs for smoke regression")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--content-tol", type=float, default=1e-9)
    parser.add_argument("--error-tol", type=float, default=1e-9)
    parser.add_argument("--normalize-weff", action="store_true", help="Treat Weff* and effW* naming as equivalent")
    parser.add_argument("--ignore-errors", action="store_true", help="Ignore histogram error differences")
    parser.add_argument(
        "--allow-extra-key",
        action="append",
        default=[],
        help="Candidate-only histogram key allowed by design. Repeat for multiple keys.",
    )
    parser.add_argument(
        "--allow-extra-prefix",
        action="append",
        default=[],
        help="Candidate-only histogram key prefix allowed by design. Repeat for multiple prefixes.",
    )
    args = parser.parse_args()

    res = compare(
        args.reference,
        args.candidate,
        args.content_tol,
        args.error_tol,
        args.normalize_weff,
        args.ignore_errors,
        args.allow_extra_key,
        args.allow_extra_prefix,
    )
    print(f"ref_keys={res.ref_keys} cand_keys={res.cand_keys} common={res.common} missing={res.missing} extra={res.extra}")
    print(f"content_diffs={res.content_diffs} max_content_diff={res.max_content_diff}")
    print(f"error_diffs={res.error_diffs} max_error_diff={res.max_error_diff}")

    ok = (res.missing == 0 and res.extra == 0 and res.content_diffs == 0 and (args.ignore_errors or res.error_diffs == 0))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
