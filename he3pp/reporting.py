from __future__ import annotations

from dataclasses import dataclass, field
import html
import json
import os
from pathlib import Path
from typing import Any

import ROOT

from . import settings as s
from .root_io import expand


@dataclass
class ReportItem:
    section: str
    title: str
    source_file: str
    object_path: str
    output_png: str
    note: str = ""
    bin_index: int | None = None


@dataclass
class RenderedItem:
    item: ReportItem
    available: bool
    status: str
    status_class: str
    note: str = ""
    metrics: dict[str, float] = field(default_factory=dict)


SECTION_DEFS = {
    "signal_tof": ("TOF Signal Extraction", "Per-bin TOF fit snapshots."),
    "signal_tpc": ("TPC-only Signal Extraction", "Per-bin TPC-only fit snapshots."),
    "raw_spectrum": ("Raw pT Spectrum", "Raw TOF spectrum before efficiency correction."),
    "efficiency": ("Efficiency", "TPC and TOF efficiencies with binomial uncertainty."),
    "pt_resolution": ("pT Resolution", "Comparison of reconstructed and generated transverse momentum."),
    "corrected_spectrum": ("Corrected Spectrum", "Final corrected spectra from systematics output."),
}
DEFAULT_SECTIONS = ["signal_tof", "signal_tpc", "raw_spectrum", "efficiency", "pt_resolution", "corrected_spectrum"]


def _mkdir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def _get(root_file: ROOT.TFile, path: str) -> Any | None:
    obj = root_file.Get(path)
    return obj if obj else None


def _species_letter(species: str) -> str:
    return "A" if species.startswith("anti") else "M"


def _preflight_inputs(signal_file_path: str, mc_file_path: str, systematics_file_path: str) -> dict[str, str]:
    required = {
        "signal": expand(signal_file_path),
        "mc": expand(mc_file_path),
        "systematics": expand(systematics_file_path),
    }
    missing = [f"{kind}: {path}" for kind, path in required.items() if not os.path.isfile(path)]
    if missing:
        raise FileNotFoundError(
            "Report preflight failed. Missing required ROOT input file(s):\n"
            + "\n".join(missing)
            + "\nRun the pipeline stages that produce these files first (e.g. task='full_chain')."
        )
    return required


def _draw_to_png(obj: Any, png_path: str) -> bool:
    if obj is None:
        return False
    _mkdir(str(Path(png_path).parent))
    c = ROOT.TCanvas("c_report", "c_report", 950, 700)
    cname = obj.ClassName()
    if cname.startswith("TH2"):
        c.SetRightMargin(0.14)
        obj.Draw("colz")
    elif cname.startswith("TH"):
        obj.Draw("hist")
    elif cname.startswith("TGraph"):
        obj.Draw("ALP")
    elif cname == "RooPlot":
        obj.Draw()
    else:
        try:
            obj.Draw()
        except Exception:
            c.Close()
            return False
    c.SaveAs(png_path)
    c.Close()
    return True


def _signal_tof_bin_paths(signal_file: ROOT.TFile, species: str) -> list[tuple[int, str]]:
    base = signal_file.Get(f"nuclei/{species}/GausExp/C_0")
    if not base:
        return []
    keys: list[tuple[int, str]] = []
    for k in base.GetListOfKeys():
        n = k.GetName()
        if n.startswith("d0_") and not n.endswith("_sideband"):
            idx = n.split("_", 1)[1]
            try:
                order = int(idx)
            except ValueError:
                order = 10_000
            keys.append((order, f"nuclei/{species}/GausExp/C_0/{n}"))
    keys.sort(key=lambda x: x[0])
    return keys


def _signal_tpc_bin_paths(signal_file: ROOT.TFile, species: str, model: str) -> list[tuple[int, str]]:
    base = signal_file.Get(f"nuclei/{species}/TPConly")
    if not base:
        return []
    out: list[tuple[int, str]] = []
    suffix = f"_{model}"
    for k in base.GetListOfKeys():
        n = k.GetName()
        if not n.startswith("TPC_d_0_") or not n.endswith(suffix):
            continue
        parts = n.split("_")
        if len(parts) < 5:
            continue
        try:
            bidx = int(parts[3])
        except ValueError:
            continue
        out.append((bidx, f"nuclei/{species}/TPConly/{n}"))
    out.sort(key=lambda x: x[0])
    return out


def _normalize_sections(sections: list[str] | None) -> list[str]:
    if not sections:
        return list(DEFAULT_SECTIONS)
    out: list[str] = []
    for sec in sections:
        sname = str(sec).strip().lower()
        if sname in SECTION_DEFS and sname not in out:
            out.append(sname)
    return out or list(DEFAULT_SECTIONS)


def _chi2_metrics(obj: Any, n_fit_parameters: int, fit_tail: str) -> dict[str, float]:
    if not obj or obj.ClassName() != "RooPlot":
        return {}
    try:
        chi2_ndf = float(obj.chiSquare("model", "data", int(n_fit_parameters)))
    except Exception:
        return {}
    if chi2_ndf <= 0:
        return {}
    try:
        data_hist = obj.getHist("data")
        n_points = int(data_hist.GetN()) if data_hist else 0
    except Exception:
        n_points = 0
    if n_points <= 0:
        return {}
    dof = max(1, n_points - int(n_fit_parameters))
    chi2 = max(0.0, chi2_ndf * dof)
    p_right = float(ROOT.TMath.Prob(chi2, dof))
    if fit_tail == "two":
        p_value = min(1.0, 2.0 * min(p_right, 1.0 - p_right))
    else:
        p_value = p_right
    return {"chi2_ndf": chi2_ndf, "dof": float(dof), "p_value": p_value}


def _status_from_metrics(available: bool, metrics: dict[str, float], alpha: float) -> tuple[str, str]:
    if not available:
        return "MISSING", "missing"
    if not metrics:
        return "UNK", "unknown"
    return ("OK", "ok") if metrics["p_value"] >= alpha else ("KO", "ko")


def _render_card(row: RenderedItem) -> str:
    it = row.item
    safe_title = html.escape(it.title)
    safe_obj = html.escape(it.object_path)
    out = [
        "<article class='plot-card'>",
        f"<header><h3>{safe_title}</h3><span class='status {row.status_class}'>{row.status}</span></header>",
        f"<p class='path'><code>{safe_obj}</code></p>",
    ]
    if row.metrics:
        out.append(f"<p class='note'>p-value={row.metrics['p_value']:.3f} ({'single-tail' if row.metrics.get('tail_single', 1.0) else 'two-tail'})</p>")
    if row.note:
        out.append(f"<p class='note'>{html.escape(row.note)}</p>")
    if row.available:
        out.append(f"<img loading='lazy' src='{html.escape(it.output_png)}' alt='{safe_title}'>")
    out.append("</article>")
    return "".join(out)


def _render_section(title: str, subtitle: str, rows: list[RenderedItem], carousel_id: str | None = None) -> str:
    safe_title = html.escape(title)
    safe_subtitle = html.escape(subtitle)
    out = [
        "<details class='report-section' open>",
        f"<summary><span class='summary-title'>{safe_title}</span><span class='summary-sub'>{safe_subtitle}</span></summary>",
    ]
    if carousel_id:
        cid = html.escape(carousel_id)
        out.append(
            f"<div class='carousel' data-carousel='{cid}'>"
            f"<button class='carousel-btn prev' type='button' aria-label='Previous' data-target='{cid}'>&lsaquo;</button>"
            f"<div class='carousel-track' id='{cid}'>"
        )
        out.extend(_render_card(r) for r in rows)
        out.append(
            "</div>"
            f"<button class='carousel-btn next' type='button' aria-label='Next' data-target='{cid}'>&rsaquo;</button>"
            "</div>"
        )
    else:
        out.append("<div class='plot-grid'>")
        out.extend(_render_card(r) for r in rows)
        out.append("</div>")
    out.append("</details>")
    return "".join(out)


def _fmt_opt(value: float | None, fmt: str) -> str:
    if value is None:
        return "-"
    return format(value, fmt)


def _render_summary_table(
    tof_rows: list[RenderedItem],
    tpc_rows: list[RenderedItem],
    tof_yield: Any,
    tpc_yield: Any,
    tpc_eff: Any,
    tof_eff: Any,
) -> str:
    if not tof_rows and not tpc_rows:
        return ""
    p_tof = {r.item.bin_index: r.metrics.get("p_value") for r in tof_rows if r.item.bin_index is not None}
    p_tpc = {r.item.bin_index: r.metrics.get("p_value") for r in tpc_rows if r.item.bin_index is not None}
    all_bins = sorted(set(p_tof) | set(p_tpc))
    if tof_yield:
        all_bins = list(range(tof_yield.GetNbinsX()))
    axis = tof_yield.GetXaxis() if tof_yield else (tpc_yield.GetXaxis() if tpc_yield else None)

    out = [
        "<details class='report-section' open>",
        "<summary><span class='summary-title'>Fit And Efficiency Summary</span><span class='summary-sub'>TOF/TPC yields, efficiencies and fit p-values by pT bin.</span></summary>",
        "<div class='table-wrap'><table class='summary-table'><thead><tr>"
        "<th>Bin</th><th>pT Range</th><th>TOF Yield</th><th>TPC Yield</th><th>TPC Eff</th><th>TOF Eff</th><th>TOF p-value</th><th>TPC p-value</th>"
        "</tr></thead><tbody>",
    ]
    for b in all_bins:
        p_range = "-"
        if axis:
            lo = axis.GetBinLowEdge(b + 1)
            hi = axis.GetBinUpEdge(b + 1)
            p_range = f"{lo:.2f} - {hi:.2f}"
        tof_y = tof_yield.GetBinContent(b + 1) if tof_yield else None
        tpc_y = tpc_yield.GetBinContent(b + 1) if tpc_yield else None
        tpc_e = tpc_eff.GetBinContent(b + 1) if tpc_eff else None
        tof_e = tof_eff.GetBinContent(b + 1) if tof_eff else None
        out.append(
            "<tr>"
            f"<td>{b}</td>"
            f"<td>{html.escape(p_range)}</td>"
            f"<td>{_fmt_opt(tof_y, '.2f')}</td>"
            f"<td>{_fmt_opt(tpc_y, '.2f')}</td>"
            f"<td>{_fmt_opt(tpc_e, '.4f')}</td>"
            f"<td>{_fmt_opt(tof_e, '.4f')}</td>"
            f"<td>{_fmt_opt(p_tof.get(b), '.3f')}</td>"
            f"<td>{_fmt_opt(p_tpc.get(b), '.3f')}</td>"
            "</tr>"
        )
    out.append("</tbody></table></div></details>")
    return "".join(out)


def generate_report(
    report_dir: str,
    signal_file_path: str,
    mc_file_path: str,
    systematics_file_path: str,
    metadata_path: str | None = None,
    species: str = "antihe3",
    sections: list[str] | None = None,
    fit_n_parameters: int = 6,
    fit_alpha: float = 0.05,
    fit_tail: str = "single",
    tpc_signal_model: str = "ExpGaus",
) -> str:
    report_dir = expand(report_dir)
    assets_dir = str(Path(report_dir) / "assets")
    _mkdir(assets_dir)

    section_order = _normalize_sections(sections)
    fit_tail = "two" if str(fit_tail).lower() == "two" else "single"

    resolved = _preflight_inputs(signal_file_path, mc_file_path, systematics_file_path)
    sig_file = ROOT.TFile(resolved["signal"])
    mc_file = ROOT.TFile(resolved["mc"])
    syst_file = ROOT.TFile(resolved["systematics"])

    items: list[ReportItem] = []

    for bidx, p in _signal_tof_bin_paths(sig_file, species):
        items.append(
            ReportItem(
                section="signal_tof",
                title=f"TOF Signal Bin {bidx}",
                source_file=signal_file_path,
                object_path=p,
                output_png=f"assets/signal_tof_bin_{bidx:02d}.png",
                bin_index=bidx,
            )
        )
    for bidx, p in _signal_tpc_bin_paths(sig_file, species, tpc_signal_model):
        items.append(
            ReportItem(
                section="signal_tpc",
                title=f"TPC Signal Bin {bidx}",
                source_file=signal_file_path,
                object_path=p,
                output_png=f"assets/signal_tpc_bin_{bidx:02d}.png",
                bin_index=bidx,
                note=f"Model: {tpc_signal_model}",
            )
        )

    raw_name = "hRawCountsA0" if species == "antihe3" else "hRawCountsM0"
    items.append(
        ReportItem(
            section="raw_spectrum",
            title="Raw pT Spectrum",
            source_file=signal_file_path,
            object_path=f"nuclei/{species}/GausExp/{raw_name}",
            output_png="assets/raw_pt_spectrum.png",
        )
    )

    items.extend(
        [
            ReportItem(
                section="efficiency",
                title="TPC Efficiency",
                source_file=mc_file_path,
                object_path="nuclei/effTPCA" if species == "antihe3" else "nuclei/effTPCM",
                output_png="assets/eff_tpc.png",
                note="Uncertainty: binomial",
            ),
            ReportItem(
                section="efficiency",
                title="TOF Efficiency",
                source_file=mc_file_path,
                object_path="nuclei/effTOFA" if species == "antihe3" else "nuclei/effTOFM",
                output_png="assets/eff_tof.png",
                note="Uncertainty: binomial",
            ),
        ]
    )

    items.extend(
        [
            ReportItem(
                section="pt_resolution",
                title="pT(rec)-pT(sim) Uncorrected",
                source_file=mc_file_path,
                object_path="nuclei/hDeltaPtHe3",
                output_png="assets/delta_pt_uncorrected.png",
            ),
            ReportItem(
                section="pt_resolution",
                title="pT(rec)-pT(sim) Corrected",
                source_file=mc_file_path,
                object_path="nuclei/hDeltaPtCorrHe3",
                output_png="assets/delta_pt_corrected.png",
            ),
        ]
    )

    items.extend(
        [
            ReportItem(
                section="corrected_spectrum",
                title="Corrected Spectrum TPC",
                source_file=systematics_file_path,
                object_path="fStatTPCA" if species == "antihe3" else "fStatTPCM",
                output_png="assets/corrected_tpc.png",
            ),
            ReportItem(
                section="corrected_spectrum",
                title="Corrected Spectrum TOF",
                source_file=systematics_file_path,
                object_path="fStatTOFA" if species == "antihe3" else "fStatTOFM",
                output_png="assets/corrected_tof.png",
            ),
        ]
    )

    file_map = {
        signal_file_path: sig_file,
        mc_file_path: mc_file,
        systematics_file_path: syst_file,
    }

    rows: list[RenderedItem] = []
    for it in items:
        f = file_map[it.source_file]
        obj = _get(f, it.object_path)
        available = _draw_to_png(obj, str(Path(report_dir) / it.output_png))
        fit_section = it.section in ("signal_tof", "signal_tpc")
        metrics = _chi2_metrics(obj, fit_n_parameters, fit_tail) if fit_section else {}
        if metrics:
            metrics["tail_single"] = 1.0 if fit_tail == "single" else 0.0
        status, status_class = _status_from_metrics(available, metrics, fit_alpha) if fit_section else (("OK", "ok") if available else ("MISSING", "missing"))
        rows.append(RenderedItem(item=it, available=available, status=status, status_class=status_class, note=it.note, metrics=metrics))

    row_by_section = {sec: [r for r in rows if r.item.section == sec] for sec in SECTION_DEFS}

    letter = _species_letter(species)
    tof_yield = _get(sig_file, f"nuclei/{species}/GausExp/{raw_name}")
    tpc_yield = _get(sig_file, f"nuclei/{species}/TPConly/hTPConly{letter}0_{tpc_signal_model}")
    tpc_eff = _get(mc_file, "nuclei/effTPCA" if species == "antihe3" else "nuclei/effTPCM")
    tof_eff = _get(mc_file, "nuclei/effTOFA" if species == "antihe3" else "nuclei/effTOFM")

    meta = {}
    if metadata_path:
        meta_path = expand(metadata_path)
        if os.path.exists(meta_path):
            with open(meta_path, "r", encoding="utf-8") as f:
                meta = json.load(f)

    html_path = Path(report_dir) / "index.html"

    style = """
<style>
:root { --bg:#f6f3ec; --panel:#fffef8; --ink:#1f2a30; --muted:#5f6c74; --accent:#005f73; --ok:#127143; --ko:#9f2a2a; --missing:#7a7a7a; --unknown:#885f00; --border:#d9d2c4; }
* { box-sizing:border-box; }
body { margin:0; font-family:"Avenir Next","Segoe UI",sans-serif; color:var(--ink); background:radial-gradient(circle at 20% 15%, #e4efe9 0%, transparent 30%),radial-gradient(circle at 85% 10%, #f4dfc8 0%, transparent 30%),var(--bg); }
.wrap { max-width:1220px; margin:0 auto; padding:28px 20px 44px; }
h1 { margin:0 0 10px; font-size:2rem; letter-spacing:0.02em; }
.meta { color:var(--muted); margin-bottom:18px; }
details { background:var(--panel); border:1px solid var(--border); border-radius:12px; padding:12px 14px; margin-bottom:18px; }
summary { cursor:pointer; }
summary::-webkit-details-marker { color:var(--accent); }
.summary-title { font-weight:700; color:var(--accent); margin-right:10px; }
.summary-sub { color:var(--muted); font-size:0.92rem; }
pre { white-space:pre-wrap; margin:10px 0 0; color:#24414c; }
.plot-grid { display:grid; grid-template-columns:repeat(auto-fit, minmax(320px,1fr)); gap:14px; margin-top:12px; }
.plot-card { background:var(--panel); border:1px solid var(--border); border-radius:14px; padding:12px 12px 10px; }
.plot-card header { display:flex; align-items:center; justify-content:space-between; gap:8px; }
.plot-card h3 { margin:0; font-size:1rem; line-height:1.25; }
.path { margin:8px 0 8px; color:#2f4858; font-size:0.9rem; }
.path code { word-break:break-word; }
.note { margin:0 0 8px; color:var(--muted); font-size:0.88rem; }
.status { font-size:0.74rem; font-weight:700; border-radius:999px; padding:4px 8px; border:1px solid; }
.status.ok { color:var(--ok); border-color:#9cd2b7; background:#e8f5ee; }
.status.ko { color:var(--ko); border-color:#e4bbbb; background:#f9ecec; }
.status.missing { color:var(--missing); border-color:#cccccc; background:#f5f5f5; }
.status.unknown { color:var(--unknown); border-color:#e3c98d; background:#fbf2df; }
img { width:100%; border-radius:10px; border:1px solid #d6dddc; background:#ffffff; }
.carousel { display:grid; grid-template-columns:42px 1fr 42px; align-items:stretch; gap:8px; margin-top:12px; }
.carousel-track { display:grid; grid-auto-flow:column; grid-auto-columns:minmax(300px,420px); gap:12px; overflow-x:auto; scroll-snap-type:x mandatory; padding:2px 2px 6px; }
.carousel-track .plot-card { scroll-snap-align:start; }
.carousel-btn { border:1px solid var(--border); border-radius:10px; background:var(--panel); color:var(--accent); font-size:1.3rem; cursor:pointer; }
.carousel-btn:hover { background:#f0ece3; }
.table-wrap { overflow-x:auto; margin-top:12px; }
.summary-table { width:100%; border-collapse:collapse; background:var(--panel); border:1px solid var(--border); border-radius:10px; overflow:hidden; }
.summary-table th, .summary-table td { padding:8px 10px; border-bottom:1px solid #e7e2d9; text-align:left; font-size:0.88rem; }
.summary-table th { background:#f1ece2; color:#24414c; font-weight:700; }
@media (max-width:760px) { .carousel { grid-template-columns:1fr; } .carousel-btn { display:none; } }
</style>
"""
    script = """
<script>
(() => {
  const scrollStep = (track) => Math.max(280, Math.floor(track.clientWidth * 0.88));
  const move = (targetId, direction) => {
    const track = document.getElementById(targetId);
    if (!track) return;
    track.scrollBy({ left: direction * scrollStep(track), behavior: "smooth" });
  };
  document.querySelectorAll(".carousel-btn").forEach((btn) => {
    const target = btn.getAttribute("data-target");
    const dir = btn.classList.contains("next") ? 1 : -1;
    btn.addEventListener("click", () => move(target, dir));
  });
})();
</script>
"""

    with open(html_path, "w", encoding="utf-8") as h:
        h.write("<html><head><meta charset='utf-8'><title>He3pp Report</title>")
        h.write(style)
        h.write("</head><body><main class='wrap'>")
        h.write("<h1>He3pp Analysis Report</h1>")
        h.write(f"<p class='meta'><b>Variant:</b> {html.escape(s.VARIANT)} | <b>Species:</b> {html.escape(species)}</p>")
        if meta:
            h.write("<details open><summary><span class='summary-title'>Run Metadata</span><span class='summary-sub'>Execution context and merged configuration.</span></summary><pre>")
            h.write(html.escape(json.dumps(meta, indent=2)))
            h.write("</pre></details>")

        h.write(_render_summary_table(row_by_section["signal_tof"], row_by_section["signal_tpc"], tof_yield, tpc_yield, tpc_eff, tof_eff))
        for sec in section_order:
            title, subtitle = SECTION_DEFS[sec]
            sec_rows = row_by_section[sec]
            if not sec_rows:
                continue
            if sec in ("signal_tof", "signal_tpc"):
                h.write(_render_section(title, subtitle, sec_rows, f"{sec}-carousel"))
            else:
                h.write(_render_section(title, subtitle, sec_rows))
        h.write(script)
        h.write("</main></body></html>")

    return str(html_path)
