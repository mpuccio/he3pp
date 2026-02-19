from __future__ import annotations

from dataclasses import dataclass, field
import html
import json
import os
from pathlib import Path
import re
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
    overlay_path: str = ""


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
    "tof_tpc_2d": ("TOF Mass vs TPC nσ", "Per-bin 2D correlation of m_TOF and nσ_TPC."),
    "efficiency": ("Efficiency", "TPC and TOF efficiencies with binomial uncertainty."),
    "pt_resolution": ("pT Resolution", "Comparison of reconstructed and generated transverse momentum."),
    "corrected_spectrum": ("Normalized Corrected Spectrum", "Normalized spectra with statistical points and systematic bands."),
}
DEFAULT_SECTIONS = ["signal_tof", "signal_tpc", "tof_tpc_2d", "efficiency", "pt_resolution", "corrected_spectrum"]


def _mkdir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def _get(root_file: ROOT.TFile, path: str) -> Any | None:
    obj = root_file.Get(path)
    return obj if obj else None


def _species_letter(species: str) -> str:
    return "A" if species.startswith("anti") else "M"


def _preflight_inputs(signal_file_path: str, mc_file_path: str, systematics_file_path: str, data_file_path: str | None) -> dict[str, str]:
    required = {
        "signal": expand(signal_file_path),
        "mc": expand(mc_file_path),
        "systematics": expand(systematics_file_path),
    }
    if data_file_path:
        required["data"] = expand(data_file_path)
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
        obj.SetStats(0)
        obj.Draw("colz")
    elif cname.startswith("TH"):
        obj.SetStats(0)
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


def _draw_overlay_to_png(
    stat_obj: Any,
    syst_obj: Any,
    png_path: str,
    y_title: str | None = None,
    y_range: tuple[float, float] | None = None,
    logy: bool = False,
) -> bool:
    if stat_obj is None or syst_obj is None:
        return False
    _mkdir(str(Path(png_path).parent))
    c = ROOT.TCanvas("c_report_overlay", "c_report_overlay", 950, 700)
    if not stat_obj.ClassName().startswith("TH1") or not syst_obj.ClassName().startswith("TH1"):
        c.Close()
        return False
    syst = syst_obj.Clone("syst_report_overlay")
    stat = stat_obj.Clone("stat_report_overlay")
    syst.SetStats(0)
    stat.SetStats(0)
    syst.SetFillColorAlpha(ROOT.kOrange + 1, 0.35)
    syst.SetLineColor(ROOT.kOrange + 2)
    syst.SetMarkerSize(0)
    stat.SetLineColor(ROOT.kBlack)
    stat.SetMarkerStyle(20)
    stat.SetMarkerSize(0.9)
    if y_title:
        syst.GetYaxis().SetTitle(y_title)
        stat.GetYaxis().SetTitle(y_title)
    if y_range:
        ymin, ymax = y_range
        syst.SetMinimum(ymin)
        syst.SetMaximum(ymax)
        stat.SetMinimum(ymin)
        stat.SetMaximum(ymax)
    if logy:
        c.SetLogy()
    syst.Draw("E2")
    stat.Draw("E1 SAME")
    leg = ROOT.TLegend(0.63, 0.76, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.AddEntry(stat, "Stat.", "lep")
    leg.AddEntry(syst, "Syst.", "f")
    leg.Draw()
    c.SaveAs(png_path)
    c.Close()
    return True


def _draw_pair_hist_to_png(
    h1_obj: Any,
    h2_obj: Any,
    png_path: str,
    y_title: str | None = None,
    labels: tuple[str, str] = ("A", "B"),
) -> bool:
    if h1_obj is None or h2_obj is None:
        return False
    if not h1_obj.ClassName().startswith("TH1") or not h2_obj.ClassName().startswith("TH1"):
        return False
    _mkdir(str(Path(png_path).parent))
    c = ROOT.TCanvas("c_report_pair", "c_report_pair", 950, 700)
    h1 = h1_obj.Clone("report_pair_h1")
    h2 = h2_obj.Clone("report_pair_h2")
    h1.SetStats(0)
    h2.SetStats(0)
    h1.SetLineColor(ROOT.kBlue + 1)
    h1.SetMarkerColor(ROOT.kBlue + 1)
    h1.SetMarkerStyle(20)
    h1.SetMarkerSize(0.9)
    h2.SetLineColor(ROOT.kRed + 1)
    h2.SetMarkerColor(ROOT.kRed + 1)
    h2.SetMarkerStyle(24)
    h2.SetMarkerSize(0.9)
    ymax = max(_hist_visible_max(h1), _hist_visible_max(h2))
    h1.SetMinimum(0.0)
    h1.SetMaximum(1.15 * ymax if ymax > 0 else 1.0)
    if y_title:
        h1.GetYaxis().SetTitle(y_title)
        h2.GetYaxis().SetTitle(y_title)
    h1.Draw("E1")
    h2.Draw("E1 SAME")
    leg = ROOT.TLegend(0.63, 0.76, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.AddEntry(h1, labels[0], "lep")
    leg.AddEntry(h2, labels[1], "lep")
    leg.Draw()
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
            try:
                order = int(n.split("_", 1)[1])
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


def _tof_tpc_2d_paths(data_file: ROOT.TFile, species: str) -> list[tuple[int, str]]:
    base = data_file.Get("nuclei")
    if not base:
        return []
    label = _species_letter(species)
    out: list[tuple[int, str]] = []
    for k in base.GetListOfKeys():
        n = k.GetName()
        if not n.startswith(f"hTOFMassVsTPCnsigma{label}_pt"):
            continue
        try:
            bidx = int(n.rsplit("pt", 1)[1])
        except ValueError:
            continue
        out.append((bidx, f"nuclei/{n}"))
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
    p_value = min(1.0, 2.0 * min(p_right, 1.0 - p_right)) if fit_tail == "two" else p_right
    return {"chi2_ndf": chi2_ndf, "dof": float(dof), "p_value": p_value}


def _status_from_metrics(available: bool, metrics: dict[str, float], alpha: float) -> tuple[str, str]:
    if not available:
        return "MISSING", "missing"
    if not metrics:
        return "UNK", "unknown"
    return ("OK", "ok") if metrics["p_value"] >= alpha else ("KO", "ko")


def _hist_visible_max(h: Any) -> float:
    if not h or not h.ClassName().startswith("TH1"):
        return 0.0
    vmax = 0.0
    for i in range(1, h.GetNbinsX() + 1):
        vmax = max(vmax, float(h.GetBinContent(i) + h.GetBinError(i)))
    return vmax


def _hist_positive_min(h: Any) -> float | None:
    if not h or not h.ClassName().startswith("TH1"):
        return None
    vmin: float | None = None
    for i in range(1, h.GetNbinsX() + 1):
        val = float(h.GetBinContent(i))
        if val > 0.0 and (vmin is None or val < vmin):
            vmin = val
    return vmin


def _render_card(row: RenderedItem) -> str:
    it = row.item
    safe_title = html.escape(it.title)
    status_text = row.status
    if row.metrics and status_text:
        status_text = f"{status_text} (p={row.metrics['p_value']:.3f})"
    path_text = html.escape(it.object_path if not it.overlay_path else f"{it.object_path} + {it.overlay_path}")
    header = f"<header><h3>{safe_title}</h3>"
    if status_text:
        header += f"<span class='status {row.status_class}'>{status_text}</span>"
    header += "</header>"
    out = ["<article class='plot-card'>", header, f"<p class='path'><code>{path_text}</code></p>"]
    if row.note:
        out.append(f"<p class='note'>{html.escape(row.note)}</p>")
    if row.available:
        out.append(f"<img loading='lazy' src='{html.escape(it.output_png)}' alt='{safe_title}'>")
    out.append("</article>")
    return "".join(out)


def _section_id(title: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
    return f"section-{slug or 'section'}"


def _render_section(title: str, subtitle: str, rows: list[RenderedItem], carousel_id: str | None = None) -> str:
    section_id = _section_id(title)
    out = [
        f"<details class='report-section' id='{html.escape(section_id)}' open>",
        f"<summary><span class='summary-title'>{html.escape(title)}</span><span class='summary-sub'>{html.escape(subtitle)}</span></summary>",
    ]
    if carousel_id:
        cid = html.escape(carousel_id)
        out.append(
            f"<div class='carousel' data-carousel='{cid}'>"
            f"<button class='carousel-btn prev' type='button' aria-label='Previous' data-target='{cid}'>&lsaquo;</button>"
            f"<div class='carousel-track' id='{cid}'>"
        )
        out.extend(_render_card(r) for r in rows)
        out.append("</div>" f"<button class='carousel-btn next' type='button' aria-label='Next' data-target='{cid}'>&rsaquo;</button>" "</div>")
    else:
        out.append("<div class='plot-grid'>")
        out.extend(_render_card(r) for r in rows)
        out.append("</div>")
    out.append("</details>")
    return "".join(out)


def _fmt_opt(value: float | None, fmt: str) -> str:
    return "-" if value is None else format(value, fmt)


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
    all_bins = list(range(tof_yield.GetNbinsX())) if tof_yield else sorted(set(p_tof) | set(p_tpc))
    axis = tof_yield.GetXaxis() if tof_yield else (tpc_yield.GetXaxis() if tpc_yield else None)
    out = [
        "<details class='report-section' id='section-fit-and-efficiency-summary'>",
        "<summary><span class='summary-title'>Fit And Efficiency Summary</span><span class='summary-sub'>TOF/TPC yields, efficiencies and fit p-values by pT bin.</span></summary>",
        "<div class='table-wrap'><table class='summary-table'><thead><tr>"
        "<th>Bin</th><th>pT Range</th><th>TOF Yield</th><th>TPC Yield</th><th>TPC Eff</th><th>TOF Eff</th><th>TOF p-value</th><th>TPC p-value</th>"
        "</tr></thead><tbody>",
    ]
    for b in all_bins:
        p_range = "-"
        if axis:
            p_range = f"{axis.GetBinLowEdge(b + 1):.2f} - {axis.GetBinUpEdge(b + 1):.2f}"
        out.append(
            "<tr>"
            f"<td>{b}</td><td>{html.escape(p_range)}</td>"
            f"<td>{_fmt_opt(tof_yield.GetBinContent(b + 1) if tof_yield else None, '.2f')}</td>"
            f"<td>{_fmt_opt(tpc_yield.GetBinContent(b + 1) if tpc_yield else None, '.2f')}</td>"
            f"<td>{_fmt_opt(tpc_eff.GetBinContent(b + 1) if tpc_eff else None, '.4f')}</td>"
            f"<td>{_fmt_opt(tof_eff.GetBinContent(b + 1) if tof_eff else None, '.4f')}</td>"
            f"<td>{_fmt_opt(p_tof.get(b), '.3f')}</td><td>{_fmt_opt(p_tpc.get(b), '.3f')}</td>"
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
    data_file_path: str | None = None,
) -> str:
    report_dir = expand(report_dir)
    _mkdir(str(Path(report_dir) / "assets"))
    ROOT.gStyle.SetOptStat(0)
    section_order = _normalize_sections(sections)
    fit_tail = "two" if str(fit_tail).lower() == "two" else "single"

    resolved = _preflight_inputs(signal_file_path, mc_file_path, systematics_file_path, data_file_path)
    sig_file = ROOT.TFile(resolved["signal"])
    mc_file = ROOT.TFile(resolved["mc"])
    syst_file = ROOT.TFile(resolved["systematics"])
    data_file = ROOT.TFile(resolved["data"]) if "data" in resolved else None

    items: list[ReportItem] = []
    for bidx, p in _signal_tof_bin_paths(sig_file, species):
        items.append(ReportItem("signal_tof", f"TOF Signal Bin {bidx}", signal_file_path, p, f"assets/signal_tof_bin_{bidx:02d}.png", bin_index=bidx))
    for bidx, p in _signal_tpc_bin_paths(sig_file, species, tpc_signal_model):
        items.append(ReportItem("signal_tpc", f"TPC Signal Bin {bidx}", signal_file_path, p, f"assets/signal_tpc_bin_{bidx:02d}.png", note=f"Model: {tpc_signal_model}", bin_index=bidx))
    if data_file:
        for bidx, p in _tof_tpc_2d_paths(data_file, species):
            items.append(ReportItem("tof_tpc_2d", f"m_TOF vs nσ_TPC Bin {bidx}", resolved["data"], p, f"assets/tof_mass_vs_tpc_nsigma_bin_{bidx:02d}.png", bin_index=bidx))

    raw_name = "hRawCountsA0" if species == "antihe3" else "hRawCountsM0"
    items.append(
        ReportItem(
            "efficiency",
            "TPC/TOF Efficiency Overlay",
            mc_file_path,
            "nuclei/effTPCA" if species == "antihe3" else "nuclei/effTPCM",
            "assets/eff_tpc_tof_overlay.png",
            note="Uncertainty: binomial (both)",
            overlay_path="nuclei/effTOFA" if species == "antihe3" else "nuclei/effTOFM",
        )
    )
    items.append(ReportItem("pt_resolution", "pT(rec)-pT(sim) Uncorrected", mc_file_path, "nuclei/hDeltaPtHe3", "assets/delta_pt_uncorrected.png"))
    items.append(ReportItem("pt_resolution", "pT(rec)-pT(sim) Corrected", mc_file_path, "nuclei/hDeltaPtCorrHe3", "assets/delta_pt_corrected.png"))
    items.append(ReportItem("corrected_spectrum", "Normalized Corrected Spectrum TPC", systematics_file_path, "fStatTPCA" if species == "antihe3" else "fStatTPCM", "assets/corrected_tpc_norm.png", overlay_path="fSystTPCA" if species == "antihe3" else "fSystTPCM"))
    items.append(ReportItem("corrected_spectrum", "Normalized Corrected Spectrum TOF", systematics_file_path, "fStatTOFA" if species == "antihe3" else "fStatTOFM", "assets/corrected_tof_norm.png", overlay_path="fSystTOFA" if species == "antihe3" else "fSystTOFM"))

    file_map = {signal_file_path: sig_file, mc_file_path: mc_file, systematics_file_path: syst_file}
    if data_file:
        file_map[resolved["data"]] = data_file

    corrected_stats = [
        _get(syst_file, "fStatTPCA" if species == "antihe3" else "fStatTPCM"),
        _get(syst_file, "fSystTPCA" if species == "antihe3" else "fSystTPCM"),
        _get(syst_file, "fStatTOFA" if species == "antihe3" else "fStatTOFM"),
        _get(syst_file, "fSystTOFA" if species == "antihe3" else "fSystTOFM"),
    ]
    corr_max = max(_hist_visible_max(h) for h in corrected_stats)
    corr_pos = [v for v in (_hist_positive_min(h) for h in corrected_stats) if v is not None]
    corr_ymin = 0.5 * min(corr_pos) if corr_pos else 1e-12
    corr_ymax = 1.15 * corr_max if corr_max > 0 else 1.0
    corr_y_title = "1/N_{ev} d^{2}#it{N}/d#it{y}d#it{p}_{T} (GeV/#it{c})^{-1}"

    rows: list[RenderedItem] = []
    for it in items:
        f = file_map[it.source_file]
        obj = _get(f, it.object_path)
        if it.section == "efficiency" and it.overlay_path:
            ok = _draw_pair_hist_to_png(
                obj,
                _get(f, it.overlay_path),
                str(Path(report_dir) / it.output_png),
                y_title="Efficiency #times Acceptance",
                labels=("TPC", "TOF"),
            )
        elif it.overlay_path:
            ok = _draw_overlay_to_png(
                obj,
                _get(f, it.overlay_path),
                str(Path(report_dir) / it.output_png),
                y_title=corr_y_title if it.section == "corrected_spectrum" else None,
                y_range=(corr_ymin, corr_ymax) if it.section == "corrected_spectrum" else None,
                logy=(it.section == "corrected_spectrum"),
            )
        else:
            ok = _draw_to_png(obj, str(Path(report_dir) / it.output_png))
        metrics = _chi2_metrics(obj, fit_n_parameters, fit_tail) if it.section in ("signal_tof", "signal_tpc") else {}
        if it.section in ("signal_tof", "signal_tpc"):
            status, status_class = _status_from_metrics(ok, metrics, fit_alpha)
        else:
            status, status_class = (("", "")) if ok else (("MISSING", "missing"))
        rows.append(RenderedItem(it, ok, status, status_class, note=it.note, metrics=metrics))

    row_by_section = {sec: [r for r in rows if r.item.section == sec] for sec in SECTION_DEFS}
    letter = _species_letter(species)
    tof_yield = _get(sig_file, f"nuclei/{species}/GausExp/{raw_name}")
    tpc_yield = _get(sig_file, f"nuclei/{species}/TPConly/hTPConly{letter}0_{tpc_signal_model}")
    tpc_eff = _get(mc_file, "nuclei/effTPCA" if species == "antihe3" else "nuclei/effTPCM")
    tof_eff = _get(mc_file, "nuclei/effTOFA" if species == "antihe3" else "nuclei/effTOFM")

    meta = {}
    if metadata_path and os.path.exists(expand(metadata_path)):
        with open(expand(metadata_path), "r", encoding="utf-8") as f:
            meta = json.load(f)

    style = """
<style>
:root { --bg:#f6f3ec; --panel:#fffef8; --ink:#1f2a30; --muted:#5f6c74; --accent:#005f73; --ok:#127143; --ko:#9f2a2a; --missing:#7a7a7a; --unknown:#885f00; --border:#d9d2c4; }
* { box-sizing:border-box; } body { margin:0; font-family:"Avenir Next","Segoe UI",sans-serif; color:var(--ink); background:radial-gradient(circle at 20% 15%, #e4efe9 0%, transparent 30%),radial-gradient(circle at 85% 10%, #f4dfc8 0%, transparent 30%),var(--bg); }
.wrap { max-width:1260px; margin:0 auto; padding:26px 20px 48px; } h1 { margin:0 0 6px; font-size:2.05rem; letter-spacing:0.02em; } .meta { color:var(--muted); margin-bottom:14px; font-size:0.98rem; }
.quick-nav { display:flex; flex-wrap:wrap; gap:8px; margin:0 0 18px; padding:10px 0 4px; border-top:1px solid #e3ddd2; border-bottom:1px solid #e3ddd2; } .quick-nav a { text-decoration:none; color:var(--accent); border:1px solid #cfe0db; background:#f2faf8; padding:5px 10px; border-radius:999px; font-size:0.82rem; font-weight:600; }
details { background:var(--panel); border:1px solid var(--border); border-radius:14px; padding:12px 15px; margin-bottom:18px; box-shadow:0 1px 0 #ece6da; } summary { cursor:pointer; list-style:none; } summary::-webkit-details-marker { display:none; }
.summary-title { font-weight:750; color:var(--accent); margin-right:10px; } .summary-sub { color:var(--muted); font-size:0.92rem; } pre { white-space:pre-wrap; margin:10px 0 0; color:#24414c; max-height:460px; overflow:auto; }
.plot-grid { display:grid; grid-template-columns:repeat(auto-fit, minmax(330px,1fr)); gap:14px; margin-top:12px; } .plot-card { background:#fffdfa; border:1px solid #ddd4c5; border-radius:14px; padding:12px 12px 10px; box-shadow:0 1px 0 #efe7da; }
.plot-card header { display:flex; align-items:center; justify-content:space-between; gap:8px; margin-bottom:3px; } .plot-card h3 { margin:0; font-size:1rem; line-height:1.25; } .path { margin:8px 0 8px; color:#2f4858; font-size:0.85rem; } .path code { word-break:break-word; } .note { margin:0 0 8px; color:var(--muted); font-size:0.86rem; }
.status { font-size:0.74rem; font-weight:700; border-radius:999px; padding:4px 8px; border:1px solid; } .status.ok { color:var(--ok); border-color:#9cd2b7; background:#e8f5ee; } .status.ko { color:var(--ko); border-color:#e4bbbb; background:#f9ecec; } .status.missing { color:var(--missing); border-color:#cccccc; background:#f5f5f5; } .status.unknown { color:var(--unknown); border-color:#e3c98d; background:#fbf2df; }
img { width:100%; border-radius:10px; border:1px solid #d6dddc; background:#ffffff; }
.carousel { display:grid; grid-template-columns:40px 1fr 40px; align-items:stretch; gap:8px; margin-top:12px; } .carousel-track { display:grid; grid-auto-flow:column; grid-auto-columns:minmax(320px,430px); gap:12px; overflow-x:auto; scroll-snap-type:x mandatory; padding:2px 2px 6px; scrollbar-width:thin; } .carousel-track .plot-card { scroll-snap-align:start; } .carousel-btn { border:1px solid #cec4b3; border-radius:10px; background:#fff7ea; color:var(--accent); font-size:1.3rem; cursor:pointer; font-weight:700; }
.table-wrap { overflow-x:auto; margin-top:12px; } .summary-table { width:100%; border-collapse:collapse; background:var(--panel); border:1px solid var(--border); border-radius:10px; overflow:hidden; } .summary-table th, .summary-table td { padding:8px 10px; border-bottom:1px solid #e7e2d9; text-align:left; font-size:0.86rem; } .summary-table th { background:#f1ece2; color:#24414c; font-weight:700; position:sticky; top:0; z-index:1; }
@media (max-width:760px) { .quick-nav { display:none; } .carousel { grid-template-columns:1fr; } .carousel-btn { display:none; } }
</style>
"""
    script = """
<script>
(() => {
  const scrollStep = (track) => Math.max(280, Math.floor(track.clientWidth * 0.88));
  const move = (targetId, direction) => { const track = document.getElementById(targetId); if (!track) return; track.scrollBy({ left: direction * scrollStep(track), behavior: "smooth" }); };
  document.querySelectorAll(".carousel-btn").forEach((btn) => { const target = btn.getAttribute("data-target"); const dir = btn.classList.contains("next") ? 1 : -1; btn.addEventListener("click", () => move(target, dir)); });
  document.querySelectorAll(".quick-nav a[href^='#']").forEach((link) => {
    link.addEventListener("click", () => {
      const target = document.querySelector(link.getAttribute("href"));
      if (target && target.tagName.toLowerCase() === "details") target.open = true;
    });
  });
})();
</script>
"""

    html_path = Path(report_dir) / "index.html"
    with open(html_path, "w", encoding="utf-8") as h:
        h.write("<html><head><meta charset='utf-8'><title>He3pp Report</title>")
        h.write(style)
        h.write("</head><body><main class='wrap'>")
        h.write("<h1>He3pp Analysis Report</h1>")
        h.write(f"<p class='meta'><b>Variant:</b> {html.escape(s.VARIANT)} | <b>Species:</b> {html.escape(species)}</p>")
        nav_items: list[str] = ["<a href='#section-fit-and-efficiency-summary'>Summary</a>"]
        for sec in section_order:
            title, _subtitle = SECTION_DEFS.get(sec, (sec, ""))
            sec_id = _section_id(title)
            nav_items.append(f"<a href='#{html.escape(sec_id)}'>{html.escape(title)}</a>")
        h.write(f"<nav class='quick-nav'>{''.join(nav_items)}</nav>")
        if meta:
            h.write("<details><summary><span class='summary-title'>Run Metadata</span><span class='summary-sub'>Execution context and merged configuration.</span></summary><pre>")
            h.write(html.escape(json.dumps(meta, indent=2)))
            h.write("</pre></details>")
        h.write(_render_summary_table(row_by_section["signal_tof"], row_by_section["signal_tpc"], tof_yield, tpc_yield, tpc_eff, tof_eff))
        for sec in section_order:
            sec_rows = row_by_section[sec]
            if not sec_rows:
                continue
            title, subtitle = SECTION_DEFS[sec]
            if sec in ("signal_tof", "signal_tpc", "tof_tpc_2d"):
                h.write(_render_section(title, subtitle, sec_rows, f"{sec}-carousel"))
            else:
                h.write(_render_section(title, subtitle, sec_rows))
        h.write(script)
        h.write("</main></body></html>")
    return str(html_path)


def generate_dual_report_index(
    report_dir: str,
    entries: list[dict[str, str]],
    metadata_path: str | None = None,
) -> str:
    report_dir = expand(report_dir)
    _mkdir(report_dir)
    meta = {}
    if metadata_path and os.path.exists(expand(metadata_path)):
        with open(expand(metadata_path), "r", encoding="utf-8") as f:
            meta = json.load(f)

    cards: list[str] = []
    for e in entries:
        label = html.escape(e.get("label", e.get("species", "species")))
        species = html.escape(e.get("species", ""))
        href = html.escape(e.get("href", ""))
        cards.append(
            "<article class='card'>"
            f"<h2>{label}</h2>"
            f"<p><b>Species:</b> {species}</p>"
            f"<p><a href='{href}'>Open Report</a></p>"
            "</article>"
        )

    style = """
<style>
:root { --bg:#f6f3ec; --panel:#fffef8; --ink:#1f2a30; --muted:#5f6c74; --accent:#005f73; --border:#d9d2c4; }
* { box-sizing:border-box; }
body { margin:0; font-family:"Avenir Next","Segoe UI",sans-serif; color:var(--ink); background:var(--bg); }
.wrap { max-width:980px; margin:0 auto; padding:26px 20px 40px; }
h1 { margin:0 0 8px; }
.meta { color:var(--muted); margin-bottom:18px; }
.grid { display:grid; grid-template-columns:repeat(auto-fit, minmax(280px,1fr)); gap:14px; }
.card { background:var(--panel); border:1px solid var(--border); border-radius:12px; padding:14px; }
.card h2 { margin:0 0 8px; color:var(--accent); font-size:1.1rem; }
a { color:var(--accent); text-decoration:none; font-weight:700; }
details { margin-top:18px; background:var(--panel); border:1px solid var(--border); border-radius:12px; padding:12px 14px; }
summary { cursor:pointer; }
pre { white-space:pre-wrap; margin:10px 0 0; }
</style>
"""

    html_path = Path(report_dir) / "index.html"
    with open(html_path, "w", encoding="utf-8") as h:
        h.write("<html><head><meta charset='utf-8'><title>He3pp Dual Report</title>")
        h.write(style)
        h.write("</head><body><main class='wrap'>")
        h.write("<h1>He3pp Dual Report</h1>")
        h.write(f"<p class='meta'><b>Variant:</b> {html.escape(s.VARIANT)} | <b>Species reports:</b> {len(entries)}</p>")
        h.write(f"<section class='grid'>{''.join(cards)}</section>")
        if meta:
            h.write("<details><summary><b>Run Metadata</b></summary><pre>")
            h.write(html.escape(json.dumps(meta, indent=2)))
            h.write("</pre></details>")
        h.write("</main></body></html>")
    return str(html_path)
