import logging
from typing import Any

import ROOT

from .root_io import build_rdf_from_ao2d, define_columns_for_data, ensure_parent, expand, write_hist
from .settings import *  # Loaded after runtime overrides in cli.run
from .tasks_common import collect_rresult_ptrs, run_graphs


LOGGER = logging.getLogger("he3pp.tasks")


def _book_data_species(df_all: Any, particle: str, skim: bool = False, tag: str = "") -> dict[str, Any]:
    p = get_particle_config(particle)
    df_base = df_all.Filter(p["base_sel"])
    df_primary = df_base.Filter(p["primary_sel"])
    df_secondary = df_base.Filter(p["secondary_sel"])
    nsigma = p["nsigma_col"]
    mass_cut = p["mass_cut_col"]
    pt_name = p["pt_col"]
    dmass = p["dmass_col"]
    nominal_mass = p["mass"]
    anti_label = str(p["label_antimatter"])
    matter_label = str(p["label_matter"])
    if int(p.get("tpc_anti_nbins", 0)) > 0:
        nbins_x = int(p["tpc_anti_nbins"])
        xmin = float(p["tpc_anti_xmin"])
        xmax = float(p["tpc_anti_xmax"])
        tpc_anti_model = ROOT.RDF.TH2DModel(f"fATPCcounts{tag}", f";#it{{p}}_{{T}}^{{rec}} (GeV/#it{{c}});{anti_label} n#sigma_{{TPC}};Counts", nbins_x, xmin, xmax, 100, -5, 5)
    else:
        tpc_anti_model = ROOT.RDF.TH2DModel(f"fATPCcounts{tag}", f";#it{{p}}_{{T}}^{{rec}} (GeV/#it{{c}});{anti_label} n#sigma_{{TPC}};Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -5, 5)
    tpc_mat_model = ROOT.RDF.TH2DModel(f"fMTPCcounts{tag}", f";#it{{p}}_{{T}}^{{rec}} (GeV/#it{{c}});{matter_label} n#sigma_{{TPC}};Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -5, 5)

    if skim:
        df_base.Filter(SKIM_SELECTION_TEMPLATE.format(nsigma=nsigma)).Snapshot("nucleiTree", "data/skimmed.root")

    matter_map = {"A": "!matter", "M": "matter"}
    h_dca_xy: dict[str, list[Any]] = {}
    h_dca_z: dict[str, list[Any]] = {}
    h_dca_xy_secondary: dict[str, list[Any]] = {}
    h_tpc: dict[str, list[Any]] = {}
    h_tof: dict[str, list[Any]] = {}

    for label, matter_sel in matter_map.items():
        h_dca_xy[label] = [df_primary.Filter(f"{matter_sel} && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"hDCAxy{label}He{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -0.2, 0.2), pt_name, "fDCAxy")]
        h_dca_z[label] = [df_primary.Filter(f"{matter_sel} && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"hDCAz{label}He{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -0.2, 0.2), pt_name, "fDCAz")]
        h_dca_xy_secondary[label] = [df_secondary.Filter(f"{matter_sel} && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"hDCAxySecondary{label}He{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -0.2, 0.2), pt_name, "fDCAxy")]

        tpc_sel = f"{matter_sel} && {mass_cut}" if p["tpc_apply_mass_cut"] else matter_sel
        tpc_model = tpc_anti_model if label == "A" else tpc_mat_model
        tpc_pt_name = p["tpc_pt_col_anti"] if label == "A" else p["tpc_pt_col_matter"]
        h_tpc[label] = [df_primary.Filter(tpc_sel).Histo2D(tpc_model, tpc_pt_name, nsigma)]

    nsigma_tof_cut = p["tof_nsigma_cut"]
    for label, matter_sel in matter_map.items():
        species_label = anti_label if label == "A" else matter_label
        h_tof[label] = [df_primary.Filter(f"{matter_sel} && std::abs({nsigma}) < {nsigma_tof_cut}").Histo2D(ROOT.RDF.TH2DModel(f"f{label}TOFsignal{tag}", f";#it{{p}}_{{T}}^{{rec}} (GeV/#it{{c}});m_{{TOF}}-m_{{{species_label}}};Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -0.9, 1.1), pt_name, dmass)]

    df_primary = df_primary.Define(f"tofMassDeltaPOI{tag}", f"tofMass - {nominal_mass}")
    tof_delta_mass_range = (float(p["tof_mass_range_min"]), float(p["tof_mass_range_max"]))
    h_tof_mass_vs_tpc_nsigma: dict[str, list[Any]] = {"A": [], "M": []}
    for i_pt in range(N_PT_BINS):
        pt_low = PT_BINS[i_pt]
        pt_high = PT_BINS[i_pt + 1]
        for label, matter_sel in matter_map.items():
            h_tof_mass_vs_tpc_nsigma[label].append(
                df_primary.Filter(f"{matter_sel} && hasTOF && std::abs({nsigma}) < 6 && {pt_name} >= {pt_low} && {pt_name} < {pt_high}").Histo2D(
                    ROOT.RDF.TH2DModel(
                        f"hTOFMassVsTPCnsigma{label}_pt{i_pt:02d}{tag}",
                        ";n#sigma_{TPC};m_{TOF}-m_{0} (GeV/#it{c}^{2});Counts",
                        120,
                        -6.0,
                        6.0,
                        120,
                        tof_delta_mass_range[0],
                        tof_delta_mass_range[1],
                    ),
                    nsigma,
                    f"tofMassDeltaPOI{tag}",
                )
            )

    i_trial = 0
    if p["trial_enabled"]:
        for cut_dca_z in CUT_NAMES["nsigmaDCAz"]:
            df_dca_z = df_base.Filter(f"std::abs(nsigmaDCAz) < {cut_dca_z}")
            for cut_tpc in CUT_NAMES["fTPCnCls"]:
                df_tpc = df_dca_z.Filter(f"fTPCnCls > {cut_tpc}")
                for cut_its in CUT_NAMES["nITScls"]:
                    df_its = df_tpc.Filter(f"nITScls >= {cut_its}")
                    for label, matter_sel in matter_map.items():
                        species_label = anti_label if label == "A" else matter_label
                        h_dca_xy[label].append(df_its.Filter(f"{matter_sel} && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"hDCAxy{label}He{i_trial}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{xy} (cm);Counts", N_PT_BINS, PT_BIN_ARRAY, 560, -0.7, 0.7), pt_name, "fDCAxy"))
                        h_dca_z[label].append(df_its.Filter(f"{matter_sel} && {nsigma} > -0.5 && {nsigma} < 3 && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"hDCAz{label}He{i_trial}{tag}", ";#it{p}_{T}^{rec} (GeV/#it{c});DCA_{z} (cm);Counts", N_PT_BINS, PT_BIN_ARRAY, 560, -0.7, 0.7), pt_name, "fDCAz"))
                        h_tpc[label].append(df_its.Filter(f"{matter_sel} && {p['trial_dca_sel']} && {mass_cut}").Histo2D(ROOT.RDF.TH2DModel(f"f{label}TPCcounts{i_trial}{tag}", f";#it{{p}}_{{T}}^{{rec}} (GeV/#it{{c}});{species_label} n#sigma_{{TPC}};Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -5, 5), pt_name, nsigma))
                        h_tof[label].append(df_its.Filter(f"{matter_sel} && {p['trial_dca_sel']} && std::abs({nsigma}) < {nsigma_tof_cut}").Histo2D(ROOT.RDF.TH2DModel(f"f{label}TOFsignal{i_trial}{tag}", f";#it{{p}}_{{T}}^{{rec}} (GeV/#it{{c}});m_{{TOF}}-m_{{{species_label}}};Counts", N_PT_BINS, PT_BIN_ARRAY, 100, -0.9, 1.1), pt_name, dmass))
                    i_trial += 1

    return {
        "particle": particle,
        "hist_suffix": str(p["mc_hist_suffix"]),
        "h_tpc": h_tpc,
        "h_tof": h_tof,
        "h_dca_xy": h_dca_xy,
        "h_dca_z": h_dca_z,
        "h_dca_xy_secondary": h_dca_xy_secondary,
        "h_tof_mass_vs_tpc_nsigma": h_tof_mass_vs_tpc_nsigma,
        "i_trial": i_trial,
    }


def _write_data_bundle(bundle: dict[str, Any], output_file: str) -> None:
    output_file = expand(output_file)
    ensure_parent(output_file)
    out = ROOT.TFile(output_file, "recreate")
    d0 = out.mkdir("nuclei")
    d0.cd()
    suffix = str(bundle.get("hist_suffix", "He3"))
    for label in ("A", "M"):
        write_hist(bundle["h_tpc"][label][0], f"f{label}TPCcounts")
        write_hist(bundle["h_tof"][label][0], f"f{label}TOFsignal")
        write_hist(bundle["h_dca_xy"][label][0], f"hDCAxy{label}{suffix}")
        write_hist(bundle["h_dca_z"][label][0], f"hDCAz{label}{suffix}")
    for label in ("M", "A"):
        write_hist(bundle["h_dca_xy_secondary"][label][0], f"hDCAxySecondary{label}{suffix}")
    for i_pt in range(N_PT_BINS):
        for label in ("A", "M"):
            write_hist(bundle["h_tof_mass_vs_tpc_nsigma"][label][i_pt], f"hTOFMassVsTPCnsigma{label}_pt{i_pt:02d}")

    for idx in range(bundle["i_trial"]):
        dtrial = out.mkdir(f"nuclei{idx}")
        dtrial.cd()
        for label in ("A", "M"):
            write_hist(bundle["h_tpc"][label][idx + 1], f"f{label}TPCcounts")
            write_hist(bundle["h_tof"][label][idx + 1], f"f{label}TOFsignal")
            write_hist(bundle["h_dca_xy"][label][idx + 1], f"hDCAxy{label}{suffix}")
            write_hist(bundle["h_dca_z"][label][idx + 1], f"hDCAz{label}{suffix}")
    out.Close()


def analyse_data(input_file: str, output_file: str, particle: str, skim: bool = False, draw: bool = False) -> None:
    LOGGER.info("analyse_data start particle=%s input=%s output=%s", particle, input_file, output_file)
    ROOT.gStyle.SetOptStat(0)
    rdf = build_rdf_from_ao2d("O2nucleitable", input_file)
    df_all = define_columns_for_data(rdf)
    bundle = _book_data_species(df_all, particle, skim=skim, tag=f"_{particle}")
    run_graphs(
        collect_rresult_ptrs(bundle["h_tpc"])
        + collect_rresult_ptrs(bundle["h_tof"])
        + collect_rresult_ptrs(bundle["h_dca_xy"])
        + collect_rresult_ptrs(bundle["h_dca_z"])
        + collect_rresult_ptrs(bundle["h_dca_xy_secondary"])
        + collect_rresult_ptrs(bundle["h_tof_mass_vs_tpc_nsigma"])
    )
    if draw:
        draw_order = [
            (bundle["h_tpc"], ("A", "M")),
            (bundle["h_tof"], ("A", "M")),
            (bundle["h_dca_xy"], ("A", "M")),
            (bundle["h_dca_z"], ("A", "M")),
            (bundle["h_dca_xy_secondary"], ("M", "A")),
        ]
        for hist_map, labels in draw_order:
            for label in labels:
                hist = hist_map[label][0]
                c = ROOT.TCanvas()
                c.SetRightMargin(0.15)
                hist.DrawClone("colz")
    _write_data_bundle(bundle, output_file)
    LOGGER.info("analyse_data done output=%s", output_file)
