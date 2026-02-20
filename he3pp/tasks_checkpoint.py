from datetime import datetime
import logging

import ROOT

from .root_io import ensure_parent, expand
from .settings import RuntimeConfig
from .tasks_common import default_tpc_model_name


LOGGER = logging.getLogger("he3pp.tasks")


def checkpoint(
    systematics_file: str,
    data_ar_file: str,
    mc_file: str,
    mc_ar_file: str,
    signal_file: str,
    output_file: str | None = None,
    particle: str = "he3",
    *,
    runtime_config: RuntimeConfig,
) -> None:
    cfg = runtime_config
    LOGGER.info("checkpoint start particle=%s output=%s", particle, output_file if output_file else "auto")
    p = cfg.get_particle_config(particle)
    anti_name = str(p["anti_name"])
    hist_suffix = str(p["mc_hist_suffix"])
    syst = ROOT.TFile(expand(systematics_file))
    data_ar = ROOT.TFile(expand(data_ar_file))
    mc = ROOT.TFile(expand(mc_file))
    mc_ar = ROOT.TFile(expand(mc_ar_file))
    sig = ROOT.TFile(expand(signal_file))

    if output_file is None:
        output_file = f"checkpoint-{datetime.now().strftime('%d%m%y')}.root"

    out_name = expand(output_file)
    ensure_parent(out_name)
    out = ROOT.TFile(out_name, "RECREATE")
    out.cd()

    syst.Get("pubStat").Clone("published_stat").Write()
    syst.Get("pubSyst").Clone("published_syst").Write()
    syst.Get("fStatTPCA").Clone("tpc_spectrum_stat").Write()
    syst.Get("fSystTPCA").Clone("tpc_spectrum_syst").Write()
    syst.Get("fStatTOFA").Clone("tof_spectrum_stat").Write()
    syst.Get("fSystTOFA").Clone("tof_spectrum_syst").Write()
    mc.Get("nuclei/effTPCA").Clone("tpc_efficiency").Write()
    mc.Get("nuclei/effTOFA").Clone("tof_efficiency").Write()

    out.mkdir("MC")
    out.cd("MC")
    mc.Get(f"nuclei/genA{hist_suffix}").Clone("generated").Write()
    mc.Get(f"nuclei/TPCA{hist_suffix}").Clone("tpc_reconstructed").Write()
    mc.Get(f"nuclei/TOFA{hist_suffix}").Clone("tof_reconstructed").Write()
    mc_ar.Get("nuclei-spectra/spectra/hRecVtxZData").Clone("events_reconstructed").Write()

    out.mkdir("Data")
    out.cd("Data")
    data_ar.Get("nuclei-spectra/spectra/hRecVtxZData").Clone("events_reconstructed").Write()
    chosen_tpc_model_name = default_tpc_model_name(cfg)
    tpc_raw = sig.Get(f"nuclei/{anti_name}/TPConly/hTPConlyA0_{chosen_tpc_model_name}")
    if not tpc_raw:
        raise RuntimeError(
            f"Missing checkpoint TPC raw counts for model '{chosen_tpc_model_name}'. "
            f"Available configured models: {', '.join(cfg.tpc_function_names)}."
        )
    tpc_raw.Clone("tpc_rawcounts").Write()
    sig.Get(f"nuclei/{anti_name}/GausExp/hRawCountsA0").Clone("tof_rawcounts").Write()

    out.Close()
    LOGGER.info("checkpoint done particle=%s output=%s", particle, output_file if output_file else "auto")
