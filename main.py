from multiprocessing import Pool

import click
import numpy as np
from tqdm import tqdm

from src import DigestParams, build_template, render_gel, simulate

DEFAULTS = DigestParams()


def run_batch(seed, template_kwargs, params, sample_times):
    rng = np.random.default_rng(seed)
    template = build_template(rng=rng, **template_kwargs)
    return simulate(template, params, sample_times, rng)


def _run_batch_star(job):
    return run_batch(*job)


@click.command(context_settings={"show_default": True})
# Chromatin template
@click.option("--n-nucleosomes", type=click.IntRange(min=1), default=100, help="Nucleosomes per DNA molecule.")
@click.option("--core-length", type=click.IntRange(min=1), default=147, help="DNA wrapped per nucleosome core (bp).")
@click.option("--linker-mean", type=click.FloatRange(min=0), default=30, help="Mean linker length (bp).")
@click.option("--linker-sd", type=click.FloatRange(min=0), default=8, help="Linker length standard deviation (bp).")
@click.option("--nfr-prob", type=click.FloatRange(0, 1), default=0.0, help="Probability a linker is a nucleosome-free region.")
@click.option("--nfr-length", type=click.IntRange(min=0), default=300, help="Nucleosome-free region length (bp).")
# Digestion kinetics
@click.option("--k-linker", type=click.FloatRange(min=0), default=DEFAULTS.k_linker, help="Endonuclease cuts per bp per min in linkers.")
@click.option("--k-nuc", type=click.FloatRange(min=0), default=DEFAULTS.k_nuc, help="Endonuclease cuts per bp per min in cores.")
@click.option("--k-exo-linker", type=click.FloatRange(min=0), default=DEFAULTS.k_exo_linker, help="bp trimmed per end per min in linkers.")
@click.option("--k-exo-nuc", type=click.FloatRange(min=0), default=DEFAULTS.k_exo_nuc, help="bp trimmed per end per min in cores.")
@click.option("--dt", type=click.FloatRange(min=0, min_open=True), default=DEFAULTS.dt, help="Simulation time step (min).")
@click.option("--min-len", type=click.IntRange(min=0), default=DEFAULTS.min_len, help="Fragments shorter than this are lost (bp).")
@click.option("--time", "times", type=click.FloatRange(min=0), multiple=True, default=(1, 3, 10, 30, 100, 300), help="Lane sample time (min); repeat for each lane.")
# Run
@click.option("--n-batches", type=click.IntRange(min=1), default=20, help="Independent batches, run in parallel.")
@click.option("--chains-per-batch", type=click.IntRange(min=1), default=100, help="DNA molecules per batch.")
@click.option("--processes", type=click.IntRange(min=1), default=None, help="Worker processes [default: all CPUs].")
@click.option("--seed", type=int, default=None, help="Random seed for reproducibility.")
@click.option("--style", type=click.Choice(["uv", "inverted"]), default="uv", help="Bright bands on black, or dark bands on white.")
@click.option("--out", type=click.Path(dir_okay=False, writable=True), default="gel.png", help="Output image path.")
def main(
    n_nucleosomes,
    core_length,
    linker_mean,
    linker_sd,
    nfr_prob,
    nfr_length,
    k_linker,
    k_nuc,
    k_exo_linker,
    k_exo_nuc,
    dt,
    min_len,
    times,
    n_batches,
    chains_per_batch,
    processes,
    seed,
    style,
    out,
):
    """Monte Carlo simulation of MNase digestion of chromatin, rendered as an
    agarose gel time course."""
    template_kwargs = {
        "n_nucleosomes": n_nucleosomes,
        "n_chains": chains_per_batch,
        "core_length": core_length,
        "linker_mean": linker_mean,
        "linker_sd": linker_sd,
        "nfr_prob": nfr_prob,
        "nfr_length": nfr_length,
    }
    params = DigestParams(
        k_linker=k_linker,
        k_nuc=k_nuc,
        k_exo_linker=k_exo_linker,
        k_exo_nuc=k_exo_nuc,
        dt=dt,
        min_len=min_len,
    )
    seeds = np.random.SeedSequence(seed).spawn(n_batches)

    with Pool(processes=processes) as pool:
        batches = list(
            tqdm(
                pool.imap_unordered(
                    _run_batch_star,
                    ((s, template_kwargs, params, times) for s in seeds),
                ),
                total=n_batches,
                desc="Digesting",
            )
        )

    lanes = {
        f"{t:g} min": np.concatenate([batch[t] for batch in batches])
        for t in times
    }
    render_gel(lanes, out, style=style, title="MNase digestion time course")
    click.echo(f"Saved gel to {out}")


if __name__ == "__main__":
    main()
