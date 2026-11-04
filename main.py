import argparse
from multiprocessing import Pool

import numpy as np
from tqdm import tqdm

from src import DigestParams, build_template, render_gel, simulate


def run_batch(seed, template_kwargs, params, sample_times):
    rng = np.random.default_rng(seed)
    template = build_template(rng=rng, **template_kwargs)
    return simulate(template, params, sample_times, rng)


def main(args):
    template_kwargs = {
        "n_nucleosomes": args.n_nucleosomes,
        "n_chains": args.chains_per_batch,
        "core_length": args.core_length,
        "linker_mean": args.linker_mean,
        "linker_sd": args.linker_sd,
        "nfr_prob": args.nfr_prob,
        "nfr_length": args.nfr_length,
    }
    params = DigestParams(
        k_linker=args.k_linker,
        k_nuc=args.k_nuc,
        k_exo_linker=args.k_exo_linker,
        k_exo_nuc=args.k_exo_nuc,
        dt=args.dt,
        min_len=args.min_len,
    )
    seeds = np.random.SeedSequence(args.seed).spawn(args.n_batches)

    with Pool(processes=args.processes) as pool:
        batches = list(
            tqdm(
                pool.imap_unordered(
                    _run_batch_star,
                    ((seed, template_kwargs, params, args.times) for seed in seeds),
                ),
                total=args.n_batches,
                desc="Digesting",
            )
        )

    lanes = {
        f"{t:g} min": np.concatenate([batch[t] for batch in batches])
        for t in args.times
    }
    render_gel(lanes, args.out, style=args.style, title="MNase digestion time course")
    print(f"Saved gel to {args.out}")


def _run_batch_star(job):
    return run_batch(*job)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Monte Carlo simulation of MNase digestion of chromatin, "
        "rendered as an agarose gel time course."
    )
    template = parser.add_argument_group("chromatin template")
    template.add_argument("--n-nucleosomes", type=int, default=100, help="nucleosomes per DNA molecule")
    template.add_argument("--core-length", type=int, default=147)
    template.add_argument("--linker-mean", type=float, default=30)
    template.add_argument("--linker-sd", type=float, default=8)
    template.add_argument("--nfr-prob", type=float, default=0.0, help="probability a linker is a nucleosome-free region")
    template.add_argument("--nfr-length", type=int, default=300)

    kinetics = parser.add_argument_group("digestion kinetics")
    defaults = DigestParams()
    kinetics.add_argument("--k-linker", type=float, default=defaults.k_linker, help="endonuclease cuts per bp per min in linkers")
    kinetics.add_argument("--k-nuc", type=float, default=defaults.k_nuc, help="endonuclease cuts per bp per min in cores")
    kinetics.add_argument("--k-exo-linker", type=float, default=defaults.k_exo_linker, help="bp trimmed per end per min in linkers")
    kinetics.add_argument("--k-exo-nuc", type=float, default=defaults.k_exo_nuc, help="bp trimmed per end per min in cores")
    kinetics.add_argument("--dt", type=float, default=defaults.dt)
    kinetics.add_argument("--min-len", type=int, default=defaults.min_len)
    kinetics.add_argument("--times", type=float, nargs="+", default=[1, 3, 10, 30, 100, 300], help="lane sample times (min)")

    run = parser.add_argument_group("run")
    run.add_argument("--n-batches", type=int, default=20)
    run.add_argument("--chains-per-batch", type=int, default=100)
    run.add_argument("--processes", type=int, default=None, help="defaults to all CPUs")
    run.add_argument("--seed", type=int, default=None)
    run.add_argument("--style", choices=["uv", "inverted"], default="uv")
    run.add_argument("--out", default="gel.png")
    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
