from typing import Dict, Literal, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

LADDER_100BP = (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000)


def migration(
    lengths: np.ndarray, min_bp: float = 80, max_bp: float = 20_000
) -> np.ndarray:
    """Relative migration distance from the well (0) to the gel bottom (1).

    Mobility is log-linear in fragment size. DNA larger than max_bp moves at the
    limiting mobility and piles up in a compression zone just below the well;
    DNA smaller than min_bp runs off the bottom.

    Args:
        lengths (np.ndarray): Fragment lengths (bp).
        min_bp (float): Size that reaches the bottom of the gel.
        max_bp (float): Size at the limiting mobility.

    Returns:
        np.ndarray: Migration distance for each fragment.
    """
    span = np.log10(max_bp) - np.log10(min_bp)
    distance = (np.log10(max_bp) - np.log10(np.maximum(lengths, 1))) / span
    return np.clip(distance, 0.0, None)


def lane_profile(
    lengths: np.ndarray,
    weights: Optional[np.ndarray] = None,
    height: int = 600,
    band_sigma: float = 2.5,
    **migration_kwargs,
) -> np.ndarray:
    """Builds the 1D intensity profile of a lane.

    By default each fragment contributes in proportion to its mass (bp), as with
    intercalating stains such as EtBr or SYBR. The histogram is then blurred
    with a Gaussian to give bands a finite width.

    Args:
        lengths (np.ndarray): Fragment lengths (bp).
        weights (Optional[np.ndarray]): Signal per fragment; defaults to its length.
        height (int): Number of pixels from well to gel bottom.
        band_sigma (float): Band width (pixels).
        **migration_kwargs: Passed to migration().

    Returns:
        np.ndarray: Intensity per pixel row, normalised to unit total loading.
    """
    lengths = np.asarray(lengths)
    weights = lengths.astype(float) if weights is None else np.asarray(weights, float)
    distance = migration(lengths, **migration_kwargs)
    profile, _ = np.histogram(
        distance, bins=height, range=(0.0, 1.0), weights=weights
    )

    radius = int(np.ceil(4 * band_sigma))
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (x / band_sigma) ** 2)
    profile = np.convolve(profile, kernel / kernel.sum(), mode="same")

    total = weights.sum()
    return profile / total if total > 0 else profile


def gel_image(
    profiles: Sequence[np.ndarray],
    lane_width: int = 40,
    gap: int = 14,
    gamma: float = 0.6,
) -> np.ndarray:
    """Lays lane profiles side by side into a 2D gel image scaled to [0, 1].

    Args:
        profiles (Sequence[np.ndarray]): Intensity profile of each lane.
        lane_width (int): Lane width (pixels).
        gap (int): Space between lanes (pixels).
        gamma (float): Display gamma; < 1 brings out faint bands, like a long exposure.

    Returns:
        np.ndarray: Image of shape (height, n_lanes * (lane_width + gap) + gap).
    """
    height = len(profiles[0])
    width = len(profiles) * (lane_width + gap) + gap

    # Soft lane edges, as from a rounded well
    x = np.linspace(-1, 1, lane_width)
    edge = np.clip(1.5 * (1 - x**8), 0, 1)

    image = np.zeros((height, width))
    for i, profile in enumerate(profiles):
        left = gap + i * (lane_width + gap)
        image[:, left : left + lane_width] = np.outer(profile, edge)

    image /= np.percentile(image[image > 0], 99.5) if np.any(image > 0) else 1
    return np.clip(image, 0, 1) ** gamma


def render_gel(
    lanes: Dict[str, np.ndarray],
    path: str,
    ladder: Sequence[int] = LADDER_100BP,
    style: Literal["uv", "inverted"] = "uv",
    height: int = 600,
    lane_width: int = 40,
    gap: int = 14,
    title: Optional[str] = None,
    **migration_kwargs,
) -> np.ndarray:
    """Renders an agarose gel of the given lanes next to a DNA ladder and saves it.

    Args:
        lanes (Dict[str, np.ndarray]): Fragment lengths per lane, keyed by label.
        path (str): Output image path.
        ladder (Sequence[int]): Ladder band sizes (bp), loaded at equal mass.
        style (Literal["uv", "inverted"]): "uv" gives bright bands on black;
            "inverted" gives dark bands on white.
        height (int): Number of pixels from well to gel bottom.
        lane_width (int): Lane width (pixels).
        gap (int): Space between lanes (pixels).
        title (Optional[str]): Figure title.
        **migration_kwargs: Passed to migration().

    Returns:
        np.ndarray: The gel image (without axes).
    """
    ladder = np.asarray(ladder)
    profiles = [
        lane_profile(ladder, np.ones(len(ladder)), height, **migration_kwargs)
    ] + [lane_profile(lengths, None, height, **migration_kwargs) for lengths in lanes.values()]
    image = gel_image(profiles, lane_width, gap)

    labels = ["Ladder"] + list(lanes.keys())
    fig, ax = plt.subplots(figsize=(1.0 + 0.55 * len(labels), 7))
    ax.imshow(
        image,
        cmap="gray" if style == "uv" else "gray_r",
        aspect="auto",
        interpolation="bilinear",
        vmin=0,
        vmax=1,
    )
    ax.set_xticks([gap + i * (lane_width + gap) + lane_width / 2 for i in range(len(labels))])
    ax.set_xticklabels(labels, rotation=45, ha="left", rotation_mode="anchor")
    ax.xaxis.tick_top()
    ax.set_yticks(migration(ladder, **migration_kwargs) * height)
    ax.set_yticklabels([f"{bp:,} bp" for bp in ladder], fontsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if title:
        ax.set_title(title, pad=40)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return image
