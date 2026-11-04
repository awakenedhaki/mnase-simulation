import numpy as np
from src import render_gel
from src.gel import lane_profile, migration


# Test larger fragments migrate less and huge DNA is clamped at the limit
def test_migration_monotonic():
    distance = migration(np.array([100, 1000, 20_000, 100_000]))
    assert distance[0] > distance[1] > distance[2]
    assert distance[2] == distance[3] == 0


# Test lane signal is mass-weighted
def test_lane_profile_mass_weighted():
    profile = lane_profile(np.array([147] * 2 + [1470]), height=100, band_sigma=1)
    small = profile[70:].sum()
    assert np.isclose(profile.sum(), 1)
    assert np.isclose(small, 2 * 147 / (2 * 147 + 1470), atol=1e-3)


# Test render_gel writes an image of the expected shape
def test_render_gel(tmp_path):
    lanes = {"a": np.array([147, 294]), "b": np.array([147])}
    path = tmp_path / "gel.png"
    image = render_gel(lanes, str(path), height=300)
    assert path.exists()
    assert image.shape == (300, 3 * (40 + 14) + 14)
    assert image.max() <= 1
