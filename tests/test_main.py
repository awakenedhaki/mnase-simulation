from click.testing import CliRunner
from main import main

# Test Constants
SMALL_RUN = [
    "--n-nucleosomes", "10",
    "--chains-per-batch", "5",
    "--n-batches", "2",
    "--processes", "1",
    "--seed", "1",
]


# Test help lists the options
def test_help():
    result = CliRunner().invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "--linker-sd" in result.output
    assert "--time" in result.output


# Test a small run writes a gel
def test_small_run(tmp_path):
    out = tmp_path / "gel.png"
    args = SMALL_RUN + ["--time", "1", "--time", "10", "--out", str(out)]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    assert out.exists()


# Test invalid values are rejected
def test_invalid_option():
    result = CliRunner().invoke(main, ["--nfr-prob", "2"])
    assert result.exit_code == 2
    assert "--nfr-prob" in result.output
