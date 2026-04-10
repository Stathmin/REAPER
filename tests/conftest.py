import os
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
if str(root) not in sys.path:
    sys.path.insert(0, str(root))


def pytest_configure(config):
    """Ensure tests run from repo root (fixes cwd-deleted failures from temp_dir chdir)."""
    os.chdir(root)


def pytest_runtest_setup(item):
    """Reset cwd before each test (previous test may have chdir'd into a deleted temp dir)."""
    os.chdir(root)

