"""
pytest configuration, lives at repo root so pytest auto-discovers it

adds scripts/ to sys.path for test session so test modules can import
stage_registry directly without packaging project or hardcoding paths
inside each test file
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "scripts"))
