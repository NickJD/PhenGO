"""Importable console entry point for PhenGO-Predict."""
from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_legacy_cli():
    script_path = Path(__file__).with_name("PhenGO-Predict.py")
    spec = importlib.util.spec_from_file_location("PhenGO.predict._legacy_cli", script_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load PhenGO-Predict CLI from {script_path}")
    module = importlib.util.module_from_spec(spec)
    module.__package__ = "PhenGO.predict"
    spec.loader.exec_module(module)
    return module


def main():
    return _load_legacy_cli().main()


__all__ = ["main"]


if __name__ == "__main__":
    main()
