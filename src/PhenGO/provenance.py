"""Reproducibility and output-safety helpers for PhenGO."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from .constants import PhenGO_VERSION


def dependency_versions(packages=(
    "numpy", "pandas", "scikit-learn", "scipy", "networkx", "tensorflow",
)) -> dict:
    versions = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            continue
    return versions


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def describe_file(path: str | None, release: str | None = None) -> dict | None:
    if not path:
        return None
    absolute = os.path.abspath(os.fspath(path))
    info = {"path": absolute, "release": release}
    if os.path.isfile(absolute):
        stat = os.stat(absolute)
        info.update({
            "size_bytes": int(stat.st_size),
            "sha256": sha256_file(absolute),
            "modified_utc": datetime.fromtimestamp(
                stat.st_mtime, tz=timezone.utc
            ).isoformat(),
        })
    else:
        info["missing"] = True
    return info


def git_commit(repo_dir: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repo_dir, check=True,
            capture_output=True, text=True,
        )
        return result.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def git_status(repo_dir: str) -> list[str]:
    try:
        result = subprocess.run(
            ["git", "status", "--short"], cwd=repo_dir, check=True,
            capture_output=True, text=True,
        )
        return result.stdout.splitlines()
    except (OSError, subprocess.SubprocessError):
        return []


def source_tree_sha256(repo_dir: str) -> str | None:
    """Hash Python source content so dirty or unpackaged code is identifiable."""
    source_dir = Path(repo_dir) / "src" / "PhenGO"
    if not source_dir.is_dir():
        return None
    digest = hashlib.sha256()
    for path in sorted(source_dir.rglob("*.py")):
        digest.update(str(path.relative_to(source_dir)).encode("utf-8"))
        digest.update(b"\0")
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def assert_safe_output_dir(path: str, protected_paths=()) -> str:
    """Reject broad or input-containing output directories before deletion."""
    resolved = Path(path).expanduser().resolve()
    forbidden = {Path("/").resolve(), Path.home().resolve(), Path.cwd().resolve()}
    if resolved in forbidden:
        raise ValueError(f"Refusing to use protected output directory: {resolved}")

    for protected in protected_paths:
        if not protected:
            continue
        protected_resolved = Path(protected).expanduser().resolve()
        if resolved == protected_resolved or resolved in protected_resolved.parents:
            raise ValueError(
                f"Output directory {resolved} contains or equals input path "
                f"{protected_resolved}; refusing destructive overwrite."
            )
    return str(resolved)


def prepare_output_dir(path: str, overwrite: bool, protected_paths=(),
                       preserve_suffixes=()) -> str:
    path = assert_safe_output_dir(path, protected_paths)
    entries = os.listdir(path) if os.path.isdir(path) else []
    removable_entries = [
        entry for entry in entries
        if not (
            os.path.isfile(os.path.join(path, entry)) and
            any(entry.endswith(suffix) for suffix in preserve_suffixes)
        )
    ]
    if removable_entries:
        if not overwrite:
            raise ValueError(
                f"Output directory '{path}' is not empty. Choose a new directory "
                "or pass -overwrite."
            )
        for entry in removable_entries:
            target = os.path.join(path, entry)
            if os.path.isdir(target) and not os.path.islink(target):
                shutil.rmtree(target)
            else:
                os.unlink(target)
    else:
        os.makedirs(path, exist_ok=True)
    return path


def build_run_manifest(*, repo_dir: str, species: str, snapshot_id: str | None,
                       strict_snapshot: bool, inputs: dict, releases: dict,
                       policies: dict, options: dict, outputs: dict,
                       counts: dict, resource_context: dict | None = None) -> dict:
    status = git_status(repo_dir)
    return {
        "schema_version": 3,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "tool": "PhenGO",
        "tool_version": PhenGO_VERSION,
        "git_commit": git_commit(repo_dir),
        "git_dirty": bool(status),
        "git_status": status,
        "source_tree_sha256": source_tree_sha256(repo_dir),
        "python": {
            "version": sys.version,
            "executable": sys.executable,
            "platform": platform.platform(),
            "machine": platform.machine(),
        },
        "command": list(sys.argv),
        "dependencies": dependency_versions(),
        "species": species,
        "snapshot_id": snapshot_id,
        "strict_snapshot": bool(strict_snapshot),
        "releases": releases,
        "resource_context": resource_context or {},
        "inputs": inputs,
        "policies": policies,
        "options": options,
        "counts": counts,
        "outputs": outputs,
    }


def write_manifest(path: str, manifest: dict) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


def load_manifest_for_arff(arff_path: str) -> dict | None:
    directory = os.path.dirname(os.path.abspath(arff_path))
    candidates = [
        os.path.join(directory, "PhenGO_manifest.json"),
        f"{arff_path}.manifest.json",
    ]
    for candidate in candidates:
        if os.path.isfile(candidate):
            with open(candidate, encoding="utf-8") as handle:
                return json.load(handle)
    return None
