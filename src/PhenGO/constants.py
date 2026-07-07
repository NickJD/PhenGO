from __future__ import annotations

import logging
import os
import sys
from importlib import resources

PhenGO_VERSION = "v0.2.0"


def _resource_path(package: str, filename: str) -> str:
	"""Return a path-like string for a bundled package data file."""
	return str(resources.files(package).joinpath(filename))

# WORM
DEFAULT_WORM_LETHAL_PHENOTYPES_FILE = _resource_path(
	"PhenGO.data.worm", "lethal_terms_traversed_2025-08-12.tsv.gz")
DEFAULT_WORM_LETHAL_GENES_FILE = _resource_path(
	"PhenGO.data.worm", "genes_direct_and_inferred_for_WBPhenotype_0000062_11-08-2025.txt.gz")

# MOUSE
DEFAULT_MOUSE_PHENOTYPES_FILE = _resource_path(
	"PhenGO.data.mouse", "mouse_lethal_terms.txt.gz")

# FLY
DEFAULT_FLY_LETHAL_GENES_FILE = _resource_path(
	"PhenGO.data.fly", "FlyBase_Lethal_Gene_IDs_2025-08-15.txt.gz")
DEFAULT_FLY_SPECIES_FIELDS_FILE = _resource_path(
	"PhenGO.data.fly", "FlyBase_Fields_2025_07_29.tsv.gz")
DEFAULT_FLY_HELPER_LINES_FILE = _resource_path(
	"PhenGO.data.fly", "FlyBase_Helper_Lines_2025-08-15.txt.gz")


def configure_logger(name: str,
					 enable_file: bool = True,
					 log_dir: str | None = None,
					 logfile_name: str | None = None,
					 level=logging.INFO,
					 file_level=logging.DEBUG):
	"""Configure and return a logger with a stdout StreamHandler and optional FileHandler.

	- name: logger name (e.g., 'PhenGO' or 'PhenGO.GO_Compare')
	- enable_file: whether to create a file handler
	- log_dir: directory to place logfile (created if missing). If None, uses current working dir.
	- logfile_name: filename for logfile; if None, defaults to '<name>.log' with dots replaced.
	- level: logging level for the stream handler
	- file_level: logging level for the file handler
	"""
	console_fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	file_fmt    = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

	logger = logging.getLogger(name)
	logger.setLevel(level)
	logger.propagate = False

	if not logger.handlers:
		sh = logging.StreamHandler(stream=sys.stdout)
		sh.setLevel(level)
		sh.setFormatter(console_fmt)
		logger.addHandler(sh)

	root = logging.getLogger()
	root.setLevel(level)
	if not any(isinstance(h, logging.StreamHandler) and
			   not isinstance(h, logging.FileHandler) for h in root.handlers):
		rsh = logging.StreamHandler(stream=sys.stdout)
		rsh.setLevel(level)
		rsh.setFormatter(console_fmt)
		root.addHandler(rsh)

	if enable_file:
		if not log_dir:
			log_dir = os.getcwd()
		os.makedirs(log_dir, exist_ok=True)
		if not logfile_name:
			logfile_name = f"{name.replace('.', '_')}.log"
		fh_path = os.path.join(log_dir, logfile_name)

		# Named logger file handler
		if not any(isinstance(h, logging.FileHandler) for h in logger.handlers):
			fh = logging.FileHandler(fh_path)
			fh.setLevel(file_level)
			fh.setFormatter(file_fmt)
			logger.addHandler(fh)

		# Root file handler so module messages also land in the log file
		if not any(isinstance(h, logging.FileHandler) for h in root.handlers):
			rfh = logging.FileHandler(fh_path)
			rfh.setLevel(file_level)
			rfh.setFormatter(file_fmt)
			root.addHandler(rfh)

	return logger
