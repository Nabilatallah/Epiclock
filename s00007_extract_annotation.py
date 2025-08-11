#!/usr/bin/env python3
# s00007_extract_annotation.py
# Author: Nabil Atallah

"""
Illumina 450k Annotation Extractor (SOLID/OOP/CLI)

This script downloads, extracts, and optionally processes Illumina Human Methylation 450k
annotation data, applying rpy2-based RDA loading and CSV export when available.

- Object-oriented architecture with single-responsibility components
- Structured logging to file + console (durations, step markers, emojis)
- Optional system resource monitoring (CPU/mem/disk deltas per step) via psutil
- Safe extraction with path traversal guard
- rpy2-enabled optional RDA loading and merge -> CSV export
- Config-aware (uses config.py if available), else sensible defaults
- CLI flags to customize URL/paths/behavior
- SOLID Principles implemented in this design

Usage:
  python s00007_extract_annotation.py \
    --process-rda \
    --url https://bioconductor.org/packages/release/data/annotation/src/contrib/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz

Output Files:
1) Downloaded archive (.tar.gz)
   <data_dir>/Illumina450k_annotation.tar.gz

2) Extracted package directory
   <data_dir>/s00007_<tag>_annotation/

3) Merged annotation CSV (when --process-rda and rpy2 is available)
   <data_dir>/s00007_<tag>_annotation_merged.csv

4) Structured log file
   <results_dir>/s00007_extract_annotation_<label>.txt
   where <label> is --geo_id if provided, else --tag
"""

# === Built-in Python modules ===
import os                   # For file path and directory operations
import argparse             # For parsing command-line arguments
import logging              # For logging output to file and console
import shutil               # For disk usage monitoring
import tarfile              # For extracting .tar.gz archives safely
import time                 # For timing code execution
from datetime import timedelta, datetime  # For formatting time durations
from collections import Counter           # For counting phenotype values (not used here but kept per request)
from typing import Optional, Tuple, Dict, Any

# === External packages ===
import requests             # For downloading remote files
import pandas as pd         # For data handling
import psutil               # For monitoring CPU and memory usage

# Optional rpy2 (RDA loading). If missing, --process-rda will be skipped.
_RPY2_AVAILABLE = True
try:
    from rpy2.robjects import r
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
except Exception:
    _RPY2_AVAILABLE = False

# === Project-specific config file (with path constants) ===
# Try to import everything requested; fall back to sensible defaults if not present.
_DEFAULT_BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
try:
    from config import soft_file, gse_object_file, data_dir, results_dir  # type: ignore
except Exception:
    soft_file = os.path.join(_DEFAULT_BASE, "data", "placeholder.soft")
    gse_object_file = os.path.join(_DEFAULT_BASE, "data", "placeholder_gse.pkl")
    data_dir = os.path.join(_DEFAULT_BASE, "data")
    results_dir = os.path.join(_DEFAULT_BASE, "results")

# Ensure base dirs exist early
os.makedirs(data_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)


# === Utility Function ===
def format_duration(seconds: float) -> str:
    """Format a duration in seconds as hh:mm:ss."""
    return str(timedelta(seconds=round(seconds)))


# === Logger Interface ===
class Logger:
    """Abstract base class for logging implementations."""
    def log(self, message: str) -> None:
        raise NotImplementedError


class FileLogger(Logger):
    """Logs messages using Python's logging module."""
    def log(self, message: str) -> None:
        logging.info(message)


# === System Resource Reporter ===
class SystemMonitor:
    """Reports system usage (CPU, memory, and disk)."""
    def __init__(self, logger: Logger):
        self.logger = logger

    def snapshot(self) -> Dict[str, Any]:
        # Use psutil for CPU/memory; shutil for disk usage
        cpu = psutil.cpu_percent(interval=None)
        mem = psutil.virtual_memory()
        disk = shutil.disk_usage(".")
        return {
            "cpu_percent": cpu,
            "mem_used": mem.used,
            "mem_total": mem.total,
            "disk_used": disk.used,
            "disk_total": disk.total,
        }

    def report(self) -> None:
        snap = self.snapshot()
        self.logger.log("🖥️ System Resource Usage:")
        self.logger.log(f"  ⚙️ CPU: {snap['cpu_percent']:.0f}%")
        self.logger.log(f"  🧠 Memory: {snap['mem_used'] / (1024**3):.2f}/{snap['mem_total'] / (1024**3):.2f} GB")
        self.logger.log(f"  💽 Disk: {snap['disk_used'] / (1024**3):.2f}/{snap['disk_total'] / (1024**3):.2f} GB\n")

    @staticmethod
    def delta(before: Dict[str, Any], after: Dict[str, Any]) -> Dict[str, Any]:
        return {
            "rss_delta": (after["mem_used"] - before["mem_used"]) / (1024**2),     # MB
            "disk_delta": (after["disk_used"] - before["disk_used"]) / (1024**2),  # MB
            "cpu_after": after["cpu_percent"],
        }


# === Step Timing Context ===
class StepTimer:
    """Context manager for step timing + resource delta logging."""
    def __init__(self, title: str, monitor: SystemMonitor, logger: Logger):
        self.title = title
        self.monitor = monitor
        self.logger = logger
        self.t0 = 0.0
        self.before = {}

    def __enter__(self):
        self.logger.log(f"🔍 START: {self.title}")
        self.t0 = time.time()
        self.before = self.monitor.snapshot()
        return self

    def __exit__(self, exc_type, exc, tb):
        after = self.monitor.snapshot()
        dt = format_duration(time.time() - self.t0)
        if exc:
            logging.exception(f"❌ FAIL: {self.title} (duration {dt})")
            return False  # re-raise
        delta = self.monitor.delta(self.before, after)
        self.logger.log(
            f"🔢 RSS Δ: {delta['rss_delta']:+.2f} MB | 💽 Disk Δ: {delta['disk_delta']:+.2f} MB | ⚙️ CPU: {delta['cpu_after']:.0f}%"
        )
        self.logger.log(f"✅ DONE: {self.title} in {dt}\n")


# === Core Functionality ===
class AnnotationInspector:
    """
    Encapsulates downloading, extracting, and (optionally) merging Illumina 450k RDA objects.
    Single Responsibility: handle annotation assets end-to-end.
    """
    def __init__(
        self,
        url: str,
        data_dir: str,
        results_dir: str,
        process_rda: bool,
        tag: str,
        timeout: int,
        logger: Logger,
        monitor: SystemMonitor,
    ):
        self.url = url
        self.data_dir = data_dir
        self.results_dir = results_dir
        self.process_rda = process_rda
        self.tag = tag
        self.timeout = timeout
        self.logger = logger
        self.monitor = monitor

        # Derived paths
        self.tar_path = os.path.join(self.data_dir, "Illumina450k_annotation.tar.gz")
        self.extract_dir = os.path.join(self.data_dir, f"s00007_{self.tag.lower()}_annotation")
        os.makedirs(self.extract_dir, exist_ok=True)
        self.output_csv = os.path.join(self.data_dir, f"s00007_{self.tag.lower()}_annotation_merged.csv")

        # rpy2 state
        self.rpy2_ok = _RPY2_AVAILABLE

        self.logger.log("🚀 Script started")
        self._log_startup_env()

    def _log_startup_env(self) -> None:
        self.logger.log(f"📂 data_dir: {self.data_dir}")
        self.logger.log(f"📂 results_dir: {self.results_dir}")
        self.logger.log(f"🌐 URL: {self.url}")
        self.logger.log(f"🧪 process_rda: {self.process_rda}")
        self.logger.log(f"🧩 rpy2_available: {self.rpy2_ok}")
        SystemMonitor(self.logger).report()

    # ---- Public API ----
    def run(self) -> None:
        with StepTimer("Step 1 - Download annotation archive", self.monitor, self.logger):
            self._download_if_needed(self.url, self.tar_path, self.timeout)

        with StepTimer("Step 2 - Extract .tar.gz archive", self.monitor, self.logger):
            self._safe_extract(self.tar_path, self.extract_dir)

        if self.process_rda and self.rpy2_ok:
            with StepTimer("Step 3 - Load and convert RDA files", self.monitor, self.logger):
                loc_path, other_path = self._discover_rda(self.extract_dir)
                locations_df, other_df = self._load_rda_to_dfs(loc_path, other_path)

            with StepTimer("Step 4 - Merge Locations and Other", self.monitor, self.logger):
                merged = self._merge_annotation_frames(locations_df, other_df)

            with StepTimer("Step 5 - Export annotation CSV", self.monitor, self.logger):
                merged.to_csv(self.output_csv, index=False)
                self.logger.log(f"📄 Saved: {self.output_csv}")
        elif self.process_rda and not self.rpy2_ok:
            self.logger.log("⚠️ rpy2 not available. Skipping RDA processing/CSV export.")

        self.logger.log("🎉 Finished successfully")

    # ---- Private helpers ----
    @staticmethod
    def _download_if_needed(url: str, out_path: str, timeout: int = 60) -> None:
        if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            return
        resp = requests.get(url, timeout=timeout)
        resp.raise_for_status()
        with open(out_path, "wb") as f:
            f.write(resp.content)

    @staticmethod
    def _safe_extract(tar_path: str, extract_dir: str) -> None:
        def is_within_directory(directory: str, target: str) -> bool:
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
            return os.path.commonprefix([abs_directory, abs_target]) == abs_directory

        with tarfile.open(tar_path, "r:gz") as tar:
            for member in tar.getmembers():
                member_path = os.path.join(extract_dir, member.name)
                if not is_within_directory(extract_dir, member_path):
                    raise Exception("Potential path traversal detected in tar file")
            tar.extractall(path=extract_dir)

    @staticmethod
    def _discover_rda(root: str) -> Tuple[str, str]:
        cand_loc, cand_other = None, None
        for dirpath, _, filenames in os.walk(root):
            for fn in filenames:
                low = fn.lower()
                if low == "locations.rda":
                    cand_loc = os.path.join(dirpath, fn)
                elif low == "other.rda":
                    cand_other = os.path.join(dirpath, fn)
            if cand_loc and cand_other:
                break
        if not (cand_loc and cand_other):
            raise FileNotFoundError("Could not locate Locations.rda and Other.rda in extracted archive.")
        return cand_loc, cand_other

    @staticmethod
    def _merge_annotation_frames(locations_df: pd.DataFrame, other_df: pd.DataFrame) -> pd.DataFrame:
        locations_df = locations_df.copy()
        other_df = other_df.copy()
        locations_df.index = locations_df.index.astype(str)
        other_df.index = other_df.index.astype(str)
        merged = locations_df.join(other_df, how="inner")
        merged = merged.reset_index().rename(columns={"index": "CpG"})
        return merged

    @staticmethod
    def _load_rda_to_dfs(locations_path: str, other_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        pandas2ri.activate()
        # --- Patch: precompute R-friendly paths to avoid f-string backslash expressions ---
        loc_r_path = locations_path.replace("\\", "/")
        other_r_path = other_path.replace("\\", "/")
        r(f"suppressPackageStartupMessages(load('{loc_r_path}'))")
        r(f"suppressPackageStartupMessages(load('{other_r_path}'))")
        # ----------------------------------------------------------------------------------
        r("locations_df <- as.data.frame(Locations)")
        r("other_df <- as.data.frame(Other)")
        locations_df = pandas2ri.rpy2py(r['locations_df'])
        other_df = pandas2ri.rpy2py(r['other_df'])
        return locations_df, other_df


# === CLI ===
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download, extract, and (optionally) process Illumina 450k annotation."
    )
    parser.add_argument(
        "--url",
        default="https://bioconductor.org/packages/release/data/annotation/src/contrib/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz",
        help="URL to the Illumina 450k annotation .tar.gz",
    )
    parser.add_argument(
        "--data_dir",
        default=data_dir,
        help="Directory for data artifacts (default: config.data_dir or ../data)",
    )
    parser.add_argument(
        "--results_dir",
        default=results_dir,
        help="Directory for logs/results (default: config.results_dir or ../results)",
    )
    parser.add_argument(
        "--tag",
        default="Illumina450k",
        help="Short tag used in output filenames/logs (e.g., Illumina450k)",
    )
    parser.add_argument(
        "--geo_id",
        default=None,
        help="Optional GEO ID label for logs; if provided, log file uses this instead of tag",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="HTTP download timeout in seconds",
    )
    parser.add_argument(
        "--process-rda",
        action="store_true",
        help="If set, attempts to load Locations.rda and Other.rda via rpy2 and export merged CSV",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # === Logging Setup ===
    os.makedirs(args.results_dir, exist_ok=True)
    label = args.geo_id if args.geo_id else args.tag
    log_file = os.path.join(args.results_dir, f"s00007_extract_annotation_{label}.txt")
    logging.basicConfig(filename=log_file, filemode="w", level=getattr(logging, args.log_level), format="%(message)s")
    console = logging.StreamHandler()
    console.setFormatter(logging.Formatter("%(message)s"))
    logging.getLogger().addHandler(console)

    file_logger = FileLogger()
    monitor = SystemMonitor(file_logger)

    # Initial environment report to both file and console
    file_logger.log(f"🕓 Started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    monitor.report()

    try:
        inspector = AnnotationInspector(
            url=args.url,
            data_dir=args.data_dir,
            results_dir=args.results_dir,
            process_rda=args.process_rda,
            tag=args.tag,
            timeout=args.timeout,
            logger=file_logger,
            monitor=monitor,
        )
        inspector.run()
    except Exception as e:
        logging.exception(f"💥 Unhandled error: {e}")
        raise
    finally:
        file_logger.log(f"🕓 Completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
