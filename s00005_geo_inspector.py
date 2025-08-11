#!/usr/bin/env python3
# s00005_geo_inspector.py
# Author: Nabil Atallah
# Contact: n.atallah@northeastern.edu

"""
GEO Dataset Downloader & Inspector (SOLID/OOP/CLI)

This script downloads, parses, and inspects an NCBI GEO dataset (default: GSE40279),
saving a cached parsed object and phenotype summary for downstream analysis.

- Object-oriented architecture with single-responsibility components
- Structured logging to file + console (durations, step markers, emojis)
- Optional system resource monitoring (CPU/mem/disk deltas per step) via psutil
- Dual parsing strategy for GEOparse (filename-GEOID + explicit GEO fallback)
- Config-aware (uses config.py if available), else sensible defaults
- CLI flags to customize GEO ID, sample index, and monitoring behavior
- SOLID Principles implemented in this design

Usage:
  python s00005_geo_inspector.py \
    --geo_id GSE40279 \
    --sample_index 0

Output Files:
1) Downloaded SOFT archive (.soft.gz)
   <data_dir>/s00005_<GEO_ID>_family.soft.gz

2) Parsed GSE object cache
   <data_dir>/s00005_<GEO_ID>_gse_object.joblib

3) Parsed phenotype CSV (auto-parsed from 'characteristics_ch1')
   <results_dir>/s00005_<GEO_ID>_phenotypes.csv

4) Structured log file
   <results_dir>/s00005_download_dataset_<GEO_ID>.txt
"""

# === Built-in Python modules ===
import os                   # For file path and directory operations
import argparse             # For parsing command-line arguments
import logging              # For logging output to file and console
import shutil               # For disk usage and file management
import time                 # For timing code execution
from datetime import timedelta, datetime  # For formatting time durations & timestamps
from collections import Counter           # For counting phenotype values
from typing import Optional, Dict, Any

# === External packages ===
import GEOparse             # For downloading and parsing GEO datasets
import joblib               # For caching/loading parsed GSE objects
import pandas as pd         # For phenotype CSV export
import psutil               # For monitoring CPU and memory usage

# === Project-specific config file (with path constants) ===
# Try to import; fall back to sensible defaults if not present.
_DEFAULT_BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
try:
    from config import data_dir, results_dir  # type: ignore
except Exception:
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
        self.before: Dict[str, Any] = {}

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
class GeoDatasetFetcher:
    """
    Handles download, safe parsing, inspection, phenotype summarization, and verification
    for a given GEO dataset.
    """
    def __init__(self, geo_id: str, logger: Logger):
        self.geo_id = geo_id
        self.logger = logger

        # Output paths (config-aware)
        self.soft_path = os.path.join(data_dir, f"s00005_{self.geo_id}_family.soft.gz")
        self.gse_object_path = os.path.join(data_dir, f"s00005_{self.geo_id}_gse_object.joblib")
        self.pheno_csv = os.path.join(results_dir, f"s00005_{self.geo_id}_phenotypes.csv")
        self.log_file = os.path.join(results_dir, f"s00005_download_dataset_{self.geo_id}.txt")

        # Suppress very verbose GEOparse DEBUG logs unless needed
        logging.getLogger("GEOparse").setLevel(logging.WARNING)

    # ---- Download ----
    def download_soft(self) -> None:
        """Download the GEO SOFT file if missing, standardize naming to soft_path."""
        if os.path.exists(self.soft_path) and os.path.getsize(self.soft_path) > 0:
            self.logger.log(f"✅ SOFT already present: {self.soft_path}")
            return

        self.logger.log(f"📥 Downloading {self.geo_id} ➜ {self.soft_path}")
        GEOparse.get_GEO(geo=self.geo_id, destdir=data_dir, how="full")
        # GEOparse writes default filenames; move/rename to our standard path if necessary
        for fn in os.listdir(data_dir):
            if fn.startswith(self.geo_id) and fn.endswith(".soft.gz"):
                src = os.path.join(data_dir, fn)
                if src != self.soft_path:
                    shutil.move(src, self.soft_path)
                break
        self.logger.log("✅ Download complete.")

    # ---- Parse or load cached GSE ----
    def load_or_parse_gse(self):
        """Load cached GSE object; else parse from SOFT using dual-approach naming/explicit GEO."""
        try:
            self.logger.log("📦 Attempting to load cached GSE object...")
            gse = joblib.load(self.gse_object_path)
            assert hasattr(gse, "gsms") and len(gse.gsms) > 0
            self.logger.log(f"⚡ Using cached GSE object — {len(gse.gsms)} samples.")
            return gse
        except Exception as e:
            self.logger.log(f"⚠️ Failed to load GSE object, will parse from SOFT... ({e})")

        # Option 2 (preferred): GEOparse expects filename starting with GEO ID
        dir_name = os.path.dirname(self.soft_path)
        friendly_name = f"{self.geo_id}_family.soft.gz"
        friendly_path = os.path.join(dir_name, friendly_name)
        if not os.path.exists(friendly_path):
            shutil.copy(self.soft_path, friendly_path)
            self.logger.log(f"📄 Created GEOparse-friendly copy: {friendly_path}")

        try:
            gse = GEOparse.get_GEO(filepath=friendly_path)
        except ValueError as ve:
            self.logger.log(f"🧯 Filename-based parse failed, trying explicit GEO: {ve}")
            # Option 1 fallback: Explicit GEO param
            gse = GEOparse.get_GEO(filepath=self.soft_path, GEO=self.geo_id)

        joblib.dump(gse, self.gse_object_path, compress=("zlib", 3))
        self.logger.log(f"💾 Re-saved parsed GSE object to {self.gse_object_path}")
        return gse

    # ---- Inspect one sample ----
    def inspect_sample(self, gse, sample_index: int = 0) -> None:
        """Log a quick view of one sample's first few metadata fields."""
        sample_id = list(gse.gsms.keys())[sample_index]
        gsm = gse.gsms[sample_id]
        self.logger.log(f"🔬 Sample {sample_index}: {sample_id}")
        for key in list(gsm.metadata.keys())[:5]:
            self.logger.log(f"  {key}: {gsm.metadata[key]}")
        self.logger.log("")

    # ---- Summarize phenotypes + CSV export ----
    def summarize_phenotypes(self, gse, field: str = "characteristics_ch1") -> None:
        """Summarize phenotype labels and export an auto-parsed CSV (one row per sample)."""
        values = []
        rows = []
        for sample_id, gsm in gse.gsms.items():
            record = {"sample_id": sample_id}
            meta = gsm.metadata.get(field)
            if meta:
                values.extend(meta)
                for item in meta:
                    if ":" in item:
                        k, v = item.split(":", 1)
                        record[k.strip().lower()] = v.strip()
            rows.append(record)

        # counts (raw)
        counter = Counter(values)
        self.logger.log(f"📊 Metadata summary: {field}")
        for label, count in counter.items():
            self.logger.log(f"  {label}: {count}")
        self.logger.log("")

        # export parsed phenotype table
        df = pd.DataFrame(rows)
        df.to_csv(self.pheno_csv, index=False)
        self.logger.log(f"📄 Parsed phenotype CSV saved to: {self.pheno_csv}\n")

    # ---- Verify cached object ----
    def verify_gse_object(self) -> None:
        """Verify that the cached object exists and can be loaded."""
        self.logger.log("🔍 Final verification of parsed GSE object")
        if os.path.exists(self.gse_object_path):
            self.logger.log(f"✅ File exists: {self.gse_object_path}")
            try:
                gse = joblib.load(self.gse_object_path)
                self.logger.log("✅ GSE object loaded.")
                self.logger.log(f"🔢 Samples loaded: {len(gse.gsms)}")
                self.logger.log(f"📑 Metadata keys: {list(gse.metadata.keys())[:5]}")
            except Exception as e:
                self.logger.log(f"❌ Failed to load GSE object: {e}")
        else:
            self.logger.log("❌ Parsed GSE object not found.")
        self.logger.log("✅ Dataset inspection complete.\n")


# === CLI / Runner ===
def configure_logging(log_file: str) -> None:
    """Configure logging to file + console, simple line-based format for easy grep/parse."""
    logging.basicConfig(
        filename=log_file,
        filemode="w",
        level=logging.INFO,
        format="%(message)s"
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter("%(message)s"))
    logging.getLogger().addHandler(console)


def main():
    parser = argparse.ArgumentParser(description="GEO Dataset Downloader & Inspector")
    parser.add_argument("--geo_id", type=str, default="GSE40279", help="GEO accession ID (e.g., GSE40279)")
    parser.add_argument("--sample_index", type=int, default=0, help="Sample index to inspect")
    parser.add_argument("--no-monitor", action="store_true", help="Disable system resource monitoring logs")
    args = parser.parse_args()

    # Configure logging with dataset-specific filename
    log_path = os.path.join(results_dir, f"s00005_download_dataset_{args.geo_id}.txt")
    configure_logging(log_path)
    logger = FileLogger()

    # Optional monitor (enabled by default)
    monitor = SystemMonitor(logger)
    logger.log(f"🕒 Run started: {datetime.now().isoformat(timespec='seconds')}")
    logger.log(f"📝 Log file: {log_path}\n")
    if not args.no_monitor:
        monitor.report()

    fetcher = GeoDatasetFetcher(args.geo_id, logger)

    # Execute steps with timing + resource deltas
    try:
        with StepTimer("Download SOFT file", monitor, logger):
            fetcher.download_soft()

        with StepTimer("Parse or load GSE object", monitor, logger):
            gse = fetcher.load_or_parse_gse()

        with StepTimer("Inspect one sample", monitor, logger):
            fetcher.inspect_sample(gse, sample_index=args.sample_index)

        with StepTimer("Summarize phenotype field", monitor, logger):
            fetcher.summarize_phenotypes(gse)

        with StepTimer("Verify GSE object", monitor, logger):
            fetcher.verify_gse_object()

    except Exception as e:
        logger.log(f"❌ Execution failed: {e}")
        raise

    logger.log("🏁 Completed all steps.")


if __name__ == "__main__":
    main()
