"""Download real GWAS summary statistics for all 4 traits."""
import subprocess
from pathlib import Path
from .traits import TRAIT_CATALOG, TraitRegistry


def download_all(data_dir: str = "data/sumstats", force: bool = False):
    """Download all GWAS summary statistics."""
    out = Path(data_dir)
    out.mkdir(parents=True, exist_ok=True)

    for tid, trait in TRAIT_CATALOG.items():
        url = trait.download_url
        suffix = ".tsv.gz" if ".tsv" in url else ".txt.gz" if ".txt" in url else ".gz"
        dest = out / f"{tid}_sumstats{suffix}"

        if dest.exists() and not force:
            print(f"[{tid}] Already exists: {dest} ({dest.stat().st_size / 1e6:.1f} MB)")
            continue

        print(f"[{tid}] Downloading {trait.name} from {trait.source}...")
        print(f"  URL: {url}")
        try:
            subprocess.run(
                ["wget", "-q", "--show-progress", "--no-check-certificate",
                 "-O", str(dest), url],
                check=True, timeout=3600,
            )
            size = dest.stat().st_size / 1e6
            print(f"  ✓ {size:.1f} MB")
        except Exception as e:
            print(f"  ✗ Failed: {e}")
            if dest.exists() and dest.stat().st_size == 0:
                dest.unlink()


if __name__ == "__main__":
    download_all()
