"""Download and validate GWAS summary statistics."""
import subprocess
import hashlib
from pathlib import Path
from .traits import TraitMetadata


class SumstatDownloader:
    """Download GWAS summary statistics with validation."""

    def __init__(self, data_dir: str | Path = "data/sumstats"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)

    def _filename(self, trait: TraitMetadata) -> Path:
        url = trait.download_url
        suffix = ".tsv.gz" if ".tsv" in url else ".txt.gz" if ".txt" in url else ".gz"
        return self.data_dir / f"{trait.trait_id}_sumstats{suffix}"

    def download(self, trait: TraitMetadata, force: bool = False) -> Path:
        dest = self._filename(trait)
        if dest.exists() and not force:
            print(f"[{trait.trait_id}] Already downloaded: {dest}")
            return dest
        print(f"[{trait.trait_id}] Downloading from {trait.source}...")
        subprocess.run(
            ["wget", "-q", "--no-check-certificate", "-O", str(dest), trait.download_url],
            check=True, timeout=1800,
        )
        size_mb = dest.stat().st_size / 1e6
        print(f"[{trait.trait_id}] Downloaded: {size_mb:.1f} MB")
        return dest

    def download_all(self, traits: list[TraitMetadata], force: bool = False) -> dict[str, Path]:
        results = {}
        for t in traits:
            try:
                results[t.trait_id] = self.download(t, force=force)
            except Exception as e:
                print(f"[{t.trait_id}] FAILED: {e}")
        return results

    def validate(self, path: Path, min_size_mb: float = 1.0) -> bool:
        if not path.exists():
            return False
        return path.stat().st_size > min_size_mb * 1e6
