"""Knowledge memory: experiment registry and failure archive."""
import json
from pathlib import Path
from dataclasses import dataclass, field, asdict
from datetime import datetime


@dataclass
class Experiment:
    name: str
    version: str
    timestamp: str = ""
    auroc: float = 0.0
    promoted: bool = False
    notes: str = ""

    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = datetime.utcnow().isoformat()


class ExperimentRegistry:
    """Persistent experiment tracking to prevent rediscovery of failed ideas."""

    def __init__(self, path: str | Path = "experiments.json"):
        self.path = Path(path)
        self.experiments: list[Experiment] = []
        if self.path.exists():
            data = json.loads(self.path.read_text())
            self.experiments = [Experiment(**e) for e in data]

    def log(self, experiment: Experiment):
        self.experiments.append(experiment)
        self._save()

    def was_tried(self, name: str) -> bool:
        return any(e.name == name for e in self.experiments)

    def failed_strategies(self) -> list[str]:
        return [e.name for e in self.experiments if not e.promoted]

    def _save(self):
        self.path.write_text(json.dumps([asdict(e) for e in self.experiments], indent=2))
