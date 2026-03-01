"""Base class for all prioritization methods."""
from abc import ABC, abstractmethod
import pandas as pd


class BaseMethod(ABC):
    """Interface all prioritization methods must implement."""

    @property
    @abstractmethod
    def name(self) -> str: ...

    @abstractmethod
    def score(self, locus_data: pd.DataFrame) -> pd.Series:
        """Return gene-level scores (higher = more likely causal)."""
        ...

    def __repr__(self):
        return f"<{self.__class__.__name__}: {self.name}>"
