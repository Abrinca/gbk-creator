"""Plain data holders used by the builder and the writer."""
from __future__ import annotations

from dataclasses import dataclass

from Bio.SeqFeature import CompoundLocation, SimpleLocation

Location = SimpleLocation | CompoundLocation

# Qualifier value: str, int, or None (value-less flag such as /pseudo).
Qualifiers = dict[str, list]


@dataclass
class Contig:
    name: str
    number: int
    sequence: str  # upper case
    topology: str  # "linear" | "circular"
    kind: str  # "scaffold" | "chromosome" | "plasmid"
    index: int  # 0-based add order


@dataclass
class Gene:
    """A locus-tagged gene: the writer emits a `gene` feature followed by the child feature."""
    locus_tag: str
    number: int
    contig: Contig
    location: Location
    feature_type: str  # CDS | tRNA | rRNA | ncRNA | tmRNA
    pseudo: bool
    gene_location: Location  # = location, except for in-order joins (PGAP: the span)
    gene_qualifiers: Qualifiers
    child_qualifiers: Qualifiers
    order: int  # global add order (stable tie-break)

    @property
    def start(self) -> int:
        return int(self.location.start)

    @property
    def has_protein(self) -> bool:
        return "translation" in self.child_qualifiers


@dataclass
class UntaggedFeature:
    contig: Contig
    feature_type: str
    location: Location
    qualifiers: Qualifiers
    order: int

    @property
    def start(self) -> int:
        return int(self.location.start)
