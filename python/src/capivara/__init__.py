"""CAPIVARA experimental Python port: explicitly versioned scientific contracts."""
from .core import IFUCube, SegmentationResult, region_spectra
from .segmentation import segment_legacy

__all__ = ["IFUCube", "SegmentationResult", "region_spectra", "segment_legacy"]
__version__ = "0.0.1.dev0"
