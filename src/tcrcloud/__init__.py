from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("TCRcloud")
except PackageNotFoundError:
    # Package isn't installed (e.g. running directly from a source checkout).
    __version__ = "0.0.0.dev"

__all__ = ["__version__"]
