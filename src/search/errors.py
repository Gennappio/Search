"""Exception hierarchy. Every error raised on purpose by Search derives from SearchError."""


class SearchError(Exception):
    """Base class."""


class ConfigError(SearchError):
    """Invalid configuration or override path."""


class DataError(SearchError):
    """Malformed or inconsistent input data."""


class SeedError(DataError):
    """Seeds missing from the data or too few seeds."""
