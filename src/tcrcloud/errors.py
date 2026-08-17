"""Package-specific exceptions for TCRcloud.

Library-style modules should raise `TCRcloudError` for expected, user-facing
failures (bad input files, unreachable AIRR repositories, empty results, ...)
instead of printing to stderr and calling ``sys.exit()``. The CLI entry point
(`tcrcloud.TCRcloud.main`) catches `TCRcloudError`, prints a clean
``TCRcloud error: ...`` message, and exits with a non-zero status, so modules
stay importable and testable as a library.
"""


class TCRcloudError(Exception):
    """An expected, user-facing error raised by TCRcloud commands.

    The message should already be complete and human-readable; the CLI adds
    the ``TCRcloud error: `` prefix when displaying it, so do not include it
    here.
    """
