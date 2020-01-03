import gzip


def open_maybe_gz(path, mode: str, **kwargs):
    """
    Like the open() built-in but transparently handles .gz files.

    Can be used as a context manager.
    """
    # Since `open` and `gzip.open` disagree on the default mode and whether just
    # 'r' implies 'b' or 't', the caller must be unambiguously specify the mode.
    # Plus, we don't support any of the write modes.
    if mode in ('rb', 'rt'):
        with open(path, 'rb') as f:
            header = f.read(2)
        if header == b'\x1f\x8b':  # GZIP magic cookie
            return gzip.open(path, mode, **kwargs)
        else:
            return open(path, mode, **kwargs)
    else:
        raise ValueError("Unsupported mode (must be 'rb' or 'rt'):", mode)
