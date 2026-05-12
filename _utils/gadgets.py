import os, shutil, gc
from tabnanny import verbose
class namespace():
    """
    A simple namespace class to hold attributes as properties.
    """
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return f"namespace({', '.join(f'{k}={v!r}' for k, v in self.__dict__.items())})"
    
def mv_symlink(src, dst):
    '''If src is a file, move it to src and create a symlink from src to dst'''
    # check if src is a file not a symlink
    if os.path.isfile(src) and not os.path.islink(src):
        try:
            # safety for copying across fs: first copy file and then remove original
            shutil.copy2(src, dst); os.remove(src)
            os.symlink(dst, src)
        except: return False
    return True

def force_gc(func):
    """
    A decorator to force garbage collection after function execution.
    """
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        gc.collect()
        return result
    return wrapper

def check_parquet(filepath, remove_broken=True):
    """
    Check a parquet file for corruption and empty rows (memory-efficient).
    
    Reads only metadata (footer), not data rows, using O(1) memory.
    Useful for validating incomplete writes from parallel jobs.
    
    Args:
        filepath (str): Path to a single .parquet file to check.
        remove_broken (bool): If True, delete the file if broken/empty; if False, only report.
        verbose (bool): If True, print results to stdout.
    
    Returns: True if file is valid, otherwise False
    """
    try:
        from pyarrow.parquet import ParquetFile
    except ImportError:
        raise ImportError("pyarrow is required for check_parquet; install via: pip install pyarrow")
    
    if not os.path.isfile(filepath): return False
    
    try:
        # ParquetFile() reads only footer/metadata, not data — O(1) memory
        pf = ParquetFile(filepath)
        if pf.metadata.num_rows == 0:
            if remove_broken: os.remove(filepath)
            return False
        else: return True
    except Exception as e:
        if remove_broken: os.remove(filepath)
        return False