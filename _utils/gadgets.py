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
    import os
    import shutil

    # check if src is a file not a symlink
    if os.path.isfile(src) and not os.path.islink(src):
        try:
            # safety for copying across fs: first copy file and then remove original
            shutil.copy2(src, dst); os.remove(src)
            os.symlink(dst, src)
        except: return False
    return True