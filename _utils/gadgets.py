class namespace():
    """
    A simple namespace class to hold attributes as properties.
    """
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return f"namespace({', '.join(f'{k}={v!r}' for k, v in self.__dict__.items())})"