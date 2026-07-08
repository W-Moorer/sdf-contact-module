from ..sdf_grid import TrilinearSDFGrid


class SDFCache:
    def __init__(self):
        self._cache = {}

    def load(self, path):
        if path not in self._cache:
            self._cache[path] = TrilinearSDFGrid(path)
        return self._cache[path]

    def clear(self):
        self._cache.clear()
