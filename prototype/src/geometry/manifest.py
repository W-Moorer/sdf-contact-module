import os


class GeometryManifest:
    def __init__(self):
        self.entries = {}

    @classmethod
    def load(cls, path):
        m = cls()
        if not os.path.exists(path):
            return m
        if path.endswith('.yaml') or path.endswith('.yml'):
            import yaml
            with open(path) as f:
                data = yaml.safe_load(f)
            for body_name, info in data.get('bodies', {}).items():
                m.entries[body_name] = {
                    'mesh': info.get('mesh'),
                    'sdf': info.get('sdf'),
                    'surface_role': info.get('surface_role', 'quadrature_and_sdf'),
                    'body_id': info.get('body_id'),
                }
        return m

    def add(self, body_name, mesh_path=None, sdf_path=None, surface_role='quadrature_and_sdf', body_id=None):
        self.entries[body_name] = {
            'mesh': mesh_path,
            'sdf': sdf_path,
            'surface_role': surface_role,
            'body_id': body_id,
        }

    def get(self, body_name):
        return self.entries.get(body_name)

    def resolve_path(self, body_name, key, base_dir='.'):
        entry = self.entries.get(body_name)
        if entry is None:
            return None
        path = entry.get(key)
        if path is None:
            return None
        full = os.path.join(base_dir, path) if not os.path.isabs(path) else path
        return full if os.path.exists(full) else None
