import re


class RMDEntity:
    def __init__(self, kind, entity_id, raw_lines, line_start):
        self.kind = kind.upper()
        self.id = entity_id
        self.raw_lines = raw_lines
        self.line_start = line_start
        self.attrs = {}
        self._parse()

    def _parse(self):
        for line in self.raw_lines:
            s = line.strip()
            if not s or s.startswith('!'):
                continue
            if s.startswith(','):
                s = s[1:].strip()
            if '=' not in s:
                continue
            eq = s.index('=')
            key = s[:eq].strip()
            val = s[eq+1:].strip()
            if key in self.attrs:
                self.attrs[key] = self.attrs[key] + ' ' + val
            else:
                self.attrs[key] = val

    def get(self, key, default=None):
        return self.attrs.get(key.upper(), default)

    def get_int(self, key, default=0):
        v = self.get(key)
        if v is None:
            return default
        try:
            return int(v)
        except (ValueError, TypeError):
            return default

    def get_float(self, key, default=0.0):
        v = self.get(key)
        if v is None:
            return default
        try:
            return float(v)
        except (ValueError, TypeError):
            return default

    def get_floats(self, key):
        v = self.get(key)
        if v is None:
            return []
        nums = []
        for part in re.split(r'[,\s]+', v.strip()):
            part = part.strip().rstrip('D')
            if not part:
                continue
            try:
                nums.append(float(part))
            except ValueError:
                pass
        return nums

    def get_bool(self, key, default=False):
        v = self.get(key)
        if v is None:
            return default
        return v.upper() in ('TRUE', '1')

    def get_str(self, key, default=''):
        v = self.get(key)
        if v is None:
            return default
        return v.strip("' ")

    def __repr__(self):
        n = len(self.raw_lines)
        return f"RMDEntity({self.kind}, id={self.id}, lines={n})"


class RMDLexer:
    def __init__(self):
        self.warnings = []

    def parse(self, text):
        lines = text.split('\n')
        return self._parse_entities(lines)

    def load(self, path):
        with open(path, encoding='utf-8', errors='replace') as f:
            text = f.read()
        return self.parse(text)

    def _parse_entities(self, lines):
        entities = []
        i = 0
        n = len(lines)

        while i < n:
            line = lines[i]
            stripped = line.strip()

            if not stripped or stripped.startswith('!'):
                i += 1
                continue

            if stripped.upper() == 'END':
                break

            match = re.match(r'(\w+)\s*/\s*(\d+)', stripped, re.IGNORECASE)
            if not match:
                m = re.match(r'(\w+)\s*/', stripped, re.IGNORECASE)
                if m:
                    kind = m.group(1)
                    ent_lines, end = self._collect(lines, i, n)
                    entities.append(RMDEntity(kind.upper(), 0, ent_lines, i))
                    i = end
                    continue
                i += 1
                continue

            kind = match.group(1).upper()
            eid = int(match.group(2))
            ent_lines, end = self._collect(lines, i, n)
            entities.append(RMDEntity(kind, eid, ent_lines, i))
            i = end

        return entities

    def _collect(self, lines, start, n):
        result = [lines[start]]
        i = start + 1
        in_patches = False

        while i < n:
            line = lines[i]
            stripped = line.strip()

            if not stripped:
                i += 1
                if not in_patches:
                    break
                continue

            if stripped.startswith('!'):
                i += 1
                continue

            if in_patches:
                if stripped.upper() == 'END':
                    break
                if re.match(r'(PART|MARKER|JOINT|GGEOM|MOTION|'
                           r'AXIAL_FORCE|EXPRESSION|GGEOMCONTACT|'
                           r'UNITS|OUTPUT|INTPAR|EQUILIBRIUM|SOLVEROPTION|'
                           r'STOPBYCONDITION|INFO|ACCGRAV)\s*/', stripped, re.IGNORECASE):
                    break
                result.append(line)
                i += 1
                continue

            if 'PATCHES' in stripped.upper():
                result.append(line)
                in_patches = True
                i += 1
                continue

            if not stripped.startswith(','):
                break

            result.append(line)
            i += 1

        return result, i
