import numpy as np
import struct

class AnalyticPlaneSDF:
    def __init__(self, z_top=0.0):
        self.z_top = z_top

    def query(self, y_local):
        g = y_local[2] - self.z_top
        raw_grad = np.array([0.0, 0.0, 1.0])
        grad_norm = 1.0
        unit_normal = np.array([0.0, 0.0, 1.0])
        valid = True
        return g, raw_grad, grad_norm, unit_normal, valid

    def query_batch(self, Y_local):
        g = Y_local[:, 2] - self.z_top
        raw_grad = np.tile([0.0, 0.0, 1.0], (len(Y_local), 1))
        grad_norm = np.ones(len(Y_local))
        unit_normal = np.tile([0.0, 0.0, 1.0], (len(Y_local), 1))
        valid = np.ones(len(Y_local), dtype=bool)
        return g, raw_grad, grad_norm, unit_normal, valid


class TrilinearSDFGrid:
    """Reads a dense .sdf binary file and provides trilinear queries."""

    def __init__(self, path):
        self.path = path
        self._load()

    def _load(self):
        with open(self.path, 'rb') as f:
            magic = f.read(4)
            if magic != b'SDFO':
                raise ValueError(f"Not a valid .sdf file (magic={magic})")

            version = struct.unpack('<I', f.read(4))[0]
            nx = struct.unpack('<I', f.read(4))[0]
            ny = struct.unpack('<I', f.read(4))[0]
            nz = struct.unpack('<I', f.read(4))[0]

            bmin = np.array(struct.unpack('<3d', f.read(24)))
            bmax = np.array(struct.unpack('<3d', f.read(24)))
            voxel = np.array(struct.unpack('<3d', f.read(24)))
            flags = struct.unpack('<I', f.read(4))
            reserved = f.read(32)

        self.nx, self.ny, self.nz = nx, ny, nz
        self.bmin, self.bmax, self.voxel = bmin, bmax, voxel
        self.flags = flags
        self.has_hessian = bool(flags[0] & 1)
        self.is_sparse = bool(flags[0] & 2)

        if self.is_sparse:
            raise NotImplementedError("Sparse .sdf not supported yet")

        N = nx * ny * nz

        with open(self.path, 'rb') as f:
            f.read(128)
            phi_bytes = N * 8
            self.phi0 = np.frombuffer(f.read(phi_bytes), dtype=np.float64).copy()
            self.phi0 = self.phi0.reshape((nz, ny, nx)).transpose((2, 1, 0))

    def _continuous_index(self, p):
        f = (p - self.bmin) / self.voxel - 0.5
        return f

    def _clamp(self, i, j, k):
        i = np.clip(i, 0, self.nx - 1)
        j = np.clip(j, 0, self.ny - 1)
        k = np.clip(k, 0, self.nz - 1)
        return int(i), int(j), int(k)

    def query(self, y_local):
        f = self._continuous_index(y_local)

        i0 = int(np.floor(f[0]))
        j0 = int(np.floor(f[1]))
        k0 = int(np.floor(f[2]))
        u = f[0] - i0
        v = f[1] - j0
        w = f[2] - k0

        if (i0 < 0 or i0 + 1 >= self.nx or
            j0 < 0 or j0 + 1 >= self.ny or
            k0 < 0 or k0 + 1 >= self.nz):
            return 0.0, np.zeros(3), 0.0, np.array([0.0, 0.0, 1.0]), False

        phi = self.phi0
        h = self.voxel

        p000 = phi[i0,   j0,   k0]
        p100 = phi[i0+1, j0,   k0]
        p010 = phi[i0,   j0+1, k0]
        p110 = phi[i0+1, j0+1, k0]
        p001 = phi[i0,   j0,   k0+1]
        p101 = phi[i0+1, j0,   k0+1]
        p011 = phi[i0,   j0+1, k0+1]
        p111 = phi[i0+1, j0+1, k0+1]

        g = (p000 * (1-u)*(1-v)*(1-w) +
             p100 * u    *(1-v)*(1-w) +
             p010 * (1-u)*v    *(1-w) +
             p110 * u    *v    *(1-w) +
             p001 * (1-u)*(1-v)*w     +
             p101 * u    *(1-v)*w     +
             p011 * (1-u)*v    *w     +
             p111 * u    *v    *w)

        inv_hx, inv_hy, inv_hz = 1.0 / h[0], 1.0 / h[1], 1.0 / h[2]

        dg_dx = ((p100 - p000) * (1-v)*(1-w) +
                 (p110 - p010) * v    *(1-w) +
                 (p101 - p001) * (1-v)*w     +
                 (p111 - p011) * v    *w) * inv_hx

        dg_dy = ((p010 - p000) * (1-u)*(1-w) +
                 (p110 - p100) * u    *(1-w) +
                 (p011 - p001) * (1-u)*w     +
                 (p111 - p101) * u    *w) * inv_hy

        dg_dz = ((p001 - p000) * (1-u)*(1-v) +
                 (p101 - p100) * u    *(1-v) +
                 (p011 - p010) * (1-u)*v     +
                 (p111 - p110) * u    *v) * inv_hz

        raw_grad = np.array([dg_dx, dg_dy, dg_dz])
        grad_norm = np.linalg.norm(raw_grad)
        if grad_norm > 1e-30:
            unit_normal = raw_grad / grad_norm
        else:
            unit_normal = np.array([0.0, 0.0, 1.0])

        return g, raw_grad, grad_norm, unit_normal, True

    def query_batch(self, Y_local):
        Y = np.asarray(Y_local, dtype=np.float64)
        N = len(Y)
        f = (Y - self.bmin[None, :]) / self.voxel[None, :] - 0.5
        i0 = np.floor(f[:, 0]).astype(np.int32)
        j0 = np.floor(f[:, 1]).astype(np.int32)
        k0 = np.floor(f[:, 2]).astype(np.int32)
        u = f[:, 0] - i0
        v = f[:, 1] - j0
        w = f[:, 2] - k0

        # Bounds check
        valid = (i0 >= 0) & (i0 + 1 < self.nx) & (j0 >= 0) & (j0 + 1 < self.ny) & (k0 >= 0) & (k0 + 1 < self.nz)

        g = np.zeros(N, dtype=np.float64)
        raw_grad = np.zeros((N, 3), dtype=np.float64)
        grad_norm = np.zeros(N, dtype=np.float64)
        unit_normal = np.zeros((N, 3), dtype=np.float64)
        unit_normal[:, 2] = 1.0

        if not np.any(valid):
            return g, raw_grad, grad_norm, unit_normal, valid

        idx = np.where(valid)[0]
        i0v, j0v, k0v = i0[idx], j0[idx], k0[idx]
        uv, vv, wv = u[idx], v[idx], w[idx]
        phi = self.phi0
        inv_h = 1.0 / self.voxel

        # Gather 8 surrounding values
        p000 = phi[i0v, j0v, k0v]
        p100 = phi[i0v + 1, j0v, k0v]
        p010 = phi[i0v, j0v + 1, k0v]
        p110 = phi[i0v + 1, j0v + 1, k0v]
        p001 = phi[i0v, j0v, k0v + 1]
        p101 = phi[i0v + 1, j0v, k0v + 1]
        p011 = phi[i0v, j0v + 1, k0v + 1]
        p111 = phi[i0v + 1, j0v + 1, k0v + 1]

        # Trilinear interpolation
        one_u, one_v, one_w = 1.0 - uv, 1.0 - vv, 1.0 - wv
        g[idx] = (p000 * one_u * one_v * one_w + p100 * uv * one_v * one_w +
                  p010 * one_u * vv * one_w + p110 * uv * vv * one_w +
                  p001 * one_u * one_v * wv + p101 * uv * one_v * wv +
                  p011 * one_u * vv * wv + p111 * uv * vv * wv)

        # Gradient
        dg_dx = ((p100 - p000) * one_v * one_w + (p110 - p010) * vv * one_w +
                 (p101 - p001) * one_v * wv + (p111 - p011) * vv * wv) * inv_h[0]
        dg_dy = ((p010 - p000) * one_u * one_w + (p110 - p100) * uv * one_w +
                 (p011 - p001) * one_u * wv + (p111 - p101) * uv * wv) * inv_h[1]
        dg_dz = ((p001 - p000) * one_u * one_v + (p101 - p100) * uv * one_v +
                 (p011 - p010) * one_u * vv + (p111 - p110) * uv * vv) * inv_h[2]

        raw_grad[idx, 0] = dg_dx
        raw_grad[idx, 1] = dg_dy
        raw_grad[idx, 2] = dg_dz

        gn = np.sqrt(dg_dx*dg_dx + dg_dy*dg_dy + dg_dz*dg_dz)
        grad_norm[idx] = gn
        mask_gn = gn > 1e-30
        if np.any(mask_gn):
            unit_normal[idx[mask_gn], 0] = dg_dx[mask_gn] / gn[mask_gn]
            unit_normal[idx[mask_gn], 1] = dg_dy[mask_gn] / gn[mask_gn]
            unit_normal[idx[mask_gn], 2] = dg_dz[mask_gn] / gn[mask_gn]

        return g, raw_grad, grad_norm, unit_normal, valid
