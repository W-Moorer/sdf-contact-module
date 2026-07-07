#!/usr/bin/env python3
"""Generate a 3D schematic of the ring-on-cube torsion scenario."""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os

plt.rcParams.update({
    'font.family': 'Times New Roman',
    'mathtext.fontset': 'stix',
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
})

def _cylinder_mesh(R_i, R_o, height, n_az=32):
    verts = []
    faces = []
    def add_band(r, z_start, z_end):
        base = len(verts)
        for k in range(2):
            z = z_start if k == 0 else z_end
            for j in range(n_az):
                theta = j * 2 * np.pi / n_az
                verts.append([r * np.cos(theta), r * np.sin(theta), z])
        for j in range(n_az):
            j1 = (j + 1) % n_az
            i0 = base + j
            i1 = base + j1
            i2 = base + n_az + j
            i3 = base + n_az + j1
            faces.append([i0, i2, i1])
            faces.append([i1, i2, i3])
    add_band(R_i, 0, height)
    add_band(R_o, 0, height)
    def add_annulus(r_inner, r_outer, z, n_rad=3):
        base = len(verts)
        for j in range(n_az):
            for i in range(n_rad + 1):
                r = r_inner + i * (r_outer - r_inner) / n_rad
                theta = j * 2 * np.pi / n_az
                verts.append([r * np.cos(theta), r * np.sin(theta), z])
        for j in range(n_az):
            j1 = (j + 1) % n_az
            for i in range(n_rad):
                i0 = base + j * (n_rad + 1) + i
                i1 = i0 + 1
                i2 = base + j1 * (n_rad + 1) + i
                i3 = i2 + 1
                faces.append([i0, i1, i2])
                faces.append([i1, i3, i2])
    add_annulus(R_i, R_o, 0)
    add_annulus(R_i, R_o, height)
    return np.array(verts), np.array(faces)

def _cube_mesh(hw, z_bottom, z_top):
    verts = []
    faces = []
    for axis, sign, ua, va in [
        (0, -1, 1, 2), (0, 1, 1, 2),
        (1, -1, 0, 2), (1, 1, 0, 2),
        (2, -1, 0, 1), (2, 1, 0, 1)
    ]:
        base = len(verts)
        for i in range(2):
            for j in range(2):
                p = [0.0] * 3
                p[axis] = (z_bottom if axis == 2 and sign == -1 else
                           z_top if axis == 2 and sign == 1 else
                           sign * hw)
                u = -hw + 2 * hw * i
                v = -hw + 2 * hw * j
                p[ua] = u
                p[va] = v
                verts.append(p)
        i0, i1, i2, i3 = base, base+1, base+2, base+3
        if sign == 1:
            faces.append([i0, i2, i1])
            faces.append([i1, i2, i3])
        else:
            faces.append([i0, i1, i2])
            faces.append([i1, i3, i2])
    return np.array(verts), np.array(faces)

def generate_schematic(path):
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=22, azim=-55)

    # --- Cube ---
    hw = 2.0
    cv, cf = _cube_mesh(hw, -4.0, 0.0)
    cube_polys = Poly3DCollection([cv[f] for f in cf],
                                  facecolor='#d5e8d4', edgecolor='#82b366',
                                  linewidth=0.3, alpha=0.7)
    ax.add_collection3d(cube_polys)

    # --- Hollow cylinder ---
    R_i, R_o = 0.4, 0.6
    height = 0.5
    penetration = 0.002
    cyl_z = -penetration
    cyv, cyf = _cylinder_mesh(R_i, R_o, height, n_az=32)
    cyv[:, 2] += cyl_z
    cyl_polys = Poly3DCollection([cyv[f] for f in cyf],
                                 facecolor='#dae8fc', edgecolor='#6c8ebf',
                                 linewidth=0.3, alpha=0.85)
    ax.add_collection3d(cyl_polys)

    # --- Contact patch (annular ring on cube top) ---
    theta = np.linspace(0, 2*np.pi, 64)
    r_ring = np.linspace(R_i, R_o, 10)
    T, R = np.meshgrid(theta, r_ring)
    patch_x = R * np.cos(T)
    patch_y = R * np.sin(T)
    patch_z = np.zeros_like(patch_x)
    ax.plot_surface(patch_x, patch_y, patch_z,
                    color='#f8cecc', alpha=0.6, edgecolor='#b85450',
                    linewidth=0.15, antialiased=True)

    # --- Dimension arrows and labels ---
    # R_i arrow
    ax.plot([0, R_i], [0, 0], [height + cyl_z + 0.15]*2,
            color='#333333', lw=1.5, marker='|', markersize=8)
    ax.text(R_i / 2, -0.12, height + cyl_z + 0.18,
            r'$R_i$', fontsize=14, color='#333333', ha='center')

    # R_o arrow
    ax.plot([0, R_o], [0, 0], [height + cyl_z + 0.3]*2,
            color='#333333', lw=1.5, marker='|', markersize=8)
    ax.text(R_o / 2, -0.12, height + cyl_z + 0.33,
            r'$R_o$', fontsize=14, color='#333333', ha='center')

    # Penetration delta arrow
    dz = penetration + 0.03
    ax.plot([0.85, 0.85], [0, 0], [-dz, 0], color='#333333', lw=1.5)
    ax.plot([0.82, 0.88], [0, 0], [-dz, -dz], color='#333333', lw=1.0)
    ax.plot([0.82, 0.88], [0, 0], [0, 0], color='#333333', lw=1.0)
    ax.plot([0.85, 0.83], [0, 0], [0, 0.012], color='#333333', lw=1.0)
    ax.plot([0.85, 0.87], [0, 0], [0, 0.012], color='#333333', lw=1.0)
    ax.text(0.95, 0, -dz / 2, r'$\delta$', fontsize=13,
            color='#333333', va='center')

    # Omega arrow (rotation)
    ax.text(0.0, 0.7, height + cyl_z + 0.1,
            r'$\Omega \, e_z$', fontsize=15, color='#a50000',
            ha='center', va='bottom')
    # Rotation arc
    for r_arc in [R_i * 1.1, R_o * 0.9]:
        arc_theta = np.linspace(np.pi/4, 3*np.pi/4, 20)
        ax.plot(r_arc * np.cos(arc_theta), r_arc * np.sin(arc_theta),
                [height + cyl_z] * 20,
                color='#a50000', lw=1.5, linestyle='--')
    ax.quiver(R_o * 0.9 * np.cos(3*np.pi/4), R_o * 0.9 * np.sin(3*np.pi/4),
              height + cyl_z,
              -R_o * 0.9 * np.sin(3*np.pi/4) * 0.08,
              R_o * 0.9 * np.cos(3*np.pi/4) * 0.08, 0,
              color='#a50000', lw=2)

    # Labels
    ax.text(-0.8, -2.0, -3.5, 'Cube (Body B)', fontsize=13,
            color='#333333', style='italic')
    ax.text(-0.4, 0.9, height/2 + cyl_z, 'Hollow\nCylinder\n(Body A)',
            fontsize=13, color='#333333', style='italic', ha='center')
    ax.text(0.0, 0.0, 0.1, 'Contact\nPatch', fontsize=11,
            color='#b85450', ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='#b85450',
                      alpha=0.8))

    # Axes
    ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.5, 2.5); ax.set_zlim(-4.5, 1.2)
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_box_aspect((1, 1, 1.2))
    # Remove grid and ticks for clean schematic look
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('none')
    ax.yaxis.pane.set_edgecolor('none')
    ax.zaxis.pane.set_edgecolor('none')
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])

    plt.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Schematic saved to: {path}")

if __name__ == '__main__':
    out = os.path.join(os.path.dirname(__file__), '..', 'assets', 'figures',
                       'scenario_schematic.png')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    generate_schematic(out)
