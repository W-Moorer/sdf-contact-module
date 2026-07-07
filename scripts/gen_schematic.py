#!/usr/bin/env python3
"""Generate a 2D cross-section schematic of the ring-on-cube torsion scenario.

Shows the xz-plane through the center of the setup.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    'font.family': 'Times New Roman',
    'mathtext.fontset': 'stix',
})

def generate_schematic(path):
    fig, ax = plt.subplots(figsize=(10, 7))

    # --- Parameters (not to scale, but visually clear) ---
    R_i, R_o = 0.4, 0.6
    hw = 2.0           # cube half-width
    penetration = 0.002
    cyl_h = 0.5

    # Cube: top at z=0, bottom at z=-4
    z_cube_bot, z_cube_top = -4.0, 0.0

    # Cylinder: bottom at z=-penetration, top at z=cyl_h-penetration
    z_cyl_bot = -penetration
    z_cyl_top = cyl_h - penetration

    # --- Cube (cross-section rectangle) ---
    cube = plt.Polygon([[-hw, z_cube_bot], [hw, z_cube_bot],
                         [hw, z_cube_top], [-hw, z_cube_top]],
                        closed=True, facecolor='#d5e8d4', edgecolor='#5a7247',
                        linewidth=2, zorder=2)
    ax.add_patch(cube)

    # --- Cylinder walls (cross-section: two rectangles) ---
    # Left wall: x from -R_o to -R_i
    left_wall = plt.Polygon([[-R_o, z_cyl_bot], [-R_i, z_cyl_bot],
                              [-R_i, z_cyl_top], [-R_o, z_cyl_top]],
                             closed=True, facecolor='#b4c7e7', edgecolor='#365897',
                             linewidth=1.5, zorder=3)
    ax.add_patch(left_wall)
    # Right wall: x from R_i to R_o
    right_wall = plt.Polygon([[R_i, z_cyl_bot], [R_o, z_cyl_bot],
                               [R_o, z_cyl_top], [R_i, z_cyl_top]],
                              closed=True, facecolor='#b4c7e7', edgecolor='#365897',
                              linewidth=1.5, zorder=3)
    ax.add_patch(right_wall)
    # Bottom face (annulus cross-section)
    bottom_left = plt.Polygon([[-R_o, z_cyl_bot], [-R_i, z_cyl_bot],
                                [-R_i, z_cyl_bot + 0.005], [-R_o, z_cyl_bot + 0.005]],
                               closed=True, facecolor='#b4c7e7', edgecolor='#365897',
                               linewidth=1.5, zorder=3)
    bottom_right = plt.Polygon([[R_i, z_cyl_bot], [R_o, z_cyl_bot],
                                 [R_o, z_cyl_bot + 0.005], [R_i, z_cyl_bot + 0.005]],
                                closed=True, facecolor='#b4c7e7', edgecolor='#365897',
                                linewidth=1.5, zorder=3)
    ax.add_patch(bottom_left)
    ax.add_patch(bottom_right)
    # Top face
    top_left = plt.Polygon([[-R_o, z_cyl_top], [-R_i, z_cyl_top],
                             [-R_i, z_cyl_top + 0.005], [-R_o, z_cyl_top + 0.005]],
                            closed=True, facecolor='#b4c7e7', edgecolor='#365897',
                            linewidth=1.5, zorder=3)
    top_right = plt.Polygon([[R_i, z_cyl_top], [R_o, z_cyl_top],
                              [R_o, z_cyl_top + 0.005], [R_i, z_cyl_top + 0.005]],
                             closed=True, facecolor='#b4c7e7', edgecolor='#365897',
                             linewidth=1.5, zorder=3)
    ax.add_patch(top_left)
    ax.add_patch(top_right)

    # --- Contact patch (highlighted) ---
    contact_left = plt.Polygon([[-R_o, z_cyl_bot - 0.005], [-R_i, z_cyl_bot - 0.005],
                                 [-R_i, z_cyl_bot], [-R_o, z_cyl_bot]],
                                closed=True, facecolor='#f8cecc', edgecolor='#b85450',
                                linewidth=2.5, zorder=4)
    contact_right = plt.Polygon([[R_i, z_cyl_bot - 0.005], [R_o, z_cyl_bot - 0.005],
                                  [R_o, z_cyl_bot], [R_i, z_cyl_bot]],
                                 closed=True, facecolor='#f8cecc', edgecolor='#b85450',
                                 linewidth=2.5, zorder=4)
    ax.add_patch(contact_left)
    ax.add_patch(contact_right)

    # --- z-axis ---
    ax.annotate('', xy=(0, 1.3), xytext=(0, -4.5),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5), zorder=10)
    ax.text(0.08, 1.35, r'$z$', fontsize=16, color='black')

    # --- x-axis ---
    ax.annotate('', xy=(2.8, 0), xytext=(-2.8, 0),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5), zorder=10)
    ax.text(2.85, 0.1, r'$x$', fontsize=16, color='black')

    # --- Dimension: R_o ---
    y_dim = z_cyl_top + 0.35
    ax.annotate('', xy=(R_o, y_dim), xytext=(0, y_dim),
                arrowprops=dict(arrowstyle='<->', color='#333333', lw=1.2), zorder=10)
    ax.plot([R_o, R_o], [y_dim - 0.05, y_dim + 0.05], color='#333333', lw=1, zorder=10)
    ax.plot([0, 0], [y_dim - 0.05, y_dim + 0.05], color='#333333', lw=1, zorder=10)
    ax.text(R_o / 2, y_dim + 0.08, r'$R_o$', fontsize=15, color='#333333', ha='center')

    # --- Dimension: R_i ---
    y_dim2 = z_cyl_top + 0.65
    ax.annotate('', xy=(R_i, y_dim2), xytext=(0, y_dim2),
                arrowprops=dict(arrowstyle='<->', color='#333333', lw=1.2), zorder=10)
    ax.plot([R_i, R_i], [y_dim2 - 0.05, y_dim2 + 0.05], color='#333333', lw=1, zorder=10)
    ax.plot([0, 0], [y_dim2 - 0.05, y_dim2 + 0.05], color='#333333', lw=1, zorder=10)
    ax.text(R_i / 2, y_dim2 + 0.08, r'$R_i$', fontsize=15, color='#333333', ha='center')

    # --- Dimension: cylinder height h ---
    x_h = R_o + 0.35
    ax.annotate('', xy=(x_h, z_cyl_top), xytext=(x_h, z_cyl_bot),
                arrowprops=dict(arrowstyle='<->', color='#365897', lw=1.2), zorder=10)
    ax.plot([x_h - 0.05, x_h + 0.05], [z_cyl_top, z_cyl_top], color='#365897', lw=1, zorder=10)
    ax.plot([x_h - 0.05, x_h + 0.05], [z_cyl_bot, z_cyl_bot], color='#365897', lw=1, zorder=10)
    ax.text(x_h + 0.12, (z_cyl_top + z_cyl_bot) / 2, r'$h$', fontsize=15, color='#365897',
            ha='left', va='center')

    # --- Dimension: penetration delta (zoomed) ---
    x_d = R_o + 0.8
    ax.annotate('', xy=(x_d, z_cyl_top), xytext=(x_d, 0),
                arrowprops=dict(arrowstyle='<->', color='#a50000', lw=1.2), zorder=10)
    ax.plot([x_d - 0.05, x_d + 0.05], [z_cyl_top, z_cyl_top], color='#a50000', lw=1, zorder=10)
    ax.plot([x_d - 0.05, x_d + 0.05], [0, 0], color='#a50000', lw=1, zorder=10)
    ax.text(x_d + 0.12, (z_cyl_top + 0) / 2, r'$\delta$', fontsize=15, color='#a50000',
            ha='left', va='center')

    # --- Rotation arrow ---
    ax.annotate('', xy=(-0.15, z_cyl_top + 0.12), xytext=(0.15, z_cyl_top + 0.12),
                arrowprops=dict(arrowstyle='->', color='#a50000', lw=2,
                                connectionstyle='arc3,rad=0.3'), zorder=10)
    ax.text(0, z_cyl_top + 0.28, r'$\Omega \, e_z$', fontsize=16, color='#a50000',
            ha='center', va='bottom')

    # --- Body labels ---
    ax.text(-hw * 0.7, (z_cube_top + z_cube_bot) / 2,
            'Cube (Body B)', fontsize=14, color='#5a7247', style='italic',
            ha='center', va='center', zorder=10,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#5a7247', alpha=0.9))
    ax.text(-R_o - 0.4, (z_cyl_top + z_cyl_bot) / 2,
            'Hollow\nCylinder\n(Body A)', fontsize=13, color='#365897', style='italic',
            ha='right', va='center', zorder=10,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#365897', alpha=0.9))

    # --- Contact patch label ---
    ax.annotate('Contact\nPatch', xy=(R_o * 0.5, z_cyl_bot), xytext=(1.8, -0.6),
                fontsize=12, color='#b85450', ha='center', va='center', zorder=10,
                arrowprops=dict(arrowstyle='->', color='#b85450', lw=1.5),
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#b85450', alpha=0.9))

    # --- Normal force arrows ---
    for xf in [R_i + 0.05, (R_i + R_o) / 2, R_o - 0.05]:
        ax.annotate('', xy=(xf, z_cyl_bot - 0.15), xytext=(xf, z_cyl_bot),
                    arrowprops=dict(arrowstyle='->', color='#d6604d', lw=1.2), zorder=10)
    ax.text(R_o * 0.5, z_cyl_bot - 0.25, r'$f_n$', fontsize=13, color='#d6604d', ha='center')

    # --- Grid and axes ---
    ax.set_xlim(-2.8, 3.2)
    ax.set_ylim(-4.8, 1.5)
    ax.set_aspect('equal')
    ax.set_axis_off()
    ax.set_title('Ring-on-Cube Torsion: Cross-Section (xz-plane)', fontsize=15, pad=20)

    plt.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {path}")

if __name__ == '__main__':
    out_dir = os.path.join(os.path.dirname(__file__), '..', 'assets', 'figures')
    os.makedirs(out_dir, exist_ok=True)
    generate_schematic(os.path.join(out_dir, 'scenario_schematic.png'))
