#!/usr/bin/env python3
"""Generate a 2D cross-section + 3D overview of the ring-on-cube torsion scenario."""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    'font.family': 'Times New Roman',
    'mathtext.fontset': 'stix',
})

# --- Global parameters ---
R_i, R_o = 0.4, 0.6
hw = 2.0
penetration = 0.002
cyl_h = 0.5
z_cube_bot, z_cube_top = -4.0, 0.0
z_cyl_bot = -penetration
z_cyl_top = cyl_h - penetration

def _draw_dimension(ax, p1, p2, label, offset, color='#333333', fontsize=14,
                    ha='center', va='center', lw=1.2):
    """Draw a dimension line with arrows and label."""
    ax.annotate('', xy=p2, xytext=p1,
                arrowprops=dict(arrowstyle='<->', color=color, lw=lw), zorder=10)
    # Tick marks
    for p in [p1, p2]:
        ax.plot([p[0] - 0.04, p[0] + 0.04], [p[1], p[1]],
                color=color, lw=1, zorder=10)
    mx = (p1[0] + p2[0]) / 2 + offset[0]
    my = (p1[1] + p2[1]) / 2 + offset[1]
    ax.text(mx, my, label, fontsize=fontsize, color=color, ha=ha, va=va)

def generate_2d(path):
    fig, ax = plt.subplots(figsize=(12, 8))

    # --- Cube ---
    cube = plt.Polygon([[-hw, z_cube_bot], [hw, z_cube_bot],
                         [hw, z_cube_top], [-hw, z_cube_top]],
                        closed=True, facecolor='#d5e8d4', edgecolor='#5a7247',
                        linewidth=2.5, zorder=2)
    ax.add_patch(cube)

    # --- Cylinder walls (left side) ---
    left_wall = plt.Polygon([[-R_o, z_cyl_bot], [-R_i, z_cyl_bot],
                              [-R_i, z_cyl_top], [-R_o, z_cyl_top]],
                             closed=True, facecolor='#b4c7e7', edgecolor='#365897',
                             linewidth=1.5, zorder=3)
    ax.add_patch(left_wall)
    # Right wall
    right_wall = plt.Polygon([[R_i, z_cyl_bot], [R_o, z_cyl_bot],
                               [R_o, z_cyl_top], [R_i, z_cyl_top]],
                              closed=True, facecolor='#b4c7e7', edgecolor='#365897',
                              linewidth=1.5, zorder=3)
    ax.add_patch(right_wall)

    # --- Contact patch (red) ---
    contact_left = plt.Polygon([[-R_o, z_cyl_bot - 0.01], [-R_i, z_cyl_bot - 0.01],
                                 [-R_i, z_cyl_bot], [-R_o, z_cyl_bot]],
                                closed=True, facecolor='#f8cecc', edgecolor='#b85450',
                                linewidth=2.5, zorder=4)
    contact_right = plt.Polygon([[R_i, z_cyl_bot - 0.01], [R_o, z_cyl_bot - 0.01],
                                  [R_o, z_cyl_bot], [R_i, z_cyl_bot]],
                                 closed=True, facecolor='#f8cecc', edgecolor='#b85450',
                                 linewidth=2.5, zorder=4)
    ax.add_patch(contact_left)
    ax.add_patch(contact_right)

    # --- Centerline (axis of symmetry, limited to cylinder height) ---
    ax.plot([0, 0], [z_cyl_bot - 0.1, z_cyl_top + 0.1],
            color='#888888', lw=0.8, ls='--', zorder=1)

    # --- z-axis (left side) ---
    ax.annotate('', xy=(-hw - 0.2, 1.5), xytext=(-hw - 0.2, -4.5),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5), zorder=10)
    ax.text(-hw - 0.25, 1.55, r'$z$', fontsize=16, color='black', ha='right')

    # --- x-axis (bottom) ---
    ax.annotate('', xy=(hw + 0.5, -0.08), xytext=(-(hw + 0.5), -0.08),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5), zorder=10)
    ax.text(hw + 0.55, -0.05, r'$x$', fontsize=16, color='black', ha='left')

    # --- Dimension: R_i (right side, lowest) ---
    y_ri = z_cyl_top + 0.3
    _draw_dimension(ax, (0, y_ri), (R_i, y_ri), r'$R_i$', offset=(0, 0.12))

    # --- Dimension: R_o (right side, above R_i) ---
    y_ro = z_cyl_top + 0.85
    _draw_dimension(ax, (0, y_ro), (R_o, y_ro), r'$R_o$', offset=(0, 0.12))

    # --- Dimension: cylinder height h (right side) ---
    x_h = R_o + 0.55
    _draw_dimension(ax, (x_h, z_cyl_bot), (x_h, z_cyl_top), r'$h$',
                    offset=(0.15, 0), color='#365897', fontsize=14)

    # --- Detail indicator on main axes (shows inset zoom region) ---
    detail_rect = plt.Rectangle(
        (0.32, -0.12), 0.36, 0.24,
        linewidth=1.5, edgecolor='#b85450', facecolor='none',
        linestyle='--', zorder=10)
    ax.add_patch(detail_rect)
    ax.text(0.50, 0.18, 'Detail', fontsize=9, color='#b85450',
            ha='center', va='bottom', zorder=10, style='italic')

    # --- δ extension lines on main diagram ---
    # Penetration is 0.002 (invisible at this scale); use extension lines
    # from cube surface and cylinder bottom to a reference column
    x_d_ref = R_o + 0.85
    ax.plot([R_o + 0.03, x_d_ref + 0.06], [0, 0],
            color='#a50000', lw=0.5, ls=':', zorder=5)
    ax.plot([R_o + 0.03, x_d_ref + 0.06], [z_cyl_bot, z_cyl_bot],
            color='#a50000', lw=0.5, ls=':', zorder=5)
    # Small δ ticks at the reference column
    ax.plot([x_d_ref - 0.03, x_d_ref + 0.03], [0, 0],
            color='#a50000', lw=1.0, zorder=10)
    ax.plot([x_d_ref - 0.03, x_d_ref + 0.03], [z_cyl_bot, z_cyl_bot],
            color='#a50000', lw=1.0, zorder=10)
    ax.text(x_d_ref + 0.1, z_cyl_bot / 2, r'$\delta$',
            fontsize=12, color='#a50000', ha='left', va='center')

    # --- Inset: zoom on contact region showing penetration δ ---
    ins_ax = fig.add_axes([0.70, 0.70, 0.26, 0.26])

    # Visual-coordinate layout (penetration exaggerated for clarity)
    vis_pen = 0.08          # exaggerated penetration in visual units
    cube_top_y = 0.0        # cube surface at y = 0
    cube_bot_y = -0.40      # bottom of visible cube region
    cyl_top_y = 0.38        # top of visible cylinder region
    cyl_bot_y = cube_top_y - vis_pen  # cylinder bottom (below surface)
    x_L, x_R = 0.12, 0.62   # cylinder wall horizontal extent
    cube_xL = 0.0            # cube extends beyond cylinder on both sides
    cube_xR = 0.80

    # Cube body (green)
    ins_ax.add_patch(plt.Rectangle(
        (cube_xL, cube_bot_y), cube_xR - cube_xL, -cube_bot_y,
        facecolor='#d5e8d4', edgecolor='#5a7247', lw=1.5, zorder=2))

    # Cylinder wall (blue) — spans from cyl_bot_y to cyl_top_y
    ins_ax.add_patch(plt.Rectangle(
        (x_L, cyl_bot_y), x_R - x_L, cyl_top_y - cyl_bot_y,
        facecolor='#b4c7e7', edgecolor='#365897', lw=1.5, zorder=3))

    # Contact patch (red thin strip at the interface)
    ins_ax.add_patch(plt.Rectangle(
        (x_L, cube_top_y - vis_pen * 0.25), x_R - x_L, vis_pen * 0.5,
        facecolor='#f8cecc', edgecolor='#b85450', lw=2.0, zorder=4))

    # z = 0 reference dashed line (cube surface)
    ins_ax.axhline(y=cube_top_y, color='#5a7247', lw=0.8, ls='--', zorder=1)
    ins_ax.text(cube_xR + 0.01, cube_top_y, r'$z=0$',
                fontsize=8, color='#5a7247', va='center', ha='left')

    # δ dimension: extension lines + double-headed arrow
    dim_x = x_R + 0.06
    ext_left = x_R + 0.01
    ext_right = dim_x + 0.03
    # Extension line from z = 0
    ins_ax.plot([ext_left, ext_right], [cube_top_y, cube_top_y],
                color='#a50000', lw=0.5, ls=':', zorder=5)
    # Extension line from cylinder bottom (penetrated depth)
    ins_ax.plot([ext_left, ext_right], [cyl_bot_y, cyl_bot_y],
                color='#a50000', lw=0.5, ls=':', zorder=5)
    # Dimension arrow
    ins_ax.annotate('', xy=(dim_x, cyl_bot_y), xytext=(dim_x, cube_top_y),
                    arrowprops=dict(arrowstyle='<->', color='#a50000', lw=1.2),
                    zorder=10)
    # Tick marks at arrow endpoints
    for yy in [cube_top_y, cyl_bot_y]:
        ins_ax.plot([dim_x - 0.015, dim_x + 0.015], [yy, yy],
                    color='#a50000', lw=1.0, zorder=10)
    # δ label
    ins_ax.text(dim_x + 0.025, (cube_top_y + cyl_bot_y) / 2, r'$\delta$',
                fontsize=14, color='#a50000', ha='left', va='center',
                fontweight='bold')

    ins_ax.set_xlim(-0.05, 0.88)
    ins_ax.set_ylim(-0.48, 0.46)
    ins_ax.tick_params(labelsize=7, pad=1)
    ins_ax.set_title(r'Detail ($\delta$ exaggerated)', fontsize=9, pad=3, color='#666666')
    for spine in ins_ax.spines.values():
        spine.set_edgecolor('#b85450')
        spine.set_linewidth(1.5)

    # --- Rotation arrow (above cylinder) ---
    r_arc = R_i * 0.7
    arc_theta = np.linspace(np.pi/3, 2*np.pi/3, 30)
    ax.plot(r_arc * np.cos(arc_theta), r_arc * np.sin(arc_theta) + z_cyl_top + 0.05,
            color='#a50000', lw=2, zorder=10)
    ax.annotate('', xy=(r_arc * np.cos(2*np.pi/3), r_arc * np.sin(2*np.pi/3) + z_cyl_top + 0.05),
                xytext=(r_arc * np.cos(2*np.pi/3 - 0.15),
                        r_arc * np.sin(2*np.pi/3 - 0.15) + z_cyl_top + 0.05),
                arrowprops=dict(arrowstyle='->', color='#a50000', lw=2), zorder=10)
    ax.text(0, z_cyl_top + 0.75, r'$\Omega \, e_z$', fontsize=15, color='#a50000',
            ha='center', va='bottom')

    # --- Body labels ---
    ax.text(-hw * 0.7, (z_cube_top + z_cube_bot) / 2,
            'Cube (Body B)', fontsize=14, color='#5a7247', style='italic',
            ha='center', va='center', zorder=10,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#5a7247', alpha=0.9))

    ax.text(-R_o - 0.5, (z_cyl_top + z_cyl_bot) / 2,
            'Hollow\nCylinder\n(Body A)', fontsize=13, color='#365897', style='italic',
            ha='right', va='center', zorder=10,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#365897', alpha=0.9))

    # --- Contact patch label ---
    ax.annotate('Contact\nPatch', xy=(R_o * 0.7, z_cyl_bot - 0.005), xytext=(2.0, -0.7),
                fontsize=12, color='#b85450', ha='center', va='center', zorder=10,
                arrowprops=dict(arrowstyle='->', color='#b85450', lw=1.5),
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#b85450', alpha=0.9))

    # --- Normal force arrows (upward: cube pushes cylinder in +z) ---
    for xf in [R_i + 0.06, (R_i + R_o) / 2, R_o - 0.06]:
        ax.annotate('', xy=(xf, 0.18), xytext=(xf, z_cyl_bot - 0.01),
                    arrowprops=dict(arrowstyle='->', color='#d6604d', lw=1.5), zorder=10)
    ax.text(0.72, 0.10, r'$f_n$', fontsize=13, color='#d6604d',
            ha='left', va='center', zorder=10)

    ax.set_xlim(-2.8, 3.5)
    ax.set_ylim(-4.8, 2.0)
    ax.set_aspect('equal')
    ax.set_axis_off()
    ax.set_title('(a) Cross-Section (xz-plane)', fontsize=15, pad=20)

    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {path}")

def generate_3d(path):
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=22, azim=-55)

    # --- Cube ---
    cv = [
        [[-hw, -hw, z_cube_bot], [ hw, -hw, z_cube_bot], [ hw,  hw, z_cube_bot], [-hw,  hw, z_cube_bot]],
        [[-hw, -hw, z_cube_top], [ hw, -hw, z_cube_top], [ hw,  hw, z_cube_top], [-hw,  hw, z_cube_top]],
        [[-hw, -hw, z_cube_bot], [ hw, -hw, z_cube_bot], [ hw, -hw, z_cube_top], [-hw, -hw, z_cube_top]],
        [[-hw,  hw, z_cube_bot], [ hw,  hw, z_cube_bot], [ hw,  hw, z_cube_top], [-hw,  hw, z_cube_top]],
        [[-hw, -hw, z_cube_bot], [-hw,  hw, z_cube_bot], [-hw,  hw, z_cube_top], [-hw, -hw, z_cube_top]],
        [[ hw, -hw, z_cube_bot], [ hw,  hw, z_cube_bot], [ hw,  hw, z_cube_top], [ hw, -hw, z_cube_top]],
    ]
    ax.add_collection3d(Poly3DCollection(cv, facecolor='#d5e8d4', edgecolor='#5a7247',
                                          linewidth=0.6, alpha=0.55))

    # --- Cylinder (wireframe) ---
    n = 60
    theta = np.linspace(0, 2*np.pi, n)
    for r in [R_i, R_o]:
        ax.plot(r*np.cos(theta), r*np.sin(theta), [z_cyl_top]*n,
                color='#365897', lw=1.0, zorder=5)
        ax.plot(r*np.cos(theta), r*np.sin(theta), [z_cyl_bot]*n,
                color='#365897', lw=1.0, zorder=5)
    for k in range(0, n, n//8):
        for r in [R_i, R_o]:
            ax.plot([r*np.cos(theta[k])]*2, [r*np.sin(theta[k])]*2,
                    [z_cyl_bot, z_cyl_top], color='#365897', lw=0.6, zorder=5)

    # --- Top face annulus ---
    for r in [R_i, R_o]:
        ax.plot(r*np.cos(theta), r*np.sin(theta), [z_cyl_top]*n,
                color='#365897', lw=0.8, zorder=5)
    for k in range(0, n, n//8):
        ax.plot([R_i*np.cos(theta[k]), R_o*np.cos(theta[k])],
                [R_i*np.sin(theta[k]), R_o*np.sin(theta[k])],
                [z_cyl_top]*2, color='#365897', lw=0.5, zorder=5)

    # --- Contact patch (red ring at bottom) ---
    for r in [R_i, R_o]:
        ax.plot(r*np.cos(theta), r*np.sin(theta), [z_cyl_bot]*n,
                color='#b85450', lw=2.5, zorder=6)
    for k in range(0, n, n//8):
        ax.plot([R_i*np.cos(theta[k]), R_o*np.cos(theta[k])],
                [R_i*np.sin(theta[k]), R_o*np.sin(theta[k])],
                [z_cyl_bot]*2, color='#b85450', lw=1.0, zorder=6)

    # --- Rotation arc ---
    r_arc = R_i * 0.6
    arc_theta = np.linspace(np.pi/4, 3*np.pi/4, 25)
    ax.plot(r_arc*np.cos(arc_theta), r_arc*np.sin(arc_theta),
            [z_cyl_top + 0.05]*25, color='#a50000', lw=2.5, zorder=10)
    # Arrowhead (manual 3-point triangle)
    tip_i = -3
    tx, ty, tz = (r_arc*np.cos(arc_theta[tip_i]),
                  r_arc*np.sin(arc_theta[tip_i]),
                  z_cyl_top + 0.05)
    ax.plot([tx, tx - 0.03], [ty, ty + 0.05], [tz]*2,
            color='#a50000', lw=2.5, zorder=10)
    ax.plot([tx, tx + 0.05], [ty, ty - 0.02], [tz]*2,
            color='#a50000', lw=2.5, zorder=10)

    # --- Labels ---
    ax.text(0, 0, z_cyl_top + 0.35, r'$\Omega \, e_z$', fontsize=16, color='#a50000',
            ha='center', va='bottom', zorder=10)

    ax.text(0, 0, z_cyl_bot - 0.25, 'Contact\nPatch', fontsize=12, color='#b85450',
            ha='center', va='top', zorder=10,
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='#b85450', alpha=0.85))

    ax.text(-hw * 0.6, -hw * 0.5, (z_cube_top + z_cube_bot) / 2,
            'Cube (Body B)', fontsize=13, color='#5a7247', style='italic', zorder=10)

    ax.text(R_o + 0.4, 0, (z_cyl_top + z_cyl_bot) / 2,
            'Hollow Cylinder\n(Body A)', fontsize=13, color='#365897', style='italic',
            zorder=10)

    # --- Axis ---
    ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.5, 2.5); ax.set_zlim(-4.5, 1.5)
    ax.set_xlabel('x', fontsize=13); ax.set_ylabel('y', fontsize=13); ax.set_zlabel('z', fontsize=13)
    ax.set_box_aspect((1, 1, 1.2))
    ax.grid(False)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False; pane.set_edgecolor('none')
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])

    ax.set_title('(b) 3D Overview', fontsize=15, pad=15)
    plt.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {path}")

if __name__ == '__main__':
    out_dir = os.path.join(os.path.dirname(__file__), '..', 'assets', 'figures')
    os.makedirs(out_dir, exist_ok=True)
    generate_2d(os.path.join(out_dir, 'scenario_schematic_2d.png'))
    generate_3d(os.path.join(out_dir, 'scenario_schematic_3d.png'))
