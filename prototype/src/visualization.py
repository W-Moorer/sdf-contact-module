import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def save_all_figures(X_q, w_q, result, params, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    R_o = params['R_o']
    R_i = params['R_i']
    mu = params['mu']
    N_numeric = result['N_numeric']
    T_z_numeric = result['T_z_numeric']
    p_hat = result['p_hat']
    forces_per_point = result['forces_per_point']
    torques_per_point = result['torques_per_point']
    gaps = result['gaps']

    fig1 = _plot_3d_pressure(X_q, p_hat, R_i, R_o)
    fig1.savefig(os.path.join(output_dir, 'contact_pressure_3d.png'), dpi=150, bbox_inches='tight')
    plt.close(fig1)

    fig2 = _plot_top_pressure(X_q, p_hat, R_i, R_o)
    fig2.savefig(os.path.join(output_dir, 'contact_pressure_top.png'), dpi=150, bbox_inches='tight')
    plt.close(fig2)

    fig3 = _plot_friction_vectors(X_q, forces_per_point, p_hat, R_i, R_o)
    fig3.savefig(os.path.join(output_dir, 'friction_vectors_top.png'), dpi=150, bbox_inches='tight')
    plt.close(fig3)

    fig4 = _plot_torque_density(X_q, torques_per_point, R_i, R_o)
    fig4.savefig(os.path.join(output_dir, 'torque_density_top.png'), dpi=150, bbox_inches='tight')
    plt.close(fig4)

def _plot_3d_pressure(X_q, p_hat, R_i, R_o):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(X_q[:, 0], X_q[:, 1], X_q[:, 2],
                    c=p_hat, cmap='viridis', s=15, alpha=0.8)
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_title('3D Contact Patch Pressure $\\hat{p}_q$')
    ax.set_box_aspect((1, 1, 0.3))
    plt.colorbar(sc, ax=ax, label='$\\hat{p}_q$')
    return fig

def _plot_top_pressure(X_q, p_hat, R_i, R_o):
    fig, ax = plt.subplots(figsize=(8, 8))
    sc = ax.scatter(X_q[:, 0], X_q[:, 1], c=p_hat,
                    cmap='viridis', s=20, alpha=0.8)
    theta = np.linspace(0, 2*np.pi, 200)
    ax.plot(R_o*np.cos(theta), R_o*np.sin(theta), 'k--', lw=1, alpha=0.5)
    ax.plot(R_i*np.cos(theta), R_i*np.sin(theta), 'k--', lw=1, alpha=0.5)
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title('Top-View Contact Pressure')
    ax.set_aspect('equal')
    plt.colorbar(sc, ax=ax, label='$\\hat{p}_q$')
    return fig

def _plot_friction_vectors(X_q, forces_per_point, p_hat, R_i, R_o):
    f_t = forces_per_point.copy()
    mask = p_hat > 0.01 * np.max(p_hat)
    X_sub = X_q[mask]
    f_sub = f_t[mask]
    fig, ax = plt.subplots(figsize=(8, 8))
    q = ax.quiver(X_sub[:, 0], X_sub[:, 1],
                  f_sub[:, 0], f_sub[:, 1],
                  np.linalg.norm(f_sub[:, :2], axis=1),
                  cmap='plasma', width=0.003)
    theta = np.linspace(0, 2*np.pi, 200)
    ax.plot(R_o*np.cos(theta), R_o*np.sin(theta), 'k--', lw=1, alpha=0.5)
    ax.plot(R_i*np.cos(theta), R_i*np.sin(theta), 'k--', lw=1, alpha=0.5)
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title('Friction Traction Vectors (top view)')
    ax.set_aspect('equal')
    plt.colorbar(q, ax=ax, label='friction magnitude')
    return fig

def _plot_torque_density(X_q, torques_per_point, R_i, R_o):
    dT_z = torques_per_point[:, 2]
    fig, ax = plt.subplots(figsize=(8, 8))
    sc = ax.scatter(X_q[:, 0], X_q[:, 1], c=dT_z,
                    cmap='RdBu_r', s=20, alpha=0.8,
                    vmin=-np.max(np.abs(dT_z)), vmax=np.max(np.abs(dT_z)))
    theta = np.linspace(0, 2*np.pi, 200)
    ax.plot(R_o*np.cos(theta), R_o*np.sin(theta), 'k--', lw=1, alpha=0.5)
    ax.plot(R_i*np.cos(theta), R_i*np.sin(theta), 'k--', lw=1, alpha=0.5)
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title('Torque Density $dT_z$ per Point')
    ax.set_aspect('equal')
    plt.colorbar(sc, ax=ax, label='$dT_z$')
    return fig

def save_convergence_plot(results, output_dir, suffix=''):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    labels = [r['label'] for r in results]
    n_pts = np.array([r['num_quadrature_points'] for r in results])
    errs = np.array([r['relative_error_Tz'] for r in results])
    rsp_ana = results[0]['Rspin_analytic']
    rsp_errs = np.abs(np.array([r['Rspin_numeric'] for r in results]) - rsp_ana) / rsp_ana

    ax1.loglog(n_pts, errs, 'o-', color='#1f77b4', lw=2, markersize=8)
    for i, l in enumerate(labels):
        ax1.annotate(l, (n_pts[i], errs[i]), textcoords="offset points",
                     xytext=(0, 10), ha='center', fontsize=8)
    ref_rate = errs[0] * (n_pts[0] / n_pts) ** 1.0
    ax1.loglog(n_pts, ref_rate, '--', color='gray', alpha=0.5, label='O(1/N) ref')
    ax1.set_xlabel('number of quadrature points')
    ax1.set_ylabel('relative torque error')
    ax1.set_title('Torque Error vs Mesh Resolution')
    ax1.legend(); ax1.grid(True, which='both', alpha=0.3)

    ax2.semilogy(n_pts, rsp_errs, 's-', color='#d62728', lw=2, markersize=8)
    for i, l in enumerate(labels):
        ax2.annotate(l, (n_pts[i], rsp_errs[i]), textcoords="offset points",
                     xytext=(0, 10), ha='center', fontsize=8)
    ax2.set_xlabel('number of quadrature points')
    ax2.set_ylabel('relative R_spin error')
    ax2.set_title('Equivalent Radius Error vs Mesh Resolution')
    ax2.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'figures', f'torque_error_vs_mesh_resolution{suffix}.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)

def save_sdf_convergence_plot(sdf_results, output_dir):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    res_list = sorted([r for r in sdf_results],
                      key=lambda x: int(x['label'].replace('sdf_res','')))
    resolutions = np.array([int(r['label'].replace('sdf_res','')) for r in res_list])
    errs = np.array([r['relative_error_Tz'] for r in res_list])
    rsp_errs = np.abs(np.array([r['Rspin_numeric'] for r in res_list])
                      - res_list[0]['Rspin_analytic']) / res_list[0]['Rspin_analytic']

    ax1.loglog(resolutions, errs, 'o-', color='#2ca02c', lw=2, markersize=8)
    ref = errs[0] * (resolutions[0] / resolutions)
    ax1.loglog(resolutions, ref, '--', color='gray', alpha=0.5, label='O(1/N) ref')
    ax1.set_xlabel('SDF resolution')
    ax1.set_ylabel('relative torque error')
    ax1.set_title('Torque Error vs SDF Resolution')
    ax1.legend(); ax1.grid(True, which='both', alpha=0.3)

    ax2.semilogy(resolutions, rsp_errs, 's-', color='#9467bd', lw=2, markersize=8)
    ax2.set_xlabel('SDF resolution')
    ax2.set_ylabel('relative R_spin error')
    ax2.set_title('Equivalent Radius Error vs SDF Resolution')
    ax2.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'figures', 'torque_error_vs_sdf_resolution.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)
