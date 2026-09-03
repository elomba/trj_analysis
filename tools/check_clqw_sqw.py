#!/usr/bin/env python3
"""
Consistency Checker & Plotter for Longitudinal Current Spectrum C_L(Q,w) and Dynamic Structure Factor S(Q,w)

Theoretical Relationships:
1. Frequency Domain:
   C_L(Q, w) = (w^2 / Q^2) * S_raw(Q, w)  or  (w^2 / Q) * S_raw(Q, w)
   C_{L,s}(Q, w) = (w^2 / Q^2) * S_{s,raw}(Q, w)
   where w is the angular frequency (2*pi*f in physical/reduced units).

2. Time Domain:
   J_{L,s}(Q, t=0) = <v_L^2> = k_B * T / m
   J_L(Q, t) = - (1 / Q^2) * d^2 F(Q, t) / dt^2

Outputs:
   - Screen plot (plt.show())
   - PDF plots: clqw_sqw_comparison.pdf and clqw_spectrum_check.pdf (no PNG)
   - Text report: check_clqw_sqw_report.txt

Usage:
   python3 check_clqw_sqw.py [data_dir]
"""

import sys
import os

# Auto-enable Conda python if numpy/matplotlib missing in default interpreter
try:
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError:
    miniconda_py = '/usr/local/miniconda3/bin/python3'
    if os.path.exists(miniconda_py) and sys.executable != miniconda_py:
        os.execv(miniconda_py, [miniconda_py] + sys.argv)
    else:
        raise ImportError("Numpy and Matplotlib are required. Please run with Conda python.")

import math
import re


def parse_dat_file(filepath):
    """Parse ASCII dat file into header lines and numpy array."""
    if not os.path.exists(filepath):
        return None, None, None

    headers = []
    col_names = []
    data = []

    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                headers.append(stripped)
                cols = re.findall(r'[^\s#]+(?:\([^\)]+\)[^\s#]*)?|\S+', stripped[1:])
                if cols:
                    col_names = cols
            else:
                try:
                    vals = [float(x) for x in stripped.split()]
                    if vals:
                        data.append(vals)
                except ValueError:
                    continue

    if not data:
        return headers, col_names, None

    return headers, col_names, np.array(data)


def extract_q_values(header_line):
    """Extract Q values from column header string."""
    if not header_line:
        return []
    matches = re.findall(r'(?:Cl|Csl|Jl|Jls|Sinel|Ssinel|S|Ss|F|Fs)\(\s*([0-9]+\.[0-9]+)', header_line)
    q_vals = []
    for m in matches:
        q = float(m)
        if q not in q_vals:
            q_vals.append(q)
    return q_vals


def analyze_and_plot(data_dir='.'):
    sqw_path = os.path.join(data_dir, 'sqw.dat')
    clqw_path = os.path.join(data_dir, 'clqw.dat')
    fqt_path = os.path.join(data_dir, 'fqt.dat')
    jqt_path = os.path.join(data_dir, 'jqt.dat')

    print("=" * 80)
    print(f"  Consistency Check & PDF Plotting for Dynamic Correlations: {os.path.abspath(data_dir)}")
    print("=" * 80)

    sqw_h, sqw_cols, sqw_d = parse_dat_file(sqw_path)
    clqw_h, clqw_cols, clqw_d = parse_dat_file(clqw_path)
    fqt_h, fqt_cols, fqt_d = parse_dat_file(fqt_path)
    jqt_h, jqt_cols, jqt_d = parse_dat_file(jqt_path)

    if sqw_d is None or clqw_d is None:
        print(f"[ERROR] Could not load required files (sqw.dat: {sqw_d is not None}, clqw.dat: {clqw_d is not None})")
        return

    freq_f = sqw_d[:, 0]
    omega = 2.0 * np.pi * freq_f
    q_clqw = extract_q_values(clqw_h[0] if clqw_h else '')

    print(f"Loaded {len(freq_f)} frequency points from clqw.dat and sqw.dat.")
    if q_clqw:
        print(f"Detected {len(q_clqw)} Q-vectors: {[round(q, 4) for q in q_clqw]}")

    report_lines = []
    report_lines.append("DYNAMIC CORRELATION CONSISTENCY REPORT")
    report_lines.append(f"Directory: {os.path.abspath(data_dir)}")
    report_lines.append("-" * 60)

    # 1. Time-Domain Check
    if jqt_d is not None:
        print("\n--- 1. Time-Domain Initial Value Limit J_{L,s}(Q, t=0) ---")
        report_lines.append("\n1. Time-Domain Initial Value Limit J_{L,s}(Q, t=0):")
        for i_q, q_val in enumerate(q_clqw):
            col_jl = 1 + 2 * i_q
            col_jls = 2 + 2 * i_q
            if col_jls < jqt_d.shape[1]:
                j_l0 = jqt_d[0, col_jl]
                j_ls0 = jqt_d[0, col_jls]
                msg = f"  Q = {q_val:6.3f} : J_L(Q,0) = {j_l0:10.5f} | J_{{L,s}}(Q,0) = {j_ls0:10.5f} (expected ~ k_B*T/m)"
                print(msg)
                report_lines.append(msg)

    # 2. Time Derivative Check
    if fqt_d is not None and jqt_d is not None:
        dt = fqt_d[1, 0] - fqt_d[0, 0]
        print("\n--- 2. Time-Domain Second Derivative -d^2 F(Q,t)/dt^2 vs Q^2 * J_L(Q,t) ---")
        report_lines.append("\n2. Time-Domain Second Derivative:")
        for i_q, q_val in enumerate(q_clqw):
            col_f = 1 + i_q
            col_jl = 1 + 2 * i_q
            if col_f < fqt_d.shape[1] and col_jl < jqt_d.shape[1]:
                f_t = fqt_d[:, col_f]
                d2f_dt2 = (f_t[2:] - 2.0 * f_t[1:-1] + f_t[:-2]) / (dt ** 2)
                calc_jl = - d2f_dt2 / (q_val ** 2)
                jl_actual = jqt_d[1:-1, col_jl]
                mae = np.mean(np.abs(calc_jl - jl_actual))
                msg = f"  Q = {q_val:6.3f} : MAE [-d^2F/dt^2 / Q^2 vs J_L(Q,t)] = {mae:10.5e}"
                print(msg)
                report_lines.append(msg)

    # 3. Frequency Domain Relations & PDF Plotting
    print("\n--- 3. Frequency-Domain Consistency Check & PDF Plot Generation ---")
    report_lines.append("\n3. Frequency-Domain Consistency Checks:")

    # Create Multi-Panel Figure for PDF Export
    n_q = len(q_clqw)
    fig, axes = plt.subplots(n_q, 2, figsize=(12, 4 * n_q), squeeze=False)

    for i_q, q_val in enumerate(q_clqw):
        col_c = 1 + 2 * i_q
        col_cs = 2 + 2 * i_q
        col_s = 1 + 2 * i_q
        col_ss = 2 + 2 * i_q

        if col_cs >= clqw_d.shape[1] or col_ss >= sqw_d.shape[1]:
            continue

        cl_actual = clqw_d[:, col_c]
        csl_actual = clqw_d[:, col_cs]

        s_val = sqw_d[:, col_s]
        ss_val = sqw_d[:, col_ss]

        # Get S(Q) = F(Q,0) scale if needed
        sq_0 = fqt_d[0, 1 + i_q] if fqt_d is not None and (1 + i_q) < fqt_d.shape[1] else 1.0

        # Un-normalized / raw S(Q,w)
        s_raw = s_val * sq_0 if "Sinel" in (sqw_h[0] if sqw_h else "") else s_val

        # Theoretical relations
        cl_calc_q2 = (omega ** 2 / (q_val ** 2)) * s_raw
        cl_calc_q1 = (omega ** 2 / q_val) * s_raw
        csl_calc = (omega ** 2 / (q_val ** 2)) * ss_val

        diff_q2 = np.mean(np.abs(cl_actual - cl_calc_q2))
        diff_q1 = np.mean(np.abs(cl_actual - cl_calc_q1))

        msg = f"  Q = {q_val:6.3f} : Mean Abs Diff [C_L vs w^2*S/Q^2] = {diff_q2:10.5e} | [C_L vs w^2*S/Q] = {diff_q1:10.5e}"
        print(msg)
        report_lines.append(msg)

        # Plot Collective C_L(Q,w) vs w^2*S(Q,w)/Q^2 and w^2*S(Q,w)/Q
        ax_col = axes[i_q, 0]
        ax_col.plot(freq_f, cl_actual, 'b-', lw=2, label=r'Measured $C_L(Q,\omega)$')
        ax_col.plot(freq_f, cl_calc_q2, 'r--', lw=1.5, label=r'$\omega^2 S(Q,\omega) / Q^2$')
        ax_col.plot(freq_f, cl_calc_q1, 'g:', lw=1.5, label=r'$\omega^2 S(Q,\omega) / Q$')
        ax_col.set_xlabel(r'Frequency $\omega / (2\pi)$')
        ax_col.set_ylabel(r'$C_L(Q,\omega)$')
        ax_col.set_title(f'Collective Longitudinal Current Spectrum ($Q = {q_val:.3f}$)')
        ax_col.legend(loc='best')
        ax_col.grid(True, alpha=0.3)
        ax_col.set_xlim(left=0, right=max(freq_f) * 0.3) # Zoom into relevant peak region

        # Plot Self C_{L,s}(Q,w) vs w^2*S_s(Q,w)/Q^2
        ax_self = axes[i_q, 1]
        ax_self.plot(freq_f, csl_actual, 'b-', lw=2, label=r'Measured $C_{L,s}(Q,\omega)$')
        ax_self.plot(freq_f, csl_calc, 'r--', lw=1.5, label=r'$\omega^2 S_s(Q,\omega) / Q^2$')
        ax_self.set_xlabel(r'Frequency $\omega / (2\pi)$')
        ax_self.set_ylabel(r'$C_{L,s}(Q,\omega)$')
        ax_self.set_title(f'Self Longitudinal Current Spectrum ($Q = {q_val:.3f}$)')
        ax_self.legend(loc='best')
        ax_self.grid(True, alpha=0.3)
        ax_self.set_xlim(left=0, right=max(freq_f) * 0.3)

    plt.tight_layout()

    # Save PDF plots (NOT PNG)
    pdf_main_path = os.path.join(data_dir, 'clqw_sqw_comparison.pdf')
    pdf_check_path = os.path.join(data_dir, 'clqw_spectrum_check.pdf')

    fig.savefig(pdf_main_path, format='pdf', bbox_inches='tight')
    fig.savefig(pdf_check_path, format='pdf', bbox_inches='tight')
    print(f"\n[OK] Vector graphics PDF plots saved to:\n  - {pdf_main_path}\n  - {pdf_check_path}")

    # Write report
    report_path = os.path.join(data_dir, 'check_clqw_sqw_report.txt')
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines) + '\n')
    print(f"[OK] Report file written to: {report_path}")

    # Plot on screen
    print("\n[DISPLAY] Rendering plot on screen...")
    try:
        if os.environ.get('DISPLAY'):
            plt.show(block=False)
            plt.pause(1.0)
        else:
            plt.show(block=False)
    except Exception as e:
        print(f"[NOTE] Screen rendering warning: {e}")

    print("=" * 80)


if __name__ == '__main__':
    target_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    analyze_and_plot(target_dir)
