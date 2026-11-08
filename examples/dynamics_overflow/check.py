#!/usr/bin/env python3
"""Check zero-time density and longitudinal-current sum rules."""

import argparse
import pathlib
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("directory", type=pathlib.Path)
parser.add_argument("--temperature", type=float, default=1.5)
parser.add_argument("--require-overflow", action="store_true")
args = parser.parse_args()
directory = args.directory

def read(name):
    return np.loadtxt(directory / name, comments="#", ndmin=2)

f = read("fqt.dat")
fs = read("fskt.dat")
j = read("jqt.dat")
sq = read("sq.dat")
headers = (directory / "fskt.dat").read_text().splitlines()[0]
q_values = np.array([float(x) for x in re.findall(r"Fs\(\s*([\d.]+),t\)", headers)])
if not len(q_values):
    raise SystemExit("No Q labels in fskt.dat")

failed = False
print("Q       Fs(0)     F(0)      S_NN      J_L(0)   J_L,s(0)  -F''(0)/Q²")
for k, q in enumerate(q_values):
    static = sq[np.argmin(abs(sq[:, 0] - q)), 1]
    f0, fs0 = f[0, k + 1], fs[0, k + 1]
    jl, jls = j[0, 2 * k + 1:2 * k + 3]
    coefficients = np.polyfit(f[:5, 0], f[:5, k + 1], 4)
    curvature = -2 * coefficients[-3] / q**2
    print(f"{q:5.2f}  {fs0:9.5f}  {f0:9.5f}  {static:9.5f}  {jl:9.5f}  {jls:9.5f}  {curvature:9.5f}")
    if abs(fs0 - 1) > 1e-5:
        failed = True
    if abs(f0 - static) > 0.02 * abs(static):
        failed = True
    if max(abs(jl - args.temperature), abs(jls - args.temperature)) > 0.005 * args.temperature:
        failed = True

if args.require_overflow:
    output = (directory / "run.log").read_text()
    match = re.search(r"(\d+) time origins", output)
    if not match:
        raise SystemExit("Cannot find origin count in run.log")
    origins = int(match.group(1))
    natoms = 256
    largest_nq = int(sq[np.argmin(abs(sq[:, 0] - q_values[-1])), -1])
    product = natoms * origins * largest_nq
    print(f"Largest integer normalization product: {product}")
    if product <= 2**31 - 1:
        failed = True

if failed:
    raise SystemExit("Zero-time dynamics regression failed")
print("Zero-time dynamics regression passed")
