#!/usr/bin/env python3
"""Scenario A forward model: 1-D exponential-decay response.

Reads parameters p1..p8 from model_in.dat, writes 15 simulated
observations to model_out.dat using the closed-form response:

    y_i = sum_{k=1..8}  k * p_k * exp(-i / 10)    for i = 1, 2, ..., 15

This mimics a transport-type response where each parameter contributes
a weighted exponential-decay component over an observation window.
"""
import math
import sys

with open("model_in.dat", "r") as fh:
    lines = [ln for ln in fh.readlines() if ln.strip() and not ln.startswith("#")]

pars = [float(ln.split()[1]) for ln in lines]
nobs = 15

with open("model_out.dat", "w") as fh:
    for i in range(1, nobs + 1):
        val = sum((k + 1) * pk * math.exp(-i / 10.0) for k, pk in enumerate(pars))
        fh.write(f"obs_{i:02d}  {val:.10E}\n")
