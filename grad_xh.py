#!/usr/bin/env python
r"""
sys.argv
--------
1) The name of the output file
2) The reference element to measure the gradients with ('h' for the individual
	[X/H] gradients).
"""

from databanks.SDSS.IPL3 import IPL3_GIANTS
from stats.mode import mode_with_skewnormal_fit
from stats.jackknife import jackknife_summary_statistic
from scipy.stats import linregress
import numpy as np
import vice
import sys
IPL3_GIANTS = vice.dataframe(IPL3_GIANTS).filter(
	"z", ">=", -0.5).filter(
	"z", "<=", 0.5)


RADIAL_BINS = list(range(3, 16))
ELEMENTS = ["al", "ca", "ce", "co", "cr", "o", "k", "s", "na",
	"mg", "mn", "fe", "si", "ni"]


with open(sys.argv[1], "w") as f:
	f.write("# Elem    grad[X/%s] [kpc^-1]    err_grad[X/%s] [kpc^-1]    " % (
		sys.argv[2].capitalize(), sys.argv[2].capitalize()))
	f.write("[X/%s] at R = 8 kpc    err_[X/%s] at R = 8 kpc\n" % (
		sys.argv[2].capitalize(), sys.argv[2].capitalize()))
	radial_bin_centers = []
	xh = dict() # or xfe, but use xh as variable name
	for elem in ELEMENTS: xh[elem] = []
	for i in range(len(RADIAL_BINS) - 1):
		sys.stdout.write("\rR = %g - %g kpc    " % (
			RADIAL_BINS[i], RADIAL_BINS[i + 1]))
		sub = IPL3_GIANTS.filter(
			"rgal", ">=", RADIAL_BINS[i]).filter(
			"rgal", "<=", RADIAL_BINS[i + 1])
		if len(sub["rgal"]) > 200:
			radial_bin_centers.append(
				(RADIAL_BINS[i] + RADIAL_BINS[i + 1]) / 2)
			for elem in ELEMENTS:
				arr = []
				if elem == sys.argv[2]:
					xh[elem].append(0)
					continue
				elif sys.argv[2] == "h":
					for _ in sub["%s_h" % (elem)]:
						if not np.isnan(_): arr.append(_)
				else:
					_xh = sub["%s_h" % (elem)]
					refh = sub["%s_h" % (sys.argv[2])]
					for i in range(len(_xh)):
						if not np.isnan(_xh[i]) and not np.isnan(refh[i]):
							arr.append(_xh[i] - refh[i])
				# xh[elem].append(mode_with_skewnormal_fit(arr))
				# xh[elem].append(np.mean(arr))
				xh[elem].append(np.median(arr))
		else:
			pass
	sys.stdout.write("\n")
	for elem in ELEMENTS:
		reg = linregress(radial_bin_centers, xh[elem])
		f.write("%s\t%.5e\t%.5e\t%.5e\t%.5e\n" % (
			elem.capitalize(),
			reg.slope, reg.stderr,
			reg.intercept, reg.intercept_stderr))
	f.close()
