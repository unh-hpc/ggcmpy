from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def plot_swdata(CDAdata, red_bzneg=False):
    nrows = 8
    ncols = 1
    for i, crd in enumerate("xyz"):
        dat = CDAdata["b%sgse" % crd]
        plt.subplot2grid((nrows, ncols), (0 + i, 0))
        if red_bzneg and crd == "z":
            vals = np.array(dat)
            bzneg = np.ma.masked_where(vals >= 0.0, vals)
            plt.plot(dat.epoch, dat, "b-", dat.epoch, bzneg, "r-")
        else:
            plt.plot(dat.epoch, dat)
        plt.ylabel(f"b{crd}")

    for i, crd in enumerate("xyz"):
        dat = CDAdata["v%sgse" % crd]
        plt.subplot2grid((nrows, ncols), (3 + i, 0))
        plt.plot(dat.epoch, dat)
        plt.ylabel(f"v{crd}")

    dat = CDAdata["np"]
    plt.subplot2grid((nrows, ncols), (6, 0))
    plt.plot(dat.epoch, dat)
    plt.ylabel("np")

    dat = CDAdata["tev"]
    plt.subplot2grid((nrows, ncols), (7, 0))
    plt.plot(dat.epoch, dat)
    plt.ylabel("T (eV)")

    plt.gcf().set_size_inches(15, 15)
    plt.savefig("sw_data.png")
    plt.show()
