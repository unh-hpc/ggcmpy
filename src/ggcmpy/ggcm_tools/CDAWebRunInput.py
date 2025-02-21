#!/usr/bin/env python
r"""Matt Gilson's awesome utility to pull solar wind data from CDAWeb

After running this utility, the standard mak.prep utility should be
able to be used on the output, although it may not even be
necessary...

You can always use `CDAWebRunInput --help` to get more info about the
program arguments.
"""

from __future__ import annotations

import argparse
import datetime
import os.path
import sys
from math import sqrt

import xarray as xr

from ggcmpy.ggcm_tools import CDAfetch as fetch
from ggcmpy.ggcm_tools import CDAWeb as cdaweb


def datetimetype(x):
    try:
        return datetime.datetime.strptime(x, "%Y:%m:%d:%H:%M:%S.%f")
    except Exception:
        l = ["%Y", "%m", "%d", "%H", "%M", "%S"]
        for i in reversed(range(len(l) + 1)):
            timefmt = ":".join(l[:i])
            print("Trying time format", timefmt)
            try:
                return datetime.datetime.strptime(x, timefmt)
            except Exception:
                pass
    sys.stderr.write("Parse Error: Cannot parse date from string: %s\n" % x)
    sys.exit(2)


def guess_satellite(parser, fname, debug=False):
    n = os.path.basename(fname).lower()
    if debug:
        print("No satellite specified -- Guessing:")
    for v in ("ac", "wi", "ge", "omni"):
        if n.startswith(v):
            if debug:
                print("\tGuessed:%s" % v)
            return v
    parser.error("Cannot guess satellite")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-d",
        "--debug",
        "--verbose",
        dest="debug",
        action="store_true",
        help="verbose output",
    )
    parser.add_argument(
        "--satellite",
        "--sat",
        "-s",
        dest="sat",
        default=None,
        help="Satellite to use (ac,wi,omni,ge ...)",
    )
    parser.add_argument(
        "--fillgap", dest="gap", type=int, default=30, help="Fill gaps at this cadence"
    )
    parser.add_argument(
        "--starttime",
        "--start",
        "-S",
        dest="starttime",
        type=datetimetype,
        default=None,
        help="Start time (YYYY:mm:dd:HH:MM:Secs)",
    )
    parser.add_argument(
        "--endtime",
        "--end",
        "-E",
        dest="endtime",
        type=datetimetype,
        default=None,
        help="End time (YYYY:mm:dd:HH:MM:Secs)",
    )
    parser.add_argument(
        "--dipoletime",
        "--dip",
        "--diptime",
        "-D",
        dest="diptime",
        type=datetimetype,
        default=None,
        help="dipole time. This is the timecode used for --f107 and --mpos",
    )
    parser.add_argument(
        "--f107",
        dest="f107",
        default=None,
        action="store_true",
        help="Attempt to fetch the solar radio flux from NOAA "
        "(uses a time halfway between startime and endtime, "
        "or --dipoletime if given)",
    )
    parser.add_argument(
        "-m",
        "--monitor-position",
        dest="mopos",
        default=None,
        action="store_true",
        help="Attempt to get monitor position also.",
    )
    parser.add_argument(
        "--no-plot",
        dest="noplot",
        default=False,
        action="store_true",
        help="Don't plot the input variables",
    )
    parser.add_argument(
        "--red_bzneg",
        dest="red_bzneg",
        action="store_true",
        help="When plotting, plot negative Bz values in red",
    )
    parser.add_argument("files", nargs="*", help="CDAWeb files used for input")
    opt = parser.parse_args()
    opt.plot = not opt.noplot

    files = opt.files
    satfetchfuncs = {
        "ac": fetch.getAceData,
        "wi": fetch.getWindData,
        "omni": fetch.getOMNIData,
        "ge": fetch.getGeotailData,
    }

    if not files:  # It looks like we're going to need to download the data
        files = []
        msg = (
            "Must provide start (--starttime), end (--endtime) and satellite "
            "(--satellite) if not providing files!"
        )
        if (not opt.starttime) or (not opt.endtime):
            parser.error(msg)
        if opt.starttime > opt.endtime:
            parser.error("Starttime must precede endtime")
        if opt.debug:
            print("Attempting to fetch necessary data from CDAWeb")
        if not opt.sat:
            parser.error(msg)
        try:
            files.extend(satfetchfuncs[opt.sat](opt.starttime, opt.endtime))
        except KeyError:
            parser.error(
                "Unrecognized satellite [%s] available: %s "
                "" % (opt.sat, str(list(satfetchfuncs.keys())))
            )

    # Files provided, but not satellite name.  Make a guess.
    if not opt.sat:
        opt.sat = guess_satellite(parser, files[0], debug=opt.debug)

    CDAdata = cdaweb.CDAWebFile()
    for f in files:
        if opt.debug:
            print("Attempting to parse file:%s" % f)
        CDAdata.parse(f)

    if opt.debug:
        print("Variables:%s" % (str(CDAdata.keys())))

    # Get the same time-cadence on all the files
    if opt.debug:
        print("Interpolating (interval:%ds)" % (opt.gap))
    for var in CDAdata.keys():
        if opt.debug:
            print("\tvariable:%s" % var)
        CDAdata[var] = CDAdata[var].t_interp(
            gap=opt.gap, starttime=opt.starttime, endtime=opt.endtime
        )

    if opt.sat == "ge":
        for v in "bxgse", "bygse", "bzgse":
            CDAdata[v] = CDAdata[v].map(lambda x: x * 0.1)

    # This portion of the script does the job of amak.derived
    for v in ("np", "rr", "pp", "temp", "tt"):
        if v in CDAdata:
            if opt.debug:
                print("%s --> max(0.001,%s)" % (v, v))
            CDAdata[v] = CDAdata[v].map(lambda x: max(0.001, x))

    if opt.debug and "vth" in CDAdata:
        print("Attempting to calculate temperature from vth")

    if opt.sat == "ac" and "vth" in CDAdata:
        CDAdata["temp"] = CDAdata["vth"].map(lambda x: 60.1366 * x * x * 1.0e-6)
    elif "vth" in CDAdata:
        CDAdata["temp"] = CDAdata["vth"].map(lambda x: 60.1366 * x * x)

    # Skipping currents.  I don't know of any solar wind monitors that
    # Are able to calculate the current ;-)

    # Get the total vector magnitudes if possible
    for v in ("btot", "vtot"):
        if v not in CDAdata:
            if opt.debug:
                print("Attempting to calculate %s" % v)
            keys = [v[0] + x + "gse" for x in ["x", "y", "z"]]
            haskeys = all(map(CDAdata.has_key, keys))
            if not haskeys:
                print("\tFailed.  Does not have all keys:%s" % (str(keys)))
            else:
                CDAdata[v] = (sum([CDAdata[k] ** 2 for k in keys])).map(sqrt)

    if (not CDAdata.has_key("tev")) and CDAdata.has_key("temp"):
        CDAdata["tev"] = CDAdata["temp"] / 11600.0

    if (not CDAdata.has_key("tkev")) and CDAdata.has_key("temp"):
        CDAdata["tkev"] = CDAdata["temp"] / 11600.0e3

    # for whatever reason, amak.derived likes to have a $SC.rr and $SC.np
    if CDAdata.has_key("np") and not CDAdata.has_key("rr"):
        CDAdata["rr"] = CDAdata["np"]

    if CDAdata.has_key("rr") and not CDAdata.has_key("np"):
        CDAdata["np"] = CDAdata["rr"]

    if not CDAdata.has_key("pp"):
        if opt.debug:
            print("Trying to calculate pp from temp and rr")
        if CDAdata.has_key("temp") and CDAdata.has_key("rr"):
            CDAdata["pp"] = 1.381e-5 * CDAdata["rr"] * CDAdata["temp"]
        else:
            print("\tFailed.")
            if not CDAdata.has_key("temp"):
                print("\tNo Temperature")
            else:
                print("\tNo number density")

    def writeout(fname):
        if not CDAdata.has_key(fname):
            return
        field = CDAdata[fname]
        epoch = field.epoch
        filename = ".".join([opt.sat, fname])
        if opt.debug:
            print("Writing:%s" % filename)
        with open(filename, "w") as f:
            for v, t in zip(field, epoch, strict=False):
                try:
                    st = t.strftime("%Y %m %d %H %M %S.%f")
                except Exception:
                    st = t.strftime("%Y %m %d %H %M %S")

                f.write("%s %f\n" % (st, v))

    outputvars = [
        "xgse",
        "ygse",
        "zgse",
        "bxgse",
        "bygse",
        "bzgse",
        "vxgse",
        "vygse",
        "vzgse",
        "pp",
        "rr",
        "np",
        "temp",
        "vth",
        "tkev",
        "tev",
        "btot",
        "vtot",
    ]
    for v in outputvars:
        writeout(v)

    dptime = opt.diptime
    if dptime is None:
        try:
            dptime = opt.starttime + (opt.endtime - opt.starttime) / 2
        except Exception:
            dptime = None

    if opt.f107:
        if dptime is None:
            sys.stderr.write("Cannot guess a suitable time for f107\n")
        else:
            f107 = fetch.get_f107(dptime)
            f = open("f107.input", "w")
            for fobj in f, sys.stdout:
                fobj.write("f107 flux: %f\n" % f107)
            f.close()

    if opt.mopos:
        if dptime is None:
            sys.stderr.write("Cannot guess a suitable time for monitor position\n")
        else:
            for crd in "xyz":
                try:
                    dat = CDAdata["%sgse" % (crd)]
                    s = "MO%s" % (crd.upper())
                    f = open(s, "w")
                    for fobj in f, sys.stdout:
                        fobj.write("%s : %f\n" % (s, dat.val_at_time(dptime)))
                    f.close()
                except Exception as e:
                    sys.stderr.write(str(e) + "\n")

    vars = {
        key: xr.DataArray(list(CDAdata[key]), coords={"time": CDAdata[key].epoch})
        for key in CDAdata.keys()
    }
    sw_data = xr.Dataset(data_vars=vars)

    if opt.plot:
        try:
            from ggcmpy.ggcm_tools import plotting
        except ImportError:
            sys.stderr.write(
                "Tried to plot, but numpy / matplotlib are not "
                "installed; no plot for you :("
            )
            return -1

        plotting.plot_swdata(sw_data, red_bzneg=opt.red_bzneg)

    return 0


if __name__ == "__main__":
    sys.exit(main())

##
## EOF
##
