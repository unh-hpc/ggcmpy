"""
Module for parsing CDAWeb data
"""
from __future__ import annotations

import datetime
import os.path

# from orbit import orbit, orbitpoint
import re as _re

import numpy as np


class CDAvar(list):
    def __init__(self, *args, **kwargs):
        list.__init__(self, *args, **kwargs)
        object.__setattr__(self, '_attributes', [])

    def val_at_time(self, time):
        """
        Get the value of this variable at time 'time'
        (time is a datetime.datetime object)
        """

        epochs = [float(x.strftime('%s')) for x in self.epoch]
        itime = float(time.strftime('%s'))
        return np.interp(itime, epochs, self)


    def t_interp(self, gap=30, starttime=None, endtime=None, epoch=None):
        """
        Interpolate in time starting at start, ending at end...
        returns a CDAvar with epochs equally spaced by
        'gap' (in seconds as anything that can cast to an integer).

        starttime and endtime are datetime.datetime objects
        epoch is an iterable of datetime.datetime objects which
              correspond to the elements of CDAvar.
        """
        if epoch is None:
            try:
                epoch = self.epoch
            except:
                raise ValueError('Cannot get epochs')
        #epochs = map(lambda x:float(x.strftime('%s')), epoch)
        epochs = [float(x.strftime('%s')) for x in epoch]
        if starttime is not None:
            st = starttime
            starttime = float(starttime.strftime('%s'))
        else:
            st = self.epoch[0]  #starttime as datetime object
            starttime = epochs[0]
        if endtime is not None:
            endtime = float(endtime.strftime('%s'))
        else:
            endtime = epochs[-1]

        t = starttime
        dt = datetime.timedelta(seconds=int(gap))
        out = []
        outepoch = []
        times = []
        while t <= endtime:
            times.append(t)
            outepoch.append(st)
            st += dt
            t += int(gap)

        out.extend(np.interp(times, epochs, self))

        output = CDAvar(out)

        for a in self._attributes:
            setattr(output, a, getattr(self, a))
        output.epoch = CDAvar(outepoch)

        return output


    def __setattr__(self, attr, value):
        object.__setattr__(self, attr, value)
        if attr not in self._attributes:
            self._attributes.append(attr)

    def __delattr__(self, attr):
        object.__delattr__(self, attr)
        self._attributes.remove(attr)

    def _add_attributes(self, other, out):
        for a in self._attributes:
            setattr(out, a, getattr(self, a))
        if hasattr(other, '_attributes'):
            try:
                for a in other._attributes:
                    setattr(out, a, getattr(other, a))
            except Exception:
                pass


    #Addition
    def __add__(self, other):
        if hasattr(other, '__iter__'):
            out = CDAvar()
            for i, v in enumerate(self):
                out.append(v+other[i])
        else:
            out = CDAvar([v+other for v in self])
        self._add_attributes(other, out)
        return out

    def __radd__(self, other):
        return self.__add__(other)

    #Multiplication
    def __mul__(self, other):
        if hasattr(other, '__iter__'):
            out = CDAvar()
            for i, v in enumerate(self):
                out.append(v*other[i])
        else:
            out = CDAvar([v*other for v in self])
        self._add_attributes(other, out)
        return out

    def __rmul__(self, other):
        return self.__mul__(other)

    #Subtraction
    def __sub__(self, other):
        if hasattr(other, '__iter__'):
            out = CDAvar()
            for i, v in enumerate(self):
                out.append(v-other[i])
        else:
            out = CDAvar([v-other for v in self])
        self._add_attributes(other, out)
        return out

    def __rsub__(self, other):
        if hasattr(other, '__iter__'):
            out = CDAvar()
            for i, v in enumerate(self):
                out.append(other[i]-v)
        else:
            out = CDAvar([other-v for v in self])
        self._add_attributes(other, out)
        return out

    #Division
    def __div__(self, other):
        if hasattr(other, '__iter__'):
            out = CDAvar()
            for i, v in enumerate(self):
                out.append(v.__div__(other[i]))
        else:
            out = CDAvar([v.__div__(other) for v in self])
        self._add_attributes(other, out)
        return out

    def __truediv__(self, other):
        if hasattr(other, '__iter__'):
            out = CDAvar()
            for i, v in enumerate(self):
                out.append(v.__truediv__(other[i]))
        else:
            out = CDAvar([v.__truediv__(other) for v in self])
        self._add_attributes(other, out)
        return out

    def __floordiv__(self, other):
        if hasattr(other, '__iter__'):
            out = CDAvar()
            for i, v in enumerate(self):
                out.append(v.__floordiv__(other[i]))
        else:
            out = CDAvar([v.__floordiv__(other) for v in self])
        self._add_attributes(other, out)
        return out

    def __rdiv__(self, other):
        if hasattr(other, '__iter__'):
            out = CDAvar()
            for i, v in enumerate(self):
                out.append(other[i]/v)
        else:
            out = CDAvar([other/v for v in self])
        self._add_attributes(other, out)
        return out

    #Raising to a power
    def __pow__(self, other):
        if hasattr(other, '__iter__'):
            out = CDAvar()
            for i, v in enumerate(self):
                out.append(v**other[i])
        else:
            out = CDAvar([v**other for v in self])
        self._add_attributes(other, out)
        return out

    def __rpow__(self, other):
        if hasattr(other, '__iter__'):
            out = CDAvar()
            for i, v in enumerate(self):
                out.append(other[i]**v)
        else:
            out = CDAvar([other**v for v in self])
        self._add_attributes(other, out)
        return out

    #It is quite useful to be able to have a different
    #map function...
    def map(self, func):
        out = CDAvar(map(func, self))
        self._add_attributes(None, out)
        return out


######################################################
class CDAWebFile:
    """
    CDAWebFile...looks a lot like a dictionary -- but it isn't.
    (perhaps it should be?)
    """
    bad_data_values = [99.9900,
                       999.990,
                     9999.99,
                     99999.9,
                     -9999.9,
                     1.00000e7,
                     -1.00000E+31,
                     -32768]
    def __init__(self, name=None):
        self.name = name
        if name is not None:
            self.parse()

    def keys(self):
        try:
            return self.attributes
        except AttributeError:
            return []

    def items(self):
        items = []
        for i in self.keys():
            items.append(i, getattr(self, i))
        return items

    def values(self):
        values = []
        for i in self.keys():
            values.append(getattr(self, i))
        return values

    def has_key(self, key):
        return key in self.keys()

    def __contains__(self, key):
        return key in self.keys()

    def name_mangle(self, name):
        ###This might be better done as a dictionary
        ###mapping???
        name = name.lower()
        for var in ('bx', 'by', 'bz',
                    'vx', 'vy', 'vz',
                    'x', 'y', 'z'):
            if var in name and 'gse' in name:
                return var+'gse'

        for var in ('np', 'temp'):
            if var in name and 'na/np' not in name:
                return var

        if 'density' in name:
            return 'np'
        if 'pressure' in name:
            return 'pp'
        if 'ae_index' in name or '1_m_ae' in name:
            return 'AE'
        if 'al_index' in name or '1_m_al' in name:
            return 'AL'
        if 'au_index' in name or '1_m_au' in name:
            return 'AU'
        if name.startswith('sym') and name.endswith('h_index'):
            return 'SYM_H'
        if name.startswith('sym') and name.endswith('d_index'):
            return 'SYM_D'
        if name.startswith('asy') and name.endswith('h_index'):
            return 'ASYM_H'
        if name.startswith('asy') and name.endswith('d_index'):
            return 'ASYM_D'
        if 'dst' in name:
            return 'DST'
        if 'vth' in name:
            return 'vth'
        #geotail
        if 'swa_ion_n' in name:
            return 'np'
        if 'swa_ion_temp' in name:
            return 'temp'
        if 'swa_ion_p' in name:
            return 'pp'


        return name

    _epoch_ids = ('CDF_EPOCH', 'EPOCH_TIME', 'UT', 'TIME_AT_CENTER_OF_HOUR')
    _epoch_ids_regex = _re.compile('|'.join(_re.escape(x) for x in _epoch_ids))
    def parse(self, filename=None):

        if not filename:
            filename = self.name

        with open(filename) as f:
            lines = f.readlines()

        if hasattr(self, 'attributes'):
            _attr = self.attributes[:]
            del self.attributes
        else:
            _attr = []

        while lines:
            l = lines.pop(0)
            l = self._epoch_ids_regex.sub('EPOCH', l)
            # Commented out for testing.  If it works, please remove next line
            #l = l.replace('CDF_EPOCH', 'EPOCH').replace('EPOCH_TIME', 'EPOCH').replace('UT', 'EPOCH').replace('TIME_AT_CENTER_OF_HOUR', 'EPOCH')
            if l.startswith('EPOCH '):
                self.attributes = [self.name_mangle(x.replace('-', '_')) for x in l.split()]
                newlist = []
                for x in self.attributes:
                    if x == 'epoch_time':
                        newlist.append('epoch')
                    else:
                        newlist.append(x)
                self.attributes = newlist

                for a in self.attributes:
                    setattr(self, a, CDAvar())

                l = lines.pop(0)
                while not l.startswith('dd-mm'):
                    l = lines.pop(0)

                ll = self._split_line(l)

                units = {}
                for a, unit in zip(self.attributes, ll, strict=False):
                    units[a] = unit
                break
        else:
            raise ValueError("Cannot parse attribute names from CDAWebFile")

        attr = self.attributes[:]
        attr.remove('epoch')

        while lines:
            try:
                l = lines.pop(0)
                l = self._split_line(l)
                self.epoch.append(datetime.datetime.strptime(l.pop(0), '%d-%m-%Y %H:%M:%S.%f'))
                for a, val in zip(attr, l, strict=False):
                    getattr(self, a).append(float(val))
            except ValueError:
                pass

        for a in attr:
            getattr(self, a).epoch = self.epoch[:]
            setattr(self, a, self.clean_bad_data(getattr(self, a)))
            getattr(self, a).unit = units[a]

        attr = self.attributes
        attr.remove('epoch')

        delattr(self, 'epoch')
        self.attributes = _attr + self.attributes
        for v in 'xyz':
            if '%sgse'%v in self.attributes:
                ll = getattr(self, '%sgse'%v)
                if 'km' in ll.unit:
                    for i in range(len(ll)):
                        ll[i] = ll[i]/6378.0
                    ll.unit = 'RE'

    @staticmethod
    def clean_bad_data(var):
        out = CDAvar()
        out.epoch = CDAvar()
        for i, v in enumerate(var):
            if v not in CDAWebFile.bad_data_values:
                out.append(v)
                out.epoch.append(var.epoch[i])
        return out


    def __getitem__(self, item):
        return getattr(self, item)

    def __setitem__(self, item, value):
        setattr(self, item, value)
        if item not in self.attributes:
            self.attributes.append(item)

    def _split_line(self, line):
        l = line.split()
        output = [' '.join(l[:2])]+l[2:]
        if len(output) == len(self.attributes):
            return output
        raise ValueError("Line Length is wrong:%s"%line)

    # def toOrbit(self, satellite = None):
    #     """
    #     Return an orbit for this satellite
    #     """
    #     for v in 'xyz':
    #         if(not hasattr(self, '%sgse'%v)):
    #             raise ValueError("Object does not have %sgse attribute"%(v))
    #     xgse = self.xgse
    #     ygse = self.ygse
    #     zgse = self.zgse
    #     time = xgse.epoch
    #     out = orbit()
    #     for x, y, z, t in zip(xgse, ygse, zgse, time):
    #         pt = orbitpoint(time = t, point = (x, y, z), satellite = satellite)
    #         out.append(pt)
    #     return out


def makprepinput(filename, field, baddata=None):
    """
    This takes a filename (as attached to a CDAWebFile)
    and turns it into something suitable for mak.prep

    This function is pretty much deprecated as the
    CDAWebRunInput utility provides much more functionality.
    """
    epoch = field.epoch
    with open(filename, 'w') as f:
        for i, v in enumerate(field):
            if v != baddata:
                f.write("%s %s\n"%(epoch[i].strftime("%Y %m %d %H %M %S.%f"), str(v)))


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        try:
            import pyProbe.probeConfigPath as pCP
        except Exception:
            import probeConfigPath as pCP

        pth = [os.path.join(pCP.pyProbePath, 'test', 'AC_H0_MFI_32204.txt')]
    else:
        pth = sys.argv[1:]

    _f = CDAWebFile()
    for _p in pth:
        _f.parse(_p)
    for _a in _f.attributes:
        print(_a, getattr(_f, _a).unit)
    # print(dir(f))
    print(len(_f.bxgse))
    print(len(_f.bxgse.t_interp(30)))
