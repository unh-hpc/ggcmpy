from __future__ import annotations

try:
    import urllib2  # python 2.x
except ImportError:
    import urllib.request as urllib2  # python 3.x

def get_f107(date, station='Penticton_Observed'):
    """
    Get the f10.7 radio flux given a particular date and station.
    This isn't actually provided by CDAWeb...But whatever.

    inputs:
    date -- either a datetime.datetime object or a time.date object.
            (anything with a year,month and day integer-like attributes)
    station -- Optional station for observation.
               defaults to 'Penticton_Observed'
    """
    year = date.year
    url = ('ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/'
           '%s/%d/%d.OBS'%(station, year, year))
    page = urllib2.urlopen(url)
    for l in page:
        if l.strip().startswith('%d'%date.day):
            values = l.split()[1:]
            return float(values[date.month-1]) / 10.0



if __name__ == "__main__":
    class dt:
        def __init__(self, year, month, day):
            self.year = year
            self.month = month
            self.day = day

    d = dt(2012, 3, 16)  #typically, we would pass datetime objects.
    print(get_f107(d))
