from urllib.request import urlopen
import os, urllib
from datetime import datetime, timedelta
from io import BytesIO
import numpy as np

urlbase = 'http://ftp.emc.ncep.noaa.gov/gc_wmb/jpeng/eemn/{0:s}/trak.{1:s}.atcfunix.{0}'
atcfdtypes = ('U2', 'U2', 'U10', 'S3', 'U4', 'u1', 'f4', 'f4', 'u1', 'u2', 'S2', 'u1', 'S3', 'u2', 'u2', 'u2', 'u2')
atcfnames = ('basin', 'num', 'time', 'a', 'code', 'fcsthour', 'lat', 'lon', 'wind', 'pres',
    'b', 'rad', 'c', 'r1', 'r2', 'r3', 'r4')
latcvt = lambda x: -int(x[:-1])/10 if x[-1] != 78 else int(x[:-1])/10
loncvt = lambda x: 360-int(x[:-1])/10 if x[-1] != 69 else int(x[:-1])/10
en_num = 25
ep_num = 25

def select(data, bnum, days=None, raw=False):
    if isinstance(data, type(None)):
        return np.array([], dtype=[('lon', np.float32), ('lat', np.float32)])
    basin, num = bnum[:2], bnum[2:4]
    if not raw:
        data = data[['basin', 'num', 'time', 'code', 'fcsthour', 'lat', 'lon', 'wind', 'pres']]
    if days:
        ret = data[np.all((data['basin']==basin, data['num']==num, data['fcsthour']<=24*days), axis=0)]
    else:
        ret = data[np.logical_and(data['basin'] == basin, data['num'] == num)]
    return np.unique(ret)

def load(time, code, prev=0):
    directory = 'tmp/{:s}'.format(time)
    filename = filename_style(code, directory, prev=prev)
    try:
        return np.load(filename)
    except FileNotFoundError:
        print('{:s} - {:s} not found.'.format(time, code))
        return None

def downall(time, prev=0):
    directory = 'tmp/{:s}'.format(time)
    if not os.path.exists(directory):
        os.makedirs(directory)
    en_codes = ['en{:02d}'.format(n) for n in range(1, en_num+1)]
    ep_codes = ['ep{:02d}'.format(n) for n in range(1, ep_num+1)]
    spec_codes = ['eemn', 'emx', 'ec00']
    prev_times = previous_times(time, prev)
    codes_downer(directory, time, en_codes)
    codes_downer(directory, time, ep_codes)
    codes_downer(directory, time, spec_codes)
    for i, t in enumerate(prev_times, 1):
        print(t)
        filename = filename_style('emx', directory, prev=i)
        downer(filename, t, 'emx')

def codes_downer(directory, time, codes):
    for code in codes:
        print(code)
        filename = filename_style(code, directory)
        try:
            downer(filename, time, code)
        except:
            print('**** failed')

def downer(filename, time, code):
    if os.path.exists(filename):
        return
    with urlopen(urlbase.format(time, code), timeout=3) as u:
        text = BytesIO(u.read())
        data = np.genfromtxt(text, dtype=atcfdtypes, names=atcfnames, delimiter=',',
        #usecols=('basin', 'num', 'time', 'code', 'fcsthour', 'lat', 'lon', 'wind', 'pres'),
                        converters={'lat':latcvt, 'lon':loncvt}, autostrip=True)
        np.save(filename, data)

def filename_style(code, directory, prev=0):
    if prev > 0:
        return os.path.join(directory, '{:s}_{:d}.npy'.format(code, prev))
    else:
        return os.path.join(directory, '{:s}.npy'.format(code))

def previous_times(time, prev):
    prevtime = []
    basedt = datetime.strptime(time, '%Y%m%d%H')
    for pt in range(1, 1+prev):
        prevtime.append((basedt - timedelta(hours=12*pt)).strftime('%Y%m%d%H'))
    return prevtime

def prettyprint():
    print('{:=^44s}'.format('PRETTYPRINT'))
    time = input('time>')
    snum = input('storm num>').upper()
    code = input('code (default EMX)>').lower()
    if code[:2] not in ('en', 'ee', 'em', 'ep'):
        print('wrong code. will fall back to emx.')
        code = 'emx'
    data = select(load(time, code), snum)
    basetime = datetime.strptime(time, '%Y%m%d%H')
    print()
    print('-----TIME:{:s}---STORM NUM:{:s}---CODE:{:s}-----'.format(basetime.strftime('%Y/%m/%d %HZ'),
                                                        snum, code.upper()))
    for d in data:
        fcsthour = int(d['fcsthour'])
        lon = round(float(d['lon']), 1)
        lat= round(float(d['lat']), 1)
        pres = int(d['pres'])
        wind = int(d['wind'])
        lonsign = 'E' if lon > 0 else 'W'
        latsign = 'N' if lat > 0 else 'S'
        t = basetime + timedelta(hours=fcsthour)
        print('{:03d}\t{:s}\t{:.1f}{:s}\t{:.1f}{:s}\t{:d} hPa\t{:d} kt'.format(fcsthour, t.strftime('%Y/%m/%d %HZ'),
                                                                      lat, latsign, lon, lonsign, pres, wind))
    print()
    
def main():
    while True:
        i = input('time>')
        if i == 'x':
            prettyprint()
        elif i == '':
            return
        else:
            downall(i)
            return

if __name__ == '__main__':
    #proxy_support = urllib.request.ProxyHandler({'http': '202.112.26.250:8080'})
    #opener = urllib.request.build_opener(proxy_support)
    #urllib.request.install_opener(opener)
    main()
    
