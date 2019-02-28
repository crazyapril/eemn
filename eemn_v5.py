import datetime
import ftplib
import os
from mpkit import fastdown, plotplus
from urllib.request import urlopen

import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.patch import geos_to_path
from matplotlib.collections import PathCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, PathPatch
from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
from pybufrkit.dataquery import DataQuerent, NodePathParser
from pybufrkit.decoder import Decoder
from shapely.geometry import LineString, MultiLineString, Point, box
from shapely.ops import cascaded_union

#EEMN Show tracks and strike probabilities of tropical cyclones based on ECMWF ensemble
#---Configuration
RES = 0.1 #spatial resolution of probabilities
RADIUS = 1 #size of strike swath, unit in lat/lon
THRESHOLD = 10 #zoom in to only show area where probabilities are greater than THRESHOLD
SHOW_PROB = 20 #only show cities where probabilities are greater than SHOW_PROB
NUM_OF_MEMS = 50 #num of ensemble members
DAYS_AHEAD = None #only draw tracks during next few days
NO_TRACK = False #no ensemble track display
NO_SUBPLOT = True #no subplot
REALTIME_POS = True #show realtime position (using JTWC fix)
TOP_CITIES = 9 #Maxmimum num of cities to be shown
ADVANCE_LINES = False #(DISABLED) Advance line based on tracks
ONLY_EXISTING_STORM = True
INVALID_MEAN_THRESHOLD = 15
INVALID_MAJORITY_THRESHOLD = 45

#---What's new
#2.0# Upgraded from EEMN v1 / Graph redesigned
#2.1# Strike probabilities added
#2.2# Intensity column added
#2.3# Probabilities of cities added
#2.4# Support for multiple storms / Improved user interface
#3.0# Shapely implement / Coastline subplot / Ensemble control / Days limit
#3.1# Wind radii / Colorbar relocated / Improved intensity column
#3.2# New colormap / Add NO_TRACK flag / Realtime fix display
#4.0# Class oriented rewrite
#5.0#

necoast = r'coastlines\10m_coastline'
nrlsector_url = 'http://tropic.ssec.wisc.edu/real-time/amsu/herndon/new_sector_file'

np.warnings.filterwarnings('ignore')


def currentstorm():
    sectors = urlopen(nrlsector_url)
    cstorms = []
    line = sectors.readline().decode(encoding='ascii')
    while line:
        sstorm = []
        data = line.split()
        sstorm.append(data[0]) #ATCFNUM
        sstorm.append(data[1]) #NAME
        sstorm.append(data[2]+data[3]) #TIME
        lat, lon = eval(data[4][:-1]), eval(data[5][:-1])
        if data[4][-1] == 'S':
            lat = -lat
        if data[5][-1] == 'W':
            lon = 360-lon
        sstorm.append(lat) #LATITUDE
        sstorm.append(lon) #LONGTITUDE
        sstorm.append(data[6]) #BASIN
        sstorm.append(int(data[7])) #WIND
        sstorm.append(int(data[8])) #PRESSURE
        line = sectors.readline().decode(encoding='ascii')
        cstorms.append(sstorm)
    sectors.close()
    return cstorms


class ECMWFDown:

    def __init__(self):
        self.connected = False
        self.storms = []
        self.existing_storms = []

    def connect(self):
        self.ftp = ftplib.FTP('data-portal.ecmwf.int', user='wmo', passwd='essential')
        self.connected = True

    def set_time(self, time):
        self.basetime = time
        self.listfile = 'tmp/{}/storms.txt'.format(time)

    def search(self):
        if not os.path.exists(self.listfile):
            self.search_from_ftp()
        else:
            with open(self.listfile) as f:
                for line in f:
                    filename = line.split(',')[-1]
                    bf = BufrFile(filename)
                    self.storms.append(bf)
                    if bf.emx_flag != 'X':
                        self.existing_storms.append(bf)
        if ONLY_EXISTING_STORM:
            return self.existing_storms
        else:
            return self.storms

    def search_from_ftp(self):
        self.connect()
        self.ftp.cwd(self.basetime + '0000')
        filenames = self.ftp.nlst()
        for fname in filenames:
            if 'tropical_cyclone_track' in fname and 'ECEP' in fname:
                bf = BufrFile(fname)
                self.storms.append(bf)
                self.existing_storms.append(bf)
        self.write_file()

    def write_file(self):
        os.makedirs(os.path.dirname(self.listfile), exist_ok=True)
        with open(self.listfile, 'w') as f:
            for storm in self.storms:
                line = ','.join([
                    storm.emx_flag,
                    str(storm.num),
                    storm.codename,
                    '{:.1f}'.format(storm.slon),
                    '{:.1f}'.format(storm.slat),
                    storm.filename
                ])
                f.write(line+'\n')
    
    def download(self, storms=None):
        if storms is None:
            storms = self.storms
        downlist = [storm for storm in storms if not storm.file_exists()]
        print(downlist)
        if not downlist:
            return storms
        if not self.connected:
            self.connect()
        dirpath = '/'+downlist[0].basetime+'0000'
        if self.ftp.pwd() != dirpath:
            self.ftp.cwd(dirpath)
        downer = fastdown.FTPFastDown(file_parallel=2)
        downer.set_ftp(self.ftp)
        downer.set_task([(s.filename, s.filepath) for s in downlist])
        downer.download()
        return storms


class BufrFile:

    CODE_LAT = '005002'
    CODE_LON = '006002'
    CODE_WIND = '011012'
    CODE_PRES = '010051'
    CODE_TIME = '004024'

    def __init__(self, filename):
        self.filename = filename
        self.loaded = False
        self._analyze_filename()

    def __repr__(self):
        return '<{}>'.format(self.codename)

    def file_exists(self):
        return os.path.exists(self.filepath)
        
    def _analyze_filename(self):
        segs = self.filename.split('_')
        self.emx_flag = 'X' if 'ECEP' not in segs[1] else 'E'
        self.num = int(segs[1][4:6])
        self.basetime = segs[4][:10]
        self.codename = segs[8]
        self.atcfname = None
        self.slon = float(segs[9][:-4].replace('p', '.'))
        self.slat = float(segs[10][:-4].replace('p', '.'))
        self.filepath = 'tmp/{}/{}.bufr'.format(self.basetime, self.codename)

    def get_fullname(self):
        if self.atcfname:
            return '{} ({})'.format(self.atcfname, self.codename)
        else:
            return '[{}]'.format(self.codename)

    def load(self):
        if self.loaded:
            return
        with open(self.filepath, 'rb') as f:
            message = Decoder().process(f.read())
        queryer = DataQuerent(NodePathParser())
        self._lons = []
        self._lats = []
        self._wind = []
        self._pres = []
        for subset in range(52):
            # lat
            try:
                values = queryer.query(message, '@[{}] > {}'.format(subset,
                    self.CODE_LAT)).all_values()
            except IndexError:
                raw_lats = np.empty(41)
                raw_lats[:] = np.nan
            else:
                raw_lats = np.array(values[0][3], dtype='float')[:,0]
                raw_lats = np.insert(raw_lats, 0, values[0][1])
            self._lats.append(raw_lats)
            # lon
            try:
                values = queryer.query(message, '@[{}] > {}'.format(subset,
                    self.CODE_LON)).all_values()
            except IndexError:
                raw_lons = np.empty(41)
                raw_lons[:] = np.nan
            else:
                raw_lons = np.array(values[0][3], dtype='float')[:,0]
                raw_lons = np.insert(raw_lons, 0, values[0][1])
            raw_lons[raw_lons<0] = raw_lons[raw_lons<0] + 360
            self._lons.append(raw_lons)
            # wind
            try:
                values = queryer.query(message, '@[{}] > {}'.format(subset,
                    self.CODE_WIND)).all_values(flat=True)
            except IndexError:
                raw_wind = np.empty(41)
                raw_wind[:] = np.nan
            else:
                raw_wind = np.array(values[0], dtype='float') * 1.94 # to kt
            self._wind.append(raw_wind)
            # pres
            try:
                values = queryer.query(message, '@[{}] > {}'.format(subset,
                    self.CODE_PRES)).all_values(flat=True)
            except IndexError:
                raw_pres = np.empty(41)
                raw_pres[:] = np.nan
            else:
                raw_pres = np.array(values[0], dtype='float') / 100 # to hPa
            self._pres.append(raw_pres)
        self.invalid_indices = []
        self.invalid_majors = []
        self._lats = self.compact_mean(self._lats)
        self._lons = self.compact_mean(self._lons)
        self._wind = self.compact_mean(self._wind)
        self._pres = self.compact_mean(self._pres)
        invalid_index = min(self.invalid_indices)
        invalid_major = min(self.invalid_majors)
        print(invalid_index, invalid_major)
        self.cut_major(self._lats, invalid_major)
        self.cut_major(self._lons, invalid_major)
        self.cut_major(self._wind, invalid_major)
        self.cut_major(self._pres, invalid_major)
        self._lats[-1, invalid_index:] = np.nan
        self._lons[-1, invalid_index:] = np.nan
        self._wind[-1, invalid_index:] = np.nan
        self._pres[-1, invalid_index:] = np.nan
        self._maxwind = np.nanmax(self._wind, axis=1)
        self._minpres = np.nanmin(self._pres, axis=1)
        #print(self._maxwind)
        #print(self._minpres)
        self.loaded = True

    def compact_mean(self, arr):
        arr = np.vstack(arr)
        means = np.nanmean(arr, axis=0)
        valid_counts = np.isnan(arr).astype(int).sum(axis=0)
        arr = np.vstack((arr, means))
        self.invalid_indices.append(self.calc_valid_index(valid_counts, INVALID_MEAN_THRESHOLD))
        self.invalid_majors.append(self.calc_valid_index(valid_counts, INVALID_MAJORITY_THRESHOLD))
        return arr

    def calc_valid_index(self, series, threshold, gap=1):
        if series[0] <= threshold:
            # Newly formed typhoon?
            first_valid_index = np.argmax(series>threshold)
            if first_valid_index == 0:
                # Never ever above the threshold -> ALL VALID
                valid_index = len(series)
            else:
                valid_index = np.argmax(series[first_valid_index:]<threshold-gap)
                #valid_index = len(series)
        else:
            valid_index = np.argmax(series<threshold)
        return valid_index

    def cut_major(self, arr, invalid_major):
        for i in range(52):
            a = arr[i, :]
            invalid_index_0 = np.argmax(~np.isnan(a))
            invalid_index = np.argmax(np.isnan(a[invalid_index_0:]))
            invalid_cut = max(invalid_index, invalid_major)
            a[invalid_cut:] = np.nan

    def set_hour_range(self, hours):
        index = hours // 6 + 1
        self._lats[:, index:] = np.nan
        self._lons[:, index:] = np.nan
        self._wind[:, index:] = np.nan
        self._pres[:, index:] = np.nan
        self._maxwind = np.nanmax(self._wind, axis=1)
        self._minpres = np.nanmin(self._pres, axis=1)

    def iter_members(self):
        for i in range(50):
            mask = np.isnan(self._lats[i, :]) | np.isnan(self._lons[i, :]) | \
                np.isnan(self._wind[i, :]) | np.isnan(self._pres[i, :])
            self.lats = self._lats[i, :][~mask]
            self.lons = self._lons[i, :][~mask]
            self.wind = self._wind[i, :][~mask]
            self.pres = self._pres[i, :][~mask]
            try:
                self.maxwind = self.wind.max()
            except ValueError:
                self.maxwind = None
            try:
                self.minpres = self.pres.min()
            except ValueError:
                self.minpres = None
            if i < 25:
                code = 'EN{:02d}'.format(i + 1)
            else:
                code = 'EP{:02d}'.format(i - 24)
            yield code

    def set_data_pointer(self, code):
        if code == 'EC00':
            i = 50
        elif code == 'EMX':
            i = 51
        elif code == 'EEMN':
            i = 52
        elif code.startswith('EN'):
            i = int(code[2:]) - 1
        elif code.startswith('EP'):
            i = int(code[2:]) + 24
        mask = np.isnan(self._lats[i, :]) | np.isnan(self._lons[i, :]) | \
                np.isnan(self._wind[i, :]) | np.isnan(self._pres[i, :])
        self.lats = self._lats[i, :][~mask]
        self.lons = self._lons[i, :][~mask]
        self.wind = self._wind[i, :][~mask]
        self.pres = self._pres[i, :][~mask]
        try:
            self.maxwind = self.wind.max()
        except ValueError:
            self.maxwind = None
        try:
            self.minpres = self.pres.min()
        except ValueError:
            self.minpres = None

    def get_georange(self):
        latmax = np.nanmax(self._lats)
        latmin = np.nanmin(self._lats)
        lonmax = np.nanmax(self._lons)
        lonmin = np.nanmin(self._lons)
        return latmin, latmax, lonmin, lonmax


def geoscale(latmin, latmax, lonmin, lonmax, scale=0.8, pad=0.):
    latmid = latmax/2 + latmin/2
    lonmid = lonmax/2 + lonmin/2
    deltalat = latmax - latmin + 2 * pad
    deltalon = lonmax - lonmin + 2 * pad
    if deltalat / deltalon > scale:
        deltalon = deltalat / scale
        lonmax = lonmid + deltalon / 2
        lonmin = lonmid - deltalon / 2
    elif deltalat / deltalon < scale:
        deltalat = deltalon * scale
        latmax = latmid + deltalat / 2
        latmin = latmid - deltalat / 2
    return latmin, latmax, lonmin, lonmax

def roundit(georange):
    latmin, latmax, lonmin, lonmax = georange
    latmin = round(latmin / RES) * RES
    latmax = round(latmax / RES) * RES
    lonmin = round(lonmin / RES) * RES
    lonmax = round(lonmax / RES) * RES
    georange = latmin, latmax, lonmin, lonmax
    return georange

def get_grids(georange):
    latmin, latmax, lonmin, lonmax = georange
    x = np.arange(lonmin, lonmax+RES, RES)
    y = np.arange(latmin, latmax+RES, RES)
    grid = np.zeros((y.shape[0], x.shape[0]))
    return x, y, grid

def getcolor(i):
    if i == None:
        return 'X', '#AAAAAA'
    s = '%d hPa' % i
    if i > 1000:
        return s, '#000000'
    if i > 990:
        return s, '#2288FF'
    if i > 970:
        return s, 'orange'
    if i > 950:
        return s, '#FF2288'
    return s, '#800000'

def a_color(p):
    txtcolor = 'k' if p < 35 else 'w'
    p = 100 - p
    p /= 100
    bgcolor = plt.cm.hot(p)
    return bgcolor, txtcolor

class EEMN:
    code_eps = ['ep%02d' % i for i in range(1, 26)]
    code_ens = ['en%02d' % i for i in range(1, 26)]
    code_sps = ['emx', 'eemn', 'ec00']

    def __init__(self, storms):
        self.storms = storms
        self.MULTIFLAG = len(self.storms) > 1

    def run(self):
        self.init()
        self.analyze()
        self.set_map()
        self.plot_probs()
        self.plot_tracks()
        self.plot_legend()
        self.plot_infos()
        if not self.MULTIFLAG:
            self.plot_city_probs()
            self.plot_columns()
            # if REALTIME_POS:
            #     self.plot_realtime_position()
            if not NO_SUBPLOT:
                self.plot_subplot()
        self.save()

    def init(self):
        ###LOAD TRACKS / CALCULATE RANGE
        print('...CALCULATING FIGURE RANGE...')
        geos = []
        for storm in self.storms:
            storm.load()
            if DAYS_AHEAD is not None:
                storm.set_hour_range(int(DAYS_AHEAD * 24))
            geos.append(storm.get_georange())
        georange = self.merge_georanges(geos)
        latmin, latmax, lonmin, lonmax = georange
        self.georange = roundit(geoscale(latmin, latmax, lonmin, lonmax))
        print(self.georange)

    def merge_georanges(self, geos):
        pad = 1
        lats = []
        lons = []
        for georange in geos:
            lats.append(georange[0])
            lats.append(georange[1])
            lons.append(georange[2])
            lons.append(georange[3])
        return min(lats)-1, max(lats)+1, min(lons)-1, max(lons)+1

    def analyze(self):
        ###ANALYZE TRACKS
        print('...ANALYZING STORM TRACKS...')
        allgrids = list()
        x, y, self.grid = get_grids(self.georange)
        self.xshape = x.shape[0]
        self.yshape = y.shape[0]
        xx, yy = np.meshgrid(x, y)
        self.xy = np.dstack((xx, yy)).reshape((self.xshape * self.yshape, 2))
        for storm in self.storms:
            storm_grid = self.add_storm(storm)
            storm_grid = storm_grid * 100 / NUM_OF_MEMS
            allgrids.append(storm_grid)
        self.prob_grid = np.amax(allgrids, axis=0)
        #probs = calc_probs(grid, georange)
        self.cities_probs = self.calc_city_probs()
        self.calc_new_georange()

    def calc_city_probs(self):
        cities = list()
        latmin, latmax, lonmin, lonmax = self.georange
        f = open('cities.txt', 'r', encoding='utf-8')
        for line in f:
            data = line.split()
            name = data[1]
            lat = round(float(data[2][:-1]) / RES) * RES
            if lat > latmax or lat < latmin:
                continue
            lon = round(float(data[3][:-1]) / RES) * RES
            if lon > lonmax or lon < lonmin:
                continue
            x = int((lon - lonmin) / RES)
            y = int((lat - latmin) / RES)
            p = self.prob_grid[y, x]
            if p > SHOW_PROB:
                cities.append([p, name])
        cities.sort(reverse=True)
        print(cities)
        f.close()
        if len(cities) > TOP_CITIES:
            cities = cities[:TOP_CITIES]
        return cities

    def set_map(self):
        ###PLOT: SET MAP
        print('...PLOTING...')
        self.p = plotplus.Plot()
        #self.p.setfamily('Segoe UI Emoji')
        self.p.setmap(projection='cyl', georange=self.ngeorange, resolution='i')
        self.p.setxy(self.georange, RES)

    def plot_probs(self):
        ###PLOT: PLOT PROBABILITIES & COLORBAR
        print('...PLOTING CONTOURF...')
        self.p.contourf(self.prob_grid, gpfcmap='strikeprob', levels=np.arange(0, 101, 2),
                        cbar=True, cbardict=dict(sidebar=True))

    def plot_tracks(self):
        ###PLOT: PLOT LINES
        print('...PLOTING LINES...')
        for storm in self.storms:
            self.intens = []
            # Deterministic
            storm.set_data_pointer('EMX')
            self.p.plot(storm.lons, storm.lats, marker='o', markersize=2, mec='none',
                        linestyle='-', lw=0.5, color='#8877CC')
            self.intens.append(('DET', storm.minpres))
            # Mean
            storm.set_data_pointer('EEMN')
            self.p.plot(storm.lons, storm.lats, marker='o', markersize=2, mec='none',
                        linestyle='-', lw=0.5, color='#99DD22')
            self.intens.append(('MEAN', storm.minpres))
            # Control
            storm.set_data_pointer('EC00')
            self.p.plot(storm.lons, storm.lats, marker='o', markersize=2, mec='none',
                        linestyle='-', lw=0.5, color='#AAAAAA')
            self.intens.append(('CTRL', storm.minpres))
            for code in storm.iter_members():
                if not NO_TRACK:
                    self.p.plot(storm.lons, storm.lats, marker=None, linestyle='-', lw=0.3,
                                color='#CCCCCC')
                self.intens.append((code, storm.minpres))

    def plot_legend(self):
        ###PLOT: PLOT LEGEND
        h_e = Line2D([], [], color='#CCCCCC', lw=0.3, marker=None, label='Ensemble Cluster')
        h_x = Line2D([], [], color='#8877CC', lw=0.5, marker='o', ms=2, mec='none', label='Deterministic')
        h_n = Line2D([], [], color='#99DD22', lw=0.5, marker='o', ms=2, mec='none', label='Ensemble Mean')
        h_c = Line2D([], [], color='#AAAAAA', lw=0.5, marker='o', ms=2, mec='none', label='Ensemble Control')
        handles = [h_e, h_x, h_n, h_c]
        self.p.legend(handles=handles, loc='upper right', framealpha=0.8)

    def plot_infos(self):
        ###PLOT: PLOT INFORMATION
        namestr = ' & '.join([storm.get_fullname() for storm in self.storms])
        self.time = self.storms[0].basetime
        hourstr = '(Within {:d} hours)'.format(int(DAYS_AHEAD * 24)) if DAYS_AHEAD else ''
        self.p.title('Strike Probabilites* of %s Based on ECMWF Ensemble %s' % (namestr, hourstr))
        self.p._timestamp('Init Time: {:s}/{:s}/{:s} {:s}Z'.format(self.time[:4], self.time[4:6],
                                                                   self.time[6:8], self.time[8:]))
        self.p.draw('meripara country province coastline')
        self.p._maxminnote('*probability that the center of the tropical cyclone will pass'
                           ' within 1 lat/lon (approx. 100~110km) of a location')

    def plot_city_probs(self):
        ###PLOT: PLOT CITY PROBABILITIES
        x = 0.02
        for item in self.cities_probs:
            prob, name = tuple(item)
            bgcolor, txtcolor = a_color(prob)
            y = -0.04
            s = '{:s}  {:.0f}%'.format(name, prob)
            a = self.p.ax.annotate(s, xy=(x, y), va='top', ha='left', xycoords='axes fraction',
                                   fontsize=6, family='Source Han Sans CN', color=txtcolor,
                                   bbox=dict(facecolor=bgcolor, edgecolor='none',
                                             boxstyle='square', alpha=0.6))
            self.p.ax.figure.canvas.draw()
            x = self.p.ax.transAxes.inverted().transform(a.get_window_extent())[1, 0] + 0.02

    def plot_columns(self):
        ###PLOT: PLOT MEMBER PRESSURE
        for i, e in enumerate(self.intens[::-1]):
            code, inten = e
            s, c = getcolor(inten)
            s = code + ' ' + s
            self.p.ax.text(1.01, i * 0.02, s, color=c, fontsize=5, family='Lato',
                           transform=self.p.ax.transAxes)

    def plot_realtime_position(self):
        basin_codes = {
            'WP': 'W',
            'EP': 'E',
            'AL': 'L',
            'SI': 'S',
            'SP': 'P'
        }
        bnum = self.bnums[0]
        current_storms = currentstorm()
        transformed_bnum = bnum[2:] + basin_codes[bnum[:2]]
        for storm in current_storms:
            if transformed_bnum == storm[0]:
                rt_time = '{}/{}Z'.format(storm[2][4:6], storm[2][6:8])
                rt_lat = storm[3]
                rt_lon = storm[4]
                rt_intensity = storm[6]
                print(rt_lat, rt_lon)
                rt_str = '{}\n{}kt'.format(rt_time, rt_intensity)
                self.p.marktext(rt_lon, rt_lat, rt_str, stroke=True)
                print('Realtime fix plotted')

    def plot_subplot(self):
        if len(self.cities_probs) == 0:
            return
        ###SUBPLOT: GET HIGHLIGHT COASTLINES
        print('...PLOTING SUBPLOT...')
        self.p.m.readshapefile(necoast, 'necoast', linewidth=0.3, color='#222222')
        highlights, cgeorange = self.get_highlight_coastline(self.p.m.necoast)
        if len(highlights[1]) == 0:
            return
        self.buffersize = (cgeorange[3] - cgeorange[2]) / 100
        cgeorange = geoscale(*cgeorange, scale=0.43)
        self.set_subplot_map(highlights, cgeorange)
        self.plot_subplot_track(highlights)

    def save(self):
        dpi = 250 if self.MULTIFLAG else 250
        self.p.setdpi(dpi)
        self.p.save('Z.png')

    def set_subplot_map(self, highlights, cgeorange):
        ###SUBPLOT: SET MAP
        self.oldax = self.p.ax
        self.p.ax = self.p.fig.add_axes([0.02,-0.39,1.04,0.44])
        self.p.setmap(projection='cyl', georange=cgeorange, resolution='i')
        self.p.draw('country coastline province city')
        ###SUBPLOT: HIGHLIGHT COASTLINES & PLOT LEGEND
        colors = ['#AAFFAA', '#FFFF44', '#FF3333', '#BB0044']
        descr = ['10~25%', '25~50%', '50~75%', '>75%']
        handles = []
        for i, clr in enumerate(colors, 0):
            patch = PathCollection(geos_to_path(MultiLineString(highlights[i]).buffer(self.buffersize)),
                                   facecolor=clr)
            self.p.ax.add_collection(patch)
            handles.append(Patch(color=clr, label=descr[i]))
        self.p.ax.text(0.98, 0.27, '中心经过1经纬度\n范围内的几率', transform=self.p.ax.transAxes,
                       va='bottom', ha='right', fontsize=6, family='Source Han Sans CN')
        self.p.legend(handles=handles, loc='lower right', framealpha=0.8)

    def plot_subplot_track(self, highlights):
        ###SUBPLOT: PLOT DETERMINISTIC TRACK
        storm = self.storms[0]
        storm.set_data_pointer('EMX')
        xlon, xlat = storm.lons, storm.lats
        self.p.plot(xlon, xlat, marker='o', markersize=2, mec='none', linestyle='-',
                    lw=0.5, color='#CCCCCC')
        self.p.ax = self.oldax

    def get_highlight_coastline(self, lines):
        latmin, latmax, lonmin, lonmax = self.georange
        highlights = [[], [], [], []] #10~25% 25~50% 50~75% 75~100%
        latmins, latmaxs, lonmins, lonmaxs = [], [], [], [] #Coastline georange
        if lonmin > 180:
            lonmin = lonmin - 360
        if lonmax > 180:
            lonmax = lonmax - 360
        boundbox = box(lonmin, latmin, lonmax, latmax)
        segs = []
        for seg in lines:
            ls = LineString(seg).intersection(boundbox)
            if isinstance(ls, LineString):
                segs.append(np.array(ls))
            elif isinstance(ls, MultiLineString):
                segs.extend([np.array(s) for s in ls.geoms])
        for s in segs:
            sr = np.round(s / RES) * RES
            xi = ((sr[:,0] - lonmin) / RES).astype(np.int)
            yi = ((sr[:,1] - latmin) / RES).astype(np.int)
            p = self.prob_grid[yi, xi].astype(np.uint8)
            if p.max() < 10:
                continue
            p[p < 10] = 0
            p[(p >= 10) & (p < 25)] = 1
            p[(p >= 25) & (p < 50)] = 2
            p[(p >= 50) & (p < 75)] = 3
            p[p >= 75] = 4
            cutindex = np.where(np.diff(p))[0] + 1
            cutseg = np.split(s, cutindex)
            cutindex = np.insert(cutindex, 0, 0)
            for cseg, cindex in zip(cutseg, cutindex):
                if p[cindex] > 0 and len(cseg) > 1:
                    if p[cindex] > 1:
                        latmins.append(cseg[:,1].min())
                        latmaxs.append(cseg[:,1].max())
                        lonmins.append(cseg[:,0].min())
                        lonmaxs.append(cseg[:,0].max())
                    highlights[p[cindex]-1].append(tuple(map(tuple, cseg)))
        if len(latmins) == 0:
            return None, None
        cgeorange = min(latmins), max(latmaxs), min(lonmins), max(lonmaxs)
        return highlights, cgeorange

    def add_storm(self, storm):
        count = self.grid.copy()
        length = NUM_OF_MEMS
        name = storm.get_fullname()
        i = 0
        for code in storm.iter_members():
            linestring = list(zip(storm.lons, storm.lats))
            if len(linestring) < 2:
                continue
            path = geos_to_path(LineString(linestring).buffer(RADIUS))[0]
            boolarr = path.contains_points(self.xy).reshape((self.yshape, self.xshape)).astype(np.uint8)
            count += boolarr
            i += 1
            print('\r{:s}: [{: <10s}] {:.0%}'.format(name, '#'*(i*10//length), i/length), end='')
        print()
        return count

    def calc_new_georange(self, pad=2):
        y, x = np.where(self.prob_grid > THRESHOLD)
        latmin, latmax, lonmin, lonmax = self.georange
        nlatmin = latmin + y.min() * RES - pad
        nlatmax = latmin + y.max() * RES + pad
        nlonmin = lonmin + x.min() * RES - pad
        nlonmax = lonmin + x.max() * RES + pad
        #print(latmin, latmax, nlatmin, nlatmax)
        self.ngeorange = geoscale(nlatmin, nlatmax, nlonmin, nlonmax)


class UserInterface:

    def interface(self):
        self.input_basetime()
        self.get_storm()
        self.input_storm()
        self.download_selected()
        self.plot()

    def input_basetime(self):
        while True:
            time_input = input('Time>')
            if time_input == 'hack':
                return self.hack_mode()
            elif time_input.isdigit() and len(time_input) == 10:
                self.set_basetime(time_input)
                break

    def set_basetime(self, basetime):
        self.basetime = basetime
        self.basetime_dt = datetime.datetime.strptime(self.basetime, '%Y%m%d%H')

    def _repr_lat_lon(self, lat, lon):
        if lat >= 0:
            lat_repr = '{:.1f}N'.format(lat)
        else:
            lat_repr = '{:.1f}S'.format(-lat)
        if lon >= 0:
            lon_repr = '{:.1f}E'.format(lon)
        else:
            lon_repr = '{:.1f}W'.format(-lon)
        return lat_repr, lon_repr

    def get_currentstorms(self):
        self.currentstorms = []
        if datetime.datetime.now() - self.basetime_dt > datetime.timedelta(days=2):
            return
        try:
            self.currentstorms = currentstorm()
        except Exception:
            raise

    def get_atcfname(self):
        for storm in self.storms:
            candidates = []
            for active_storm in self.currentstorms:
                alat = active_storm[3]
                alon = active_storm[4]
                latdelta = abs(alat - storm.slat)
                londelta = abs(alon - storm.slon)
                if latdelta > 5 or londelta > 5:
                    continue
                candidates.append((latdelta**2+londelta**2, active_storm[0], active_storm[1]))
            if candidates:
                candidates.sort()
                storm.atcfname = candidates[0][1]
                if storm.codename[0].isdigit():
                    storm.codename = candidates[0][2]

    def get_storm(self):
        os.makedirs('tmp/{}'.format(self.basetime), exist_ok=True)
        self.downer = ECMWFDown()
        self.downer.set_time(self.basetime)
        self.storms = self.downer.search()
        #self.get_currentstorms()
        #self.get_atcfname()

    def input_storm(self):
        for i, storm in enumerate(self.storms, 1):
            lat_repr, lon_repr = self._repr_lat_lon(storm.slat, storm.slon)
            print('{}. {} {} {}'.format(i, storm.get_fullname(), lat_repr, lon_repr))
        while True:
            index = input('Select a storm>')
            if index.isdigit() and 0 < int(index) <= len(self.storms):
                self.set_selected([self.storms[int(index)-1]])
                break
    
    def set_selected(self, storms):
        self.selected = storms

    def download_selected(self):
        self.downer.download(self.selected)

    def plot(self):
        EEMN(self.selected).run()

    def hack_mode(self):
        pass


if __name__ == '__main__':
    UserInterface().interface()
