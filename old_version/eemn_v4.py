from mpkit import plotplus

import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.patch import geos_to_path
from matplotlib.collections import PathCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch, Circle
from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
from shapely.geometry import LineString, MultiLineString, box, Point
from shapely.ops import cascaded_union

import trackdown

#EEMN Show tracks and strike probabilities of tropical cyclones based on ECMWF ensemble
#---Configuration
RES = 0.1 #spatial resolution of probabilities
RADIUS = 1 #size of strike swath, unit in lat/lon
THRESHOLD = 5 #zoom in to only show area where probabilities are greater than THRESHOLD
SHOW_PROB = 20 #only show cities where probabilities are greater than SHOW_PROB
NUM_OF_MEMS = 50 #num of ensemble members
DAYS_AHEAD = None #only draw tracks during next few days
NO_TRACK = False #no ensemble track display
NO_SUBPLOT = False #no subplot
REALTIME_POS = True #show realtime position (using JTWC fix)
TOP_CITIES = 9 #Maxmimum num of cities to be shown
ADVANCE_LINES = False #Advance line based on tracks

#---What's new
#2.0# Upgraded from EEMN v1
#2.1# Strike probabilities added
#2.2# Intensity column added
#2.3# Probabilities of cities added
#2.4# Support for multiple storms / Improved user interface
#3.0# Shapely implement / Coastline subplot / Ensemble control / Days limit
#3.1# Wind radii / Colorbar relocated / Improved intensity column
#3.2# New colormap / Add NO_TRACK flag / Realtime fix display
#4.0# Class oriented refactor

necoast = r'C:\Users\hao\AppData\Roaming\SPB_16.6\.local\share\cartopy\shapefiles\natural_earth\physical\10m_coastline'

def loadall(time, bnum, codes, days=None):
    return [trackdown.select(trackdown.load(time, code), bnum, days=days) for code in codes]

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

def npmin(x):
    try:
        return np.amin(x)
    except:
        return 0

def getcolor(i):
    if i == 0:
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

def calculate_point(lat, lon, distance, begin=0, end=90, res=1):
    R = 6378.1
    distance *= 1.852 #nmile to km
    if distance == 0:
        return [(lon, lat)]
    bearing = np.radians(np.arange(begin, end+res, res))
    lat = np.radians(lat)
    lon = np.radians(lon)
    lat2 = np.arcsin(np.sin(lat) * np.cos(distance/R) + np.cos(lat) * 
        np.sin(distance/R) * np.cos(bearing))
    lon2 = lon + np.arctan2(np.sin(bearing) * np.sin(distance/R) *
        np.cos(lat), np.cos(distance/R) - np.sin(lat) * np.sin(lat2))
    lat2 = np.degrees(lat2)
    lon2 = np.degrees(lon2)
    return list(zip(lon2, lat2))

def wind_radii_path(lat, lon, radii):
    ne, se, sw, nw = radii
    nev = calculate_point(lat, lon, ne)
    sev = calculate_point(lat, lon, se, begin=90, end=180)
    swv = calculate_point(lat, lon, sw, begin=180, end=270)
    nwv = calculate_point(lat, lon, nw, begin=270, end=360)
    ver = nev + sev + swv + nwv
    ver.append(nev[0])
    cod = [Path.LINETO for i in ver]
    cod[0] = Path.MOVETO
    cod[-1] = Path.CLOSEPOLY
    return Path(ver, cod)

class EEMN:
    code_eps = ['ep%02d' % i for i in range(1, 26)]
    code_ens = ['en%02d' % i for i in range(1, 26)]
    code_sps = ['emx', 'eemn', 'ec00']

    def run(self):
        self.interface()
        self.init()
        self.analyze()
        self.set_map()
        self.plot_probs()
        self.plot_tracks()
        if ADVANCE_LINES:
            self.plot_advance_lines()
        self.plot_legend()
        self.plot_infos()
        if not self.MULTIFLAG:
            self.plot_city_probs()
            self.plot_columns()
            if REALTIME_POS:
                self.plot_realtime_position()
            if not NO_SUBPLOT:
                self.plot_subplot()
        self.save()

    def interface(self):
        ##INPUT
        self.time = input('Time>')
        raw_tracks = trackdown.load(self.time, 'EMX')
        opts = np.unique(np.core.defchararray.add(raw_tracks['basin'], raw_tracks['num']))
        optdict = dict()
        for i, o in enumerate(opts, 1):
            print('%d.%s' % (i, o), end=' ')
            optdict.update({i:o})
        if len(opts) == 1:
            choices = [opts[0]]
            print()
        else:
            rawstr = input('\nSelect>')
            if not rawstr:
                exit(0)
            choices = [optdict[int(s)] for s in rawstr.split()]
        names = list()
        for ch in choices:
            name = input('Input name for %s>' % ch)
            if name == '':
                name = 'INVEST'
            names.append(name.upper())
        self.bnums = choices
        self.names = names
        self.MULTIFLAG = len(self.bnums) > 1

    def init(self):
        ###LOAD TRACKS / CALCULATE RANGE
        print('...CALCULATING FIGURE RANGE...')
        self.data = list()
        geos = list()
        for bnum in self.bnums:
            eps = loadall(self.time, bnum, self.code_eps, days=DAYS_AHEAD)
            ens = loadall(self.time, bnum, self.code_ens, days=DAYS_AHEAD)
            sps = loadall(self.time, bnum, self.code_sps, days=DAYS_AHEAD)
            georange = self.calc_georange(eps, ens, sps)
            self.data.append((eps, ens, sps))
            geos.append(georange)
        georange = self.merge_georanges(geos)
        latmin, latmax, lonmin, lonmax = georange
        self.georange = roundit(geoscale(latmin, latmax, lonmin, lonmax))

    def analyze(self):
        ###ANALYZE TRACKS
        print('...ANALYZING STORM TRACKS...')
        allgrids = list()
        x, y, self.grid = get_grids(self.georange)
        self.xshape = x.shape[0]
        self.yshape = y.shape[0]
        xx, yy = np.meshgrid(x, y)
        self.xy = np.dstack((xx, yy)).reshape((self.xshape * self.yshape, 2))
        for i, bnum in enumerate(self.bnums):
            eps, ens, sps = self.data[i]
            storm_grid = self.add_storm(eps + ens, bnum)
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
        for eps, ens, sps in self.data:
            self.intens = list()
            # Deterministic
            emx = sps[0]
            self.p.plot(emx['lon'], emx['lat'], marker='o', markersize=2, mec='none',
                        linestyle='-', lw=0.5, color='#8877CC')
            self.intens.append(('DET', npmin(emx['pres'])))
            # Mean
            eemn = sps[1]
            self.p.plot(eemn['lon'], eemn['lat'], marker='o', markersize=2, mec='none',
                        linestyle='-', lw=0.5, color='#99DD22')
            self.intens.append(('MEAN', npmin(eemn['pres'])))
            # Control
            ec00 = sps[2]
            self.p.plot(ec00['lon'], ec00['lat'], marker='o', markersize=2, mec='none',
                        linestyle='-', lw=0.5, color='#AAAAAA')
            self.intens.append(('CTRL', npmin(ec00['pres'])))
            # Members
            for i, e in enumerate(eps, 1):
                if not NO_TRACK:
                    self.p.plot(e['lon'], e['lat'], marker=None, linestyle='-', lw=0.3,
                                color='#CCCCCC')
                self.intens.append(('EP%02d' % i, npmin(e['pres'])))
            for i, e in enumerate(ens, 1):
                if not NO_TRACK:
                    self.p.plot(e['lon'], e['lat'], marker=None, linestyle='-', lw=0.3,
                                color='#CCCCCC')
                self.intens.append(('EN%02d' % i, npmin(e['pres'])))
    
    def plot_advance_lines(self):
        return
        ##PLOT: PLOT ADVANCE LINE
        print('...PLOTING ADVANCE LINES...')
        for i, storm in enumerate(self.data):
            tracks = np.concatenate(storm[0] + storm[1])
            for fcsthour in [24, 48, 72]:
                points = tracks[tracks['fcsthour']==fcsthour]
                lons = points['lon']
                lats = points['lat']
                boundaries = []
                if len(lons) == 0:
                    continue
                for lon, lat in zip(lons, lats):
                    self.p.plot(lon, lat, marker='x', ms=3, mec='k')
                    boundaries.append(Point(lon, lat).buffer(RADIUS/2))
                path = geos_to_path(cascaded_union(boundaries))[0]
                self.p.ax.add_patch(PathPatch(path, ec='k', fc='none'))
                boolarr = path.contains_points(self.xy).reshape((self.yshape, self.xshape)).astype(np.uint8)
                np.save('data', boolarr)
                np.save('prob', self.prob_grid)
                from scipy.ndimage.morphology import binary_fill_holes
                from skimage.morphology import binary_erosion, binary_dilation, skeletonize_3d, disk
                #binary_fill_holes(boolarr, output=boolarr)
                core = disk(5)
                boolarr = binary_dilation(binary_erosion(boolarr, core), core)
                boolarr = skeletonize_3d(boolarr)
                ys, xs = np.where(boolarr==255)
                latmin, latmax, lonmin, lonmax = self.georange
                for y, x in zip(ys, xs):
                    self.p.ax.plot(lonmin+x*RES, latmin+y*RES, marker='o', ms=1, color='k')

    def plot_legend(self):
        ###PLOT: PLOT LEGEND
        h_e = Line2D([], [], color='#CCCCCC', lw=0.3, marker=None, label='Ensemble Cluster')
        h_x = Line2D([], [], color='#8877CC', lw=0.5, marker='o', ms=2, mec='none', label='Deterministic')
        h_n = Line2D([], [], color='#99DD22', lw=0.5, marker='o', ms=2, mec='none', label='Ensemble Mean')
        h_c = Line2D([], [], color='#AAAAAA', lw=0.5, marker='o', ms=2, mec='none', label='Ensemble Control')
        handles = [h_e, h_x, h_n, h_c]
        self.p.legend(handles=handles, loc='lower right', framealpha=0.8)

    def plot_infos(self):
        ###PLOT: PLOT INFORMATION
        namestr = '%s (%s)' % (self.bnums[0], self.names[0])
        for b, n in zip(self.bnums[1:], self.names[1:]):
            namestr = namestr + ' & %s (%s)' % (b, n)
        hourstr = '(Within {:d} hours)'.format(int(DAYS_AHEAD * 24)) if DAYS_AHEAD else ''
        self.p.title('Strike Probabilites* of %s Based on ECMWF Ensemble %s' % (namestr, hourstr), nasdaq=False)
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
            self.p.ax.text(1.01, i * 0.02, s, color=c, fontsize=5, family='Segoe UI Emoji',
                           transform=self.p.ax.transAxes)

    def plot_realtime_position(self):
        import pybst
        basin_codes = {
            'WP': 'W',
            'EP': 'E',
            'AL': 'L',
            'SI': 'S',
            'SP': 'P'
        }
        bnum = self.bnums[0]
        current_storms = pybst.currentstorm()
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
        print(cgeorange)
        cgeorange = geoscale(*cgeorange, scale=0.43)
        self.set_subplot_map(highlights, cgeorange)
        self.plot_subplot_track_and_wind_radii(highlights)

    def save(self):
        dpi = 250 if self.MULTIFLAG else 200
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

    def plot_subplot_track_and_wind_radii(self, highlights):
        ###SUBPLOT: PLOT DETERMINISTIC TRACK
        emx = self.data[0][-1][0]
        xlon, xlat = emx['lon'], emx['lat']
        self.p.plot(xlon, xlat, marker='o', markersize=2, mec='none', linestyle='-',
                    lw=0.5, color='#CCCCCC')
        ###SUBPLOT: PLOT WIND RADII
        highcoasts = MultiLineString([j for i in highlights[2:] for j in i])
        xraw = trackdown.select(trackdown.load(self.time, 'emx'), self.bnums[0],
                                days=DAYS_AHEAD, raw=True)
        radii_record = ['r1', 'r2', 'r3', 'r4']
        for i in range(len(xlon)-1):
            lineseg = LineString([(xlon[i], xlat[i]), (xlon[i+1], xlat[i+1])])
            if True:
            #if lineseg.intersects(highcoasts):
                fh = emx['fcsthour'][i]
                xr = xraw[np.logical_and(xraw['fcsthour']==fh, xraw['rad']==34)]
                if len(xr) > 0 and xr['wind'].item() >= 34:
                    radii = xr[radii_record].item()
                    radiipath = wind_radii_path(xlat[i], xlon[i], radii)
                    self.p.ax.add_patch(PathPatch(radiipath, fc='none', lw=0.5, ec='#CCCCCC'))
                    self.p.plot(xlon[i], xlat[i], color='k', marker='o', mec='none', ms=2)
                xr = xraw[np.logical_and(xraw['fcsthour']==fh, xraw['rad']==50)]
                if len(xr) > 0 :
                    radii = xr[radii_record].item()
                    radiipath = wind_radii_path(xlat[i], xlon[i], radii)
                    self.p.ax.add_patch(PathPatch(radiipath, fc='none', lw=0.5, ec='#FF9999'))
                xr = xraw[np.logical_and(xraw['fcsthour']==fh, xraw['rad']==64)]
                if len(xr) > 0 :
                    radii = xr[radii_record].item()
                    radiipath = wind_radii_path(xlat[i], xlon[i], radii)
                    self.p.ax.add_patch(PathPatch(radiipath, fc='none', lw=0.5, ec='brown', hatch=':'))
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

    def calc_georange(self, *args):
        pad = 2
        ranges = list()
        for i in args:
            for e in i:
                if e['lat'].shape[0] != 0:
                    ranges.append([np.amin(e['lat']), np.amax(e['lat']), np.amin(e['lon']), np.amax(e['lon'])])
        ranges = np.ma.masked_invalid(np.array(ranges))
        latmin = np.amin(ranges[:,0]) - pad
        latmax = np.amax(ranges[:,1]) + pad
        lonmin = np.amin(ranges[:,2]) - pad
        lonmax = np.amax(ranges[:,3]) + pad
        return latmin, latmax, lonmin, lonmax
    
    def merge_georanges(self, geos):
        if len(geos) == 1:
            return geos[0]
        latmins = [g[0] for g in geos]
        latmaxs = [g[1] for g in geos]
        lonmins = [g[2] for g in geos]
        lonmaxs = [g[3] for g in geos]
        return min(latmins), max(latmaxs), min(lonmins), max(lonmaxs)

    def add_storm(self, tracks, bnum):
        count = self.grid.copy()
        length = len(tracks)
        for i, t in enumerate(tracks, 1):
            linestring = list(zip(t['lon'], t['lat']))
            if len(linestring) < 2:
                continue
            path = geos_to_path(LineString(linestring).buffer(RADIUS))[0]
            boolarr = path.contains_points(self.xy).reshape((self.yshape, self.xshape)).astype(np.uint8)
            count += boolarr
            print('\r{:s}: [{: <10s}] {:.0%}'.format(bnum, '#'*(i*10//length), i/length), end='')
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

if __name__ == '__main__':
    EEMN().run()
