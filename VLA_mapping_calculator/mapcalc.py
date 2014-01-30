# https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/mosaicking
#
from astropy import units as u
from numpy import pi,sqrt,log,sin,cos,ceil
import numpy as np
from astropy import coordinates

def united(qty, unit):
    if isinstance(qty,u.Quantity):
        return qty.to(unit)
    else:
        return qty*u.Unit(unit)

def uvalue(qty, unit):
    return united(qty, unit).value


class MappingObservation(object):

    def __init__(self, xlen, ylen, freq, 
                 pointing_center=None,
                 integration_time_per_beam=0*u.s,
                 overhead_per_beam=0*u.s,
                 turnaround_overhead=30*u.s):

        self.xlen = united(xlen,u.arcsec)
        self.ylen = united(ylen,u.arcsec)
        self.freq = united(freq,u.GHz)
        if pointing_center is None:
            self.poining_center = coordinates.Galactic(0,0,unit=(u.deg,u.deg))
        else:
            self.pointing_center = pointing_center

        self.integration_time_per_beam = integration_time_per_beam
        self.overhead_per_beam = overhead_per_beam

    def summary(self):
        print "Scan length:      {:10}   {:10}".format(self.scan_length,self.scan_length.to(u.minute))
        print "Number of scans:  {:10}        ".format(self.nscans)
        print "Total time:       {:10}   {:10}   {:10.1f}".format(self.total_time, self.total_time.to(u.minute), self.total_time.to(u.hour))
        print "Pixels per beam:  {:10}".format(self.pixels_per_beam)
        print "Time per beam:    {:10}".format(self.time_per_beam)

    @property
    def hex_pointings(self):
        spacing = self.beam/sqrt(2)
        #angles = range(0,360,60)
        npt_vertical   = (self.ylen / (spacing*sin(60/180.*pi))).decompose().value
        npt_horizontal = (self.xlen / spacing).decompose().value
        pointing_centers = [
                self.pointing_center.__class__(
                    self.pointing_center.lonangle + (spacing * ii) + (spacing * (jj % 2) * cos(60/180.*pi)),
                    self.pointing_center.latangle + (spacing * sin(60/180.*pi) * jj))
                for ii in np.arange(-ceil(npt_horizontal/2.),ceil(npt_horizontal/2.)+1)
                for jj in np.arange(-ceil(npt_vertical/2.),ceil(npt_vertical/2.)+1)]

        return pointing_centers

    def hex_pointings_regions(self, outfile=None):

        beam = self.beam

        regions = ["circle({0:f},{1:f},{2:f}')".format(c.icrs.ra.deg,
                                                      c.icrs.dec.deg,
                                                      (beam / 2).to(u.arcmin).value)
                    for c in self.hex_pointings]

        if outfile:
            with open(outfile,'w') as outf:
                outf.write('fk5\n')
                outf.writelines("\n".join(regions))

        return regions

    @property
    def beam(self):
        return (45 * u.arcmin * (1*u.GHz/self.freq)).to(u.arcmin)

    @property
    def scan_length(self):
        if self.scan_direction == 'x':
            return ((self.xlen / self.pixelsize) * self.samples_per_pixel * self.dumpeach).to(u.s)
        if self.scan_direction == 'y':
            return ((self.ylen / self.pixelsize) * self.samples_per_pixel * self.dumpeach).to(u.s)

    @property
    def nscans(self):
        if self.scan_direction == 'x':
            return self.ylen / self.pixelsize * self.samples_per_pixel
        if self.scan_direction == 'y':
            return self.xlen / self.pixelsize * self.samples_per_pixel

    @property
    def total_time(self):
        return (self.n_beams * (self.integration_time_per_beam+self.overhead_per_beam))

    @property
    def total_time_onsource(self):
        return (self.n_beams * self.integration_time_per_beam)

    @property
    def beam_area(self):
        return pi * (self.beam/sqrt(8*log(2)))**2

    @property
    def n_beams(self):
        return (self.xlen * self.ylen / self.beam_area).decompose().value

    @property
    def pixels_per_beam(self):
        return (self.beam_area / self.pixelsize**2).decompose()

def VLA_H2CO_CMZ():
    M1 = MappingObservation(90*u.arcmin,30*u.arcmin,
                           4.829*u.GHz,
                           pointing_center=coordinates.Galactic(0*u.deg,0*u.deg),
                           integration_time_per_beam=17.5*u.min,
                           overhead_per_beam=4.5*u.min
                           )
    print "Total time: ",M1.total_time.to(u.hour)
    print "Total time ON: ",M1.total_time_onsource.to(u.hour)
    ptgs = M1.hex_pointings
    print "N(beams): ",M1.n_beams
    print "N(pointings): ",len(ptgs)
    print

    M2 = MappingObservation(53*u.arcmin,11*u.arcmin,
                            14.488*u.GHz,
                            pointing_center=coordinates.Galactic(0.293*u.deg,-0.058*u.deg),
                            integration_time_per_beam=15*u.min,
                            overhead_per_beam=13.35*u.min)
    print "Total time: ",M2.total_time.to(u.hour)
    print "Total time ON: ",M2.total_time_onsource.to(u.hour)
    ptgs = M2.hex_pointings
    print "N(beams): ",M2.n_beams
    print "N(pointings): ",len(ptgs)
    print

    M1c = MappingObservation(53*u.arcmin,11*u.arcmin,
                            4.829*u.GHz,
                            pointing_center=coordinates.Galactic(0.293*u.deg,-0.058*u.deg),
                            integration_time_per_beam=1.5*u.hour,
                            overhead_per_beam=0.25*1.5*u.hour)
    print "Total time: ",M1c.total_time.to(u.hour)
    print "Total time ON: ",M1c.total_time_onsource.to(u.hour)
    ptgs = M1c.hex_pointings
    print "N(beams): ",M1c.n_beams
    print "N(pointings): ",len(ptgs)
    print

    M1b = MappingObservation(90*u.arcmin,30*u.arcmin,
                             4.829*u.GHz,
                             pointing_center=coordinates.Galactic(0*u.deg,0*u.deg),
                             integration_time_per_beam=20*u.min,
                             overhead_per_beam=5*u.min
                             )
    print "Total time: ",M1b.total_time.to(u.hour)
    print "Total time ON: ",M1b.total_time_onsource.to(u.hour)
    ptgs = M1b.hex_pointings
    print "N(beams): ",M1b.n_beams
    print "N(pointings): ",len(ptgs)

    return M1,M2

M1,M2 = VLA_H2CO_CMZ()
M1.hex_pointings_regions(outfile='/Users/adam/work/gc/vla_largearea_hexpointings.reg')
