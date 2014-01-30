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
                 turnaround_overhead=30*u.s):

        self.xlen = united(xlen,u.arcsec)
        self.ylen = united(ylen,u.arcsec)
        self.freq = united(freq,u.GHz)
        if pointing_center is None:
            self.poining_center = coordinates.Galactic(0,0,unit=(u.deg,u.deg))
        else:
            self.pointing_center = pointing_center

        self.integration_time_per_beam = integration_time_per_beam

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
        npt_vertical   = (self.ylen / (spacing*sin(60/360.*pi))).decompose().value
        npt_horizontal = (self.xlen / spacing).decompose().value
        pointing_centers = [
                self.pointing_center.__class__(
                    self.pointing_center.lonangle + (spacing * ii),
                    self.pointing_center.latangle + (spacing * sin(60/360.*pi) * jj))
                for ii in np.arange(-ceil(npt_horizontal/2.),ceil(npt_horizontal/2.)+1)
                for jj in np.arange(-ceil(npt_vertical/2.),ceil(npt_vertical/2.)+1)]

        return pointing_centers

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
        return (self.n_beams * self.integration_time_per_pixel)

    @property
    def beam_area(self):
        return pi * (self.beam/sqrt(8*log(2)))**2

    @property
    def n_beams(self):
        return (self.xlen * self.ylen / self.beam_area).decompose().value

    @property
    def pixels_per_beam(self):
        return (self.beam_area / self.pixelsize**2).decompose()

