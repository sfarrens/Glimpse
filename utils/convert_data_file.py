#!/usr/bin/python
# This script converts txt files to fits format for use with glimpse
import pyfits as fits
from pyfits import Column
import sys
import numpy as np

def main(argv):
    data = np.loadtxt(argv[0])

    ra   = Column(name='ra_gal',    format='D', array=data[:, 0])
    dec  = Column(name='dec_gal',   format='D', array=data[:, 1])
    z    = Column(name='z',         format='D', array=data[:, 2])
    zphot= Column(name='z_phot',    format='D', array=data[:, 3])
    zsig = Column(name='z_sig',     format='D', array=data[:, 4])
    e1   = Column(name='e1_gal',    format='D', array=data[:, 5])
    e2   = Column(name='e2_gal',    format='D', array=data[:, 6])
    f1   = Column(name='f1_gal',    format='D', array=data[:, 7])
    f2   = Column(name='f2_gal',    format='D', array=data[:, 8])   
    e1n  = Column(name='e1n_gal',   format='D', array=data[:, 9])
    e2n  = Column(name='e2n_gal',   format='D', array=data[:, 10])
    f1n  = Column(name='f1n_gal',   format='D', array=data[:, 11])
    f2n  = Column(name='f2n_gal',   format='D', array=data[:, 12])

    tbhdu = fits.BinTableHDU.from_columns([ra, dec, z, zphot, zsig, e1, e2, f1, f2, e1n, e2n, f1n, f2n])

    tbhdu.writeto(argv[1], clobber=True)

if __name__ == "__main__":
    main(sys.argv[1:])
