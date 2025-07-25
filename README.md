# QTelescopePositionPlot

Example code:

    import sys
    from PyQt5.QtWidgets import QApplication
    from astropy.coordinates import EarthLocation, AltAz, SkyCoord
    import astropy.units as u
    
    from qtelescopepositionplot import QTelescopePositionPlot
    
    app = QApplication(sys.argv)
    plot = QTelescopePositionPlot(location=EarthLocation(lat=9 * u.deg, lon=52 * u.deg, height=200 * u.m))
    plot.set_current_position(AltAz(az=23.0 * u.deg, alt=60.0 * u.deg))
    plot.set_target_position(SkyCoord(ra=134.0 * u.deg, dec=78.0 * u.deg, frame="icrs"))
    # plot.set_target_position(SkyCoord(az=45.0 * u.deg, alt=45.0 * u.deg, frame="altaz"))
    plot.show()
    sys.exit(app.exec_())
