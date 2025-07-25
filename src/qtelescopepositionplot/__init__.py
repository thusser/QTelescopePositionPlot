from typing import cast

from PyQt5.QtGui import QPalette
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.units as u
import numpy as np
from matplotlib.markers import MarkerStyle
from matplotlib.projections.polar import PolarAxes


class QTelescopePositionPlot(FigureCanvasQTAgg):
    def __init__(self, location: EarthLocation, rticks: tuple[float] | None = None):
        # create plot
        self.fig, self.ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(6, 6), dpi=60)

        # init canvas
        FigureCanvasQTAgg.__init__(self, self.fig)
        self.setContentsMargins(0, 0, 0, 0)

        # store stuff
        self._current_position = AltAz(az=0.0 * u.deg, alt=90.0 * u.deg)
        self._target_position = SkyCoord(az=0.0 * u.deg, alt=90.0 * u.deg, frame="altaz")
        self._telescope_location = location

        # ticks along altitude axis
        self._rticks = rticks if rticks is not None else (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)
        self._rtick_labels = tuple(
            [
                f"{r:.0f}°" if i in [0, len(self._rticks) // 2, len(self._rticks) - 1] else ""
                for i, r in enumerate(self._rticks)
            ]
        )

        # plot!
        self._update_plot()

    def set_current_position(self, coords: AltAz) -> None:
        self._current_position = coords
        self._update_plot()

    def set_target_position(self, coords: SkyCoord) -> None:
        self._target_position = coords
        self._update_plot()

    def _update_plot(self) -> None:
        # get some values and clear plot
        self.ax.clear()
        now = Time.now()
        pal = QPalette()

        # draw current position
        # needs az is in radians, and zenith distance in degrees
        self.ax.scatter(
            np.radians(self._current_position.az.degree),
            90 - self._current_position.alt.degree,
            label="Telescope",
            marker=MarkerStyle("o"),
            s=150,
            facecolors="none",
            edgecolors=pal.highlight().color().name(),
        )

        # is the current target position AltAz or RaDec?
        if hasattr(self._target_position, "alt"):
            target = self._target_position
        else:
            target = self._target_position.transform_to(AltAz(obstime=Time.now(), location=self._telescope_location))

        # plot current target position
        self.ax.scatter(
            np.pi / 180 * target.az.degree,
            90 - target.alt.degree,
            label="Target",
            color=pal.light().color().name(),
            marker=MarkerStyle("+"),
            s=200,
        )

        # are target position given in RaDec? then plot projection
        if hasattr(self._target_position, "ra"):
            # calculate positions for 5 hours into the future
            times = now + np.linspace(0, 5, 100) * u.hour
            target = self._target_position.transform_to(AltAz(obstime=times, location=self._telescope_location))

            # plot it
            self.ax.plot(
                np.pi / 180 * target.az.value,
                90 - target.alt.value,
                label="Path +5h",
                color=pal.light().color().name(),
                alpha=0.4,
            )

        # plot settings
        ax = cast(PolarAxes, self.ax)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(1)
        ax.set_ylim(0, 90)  # Altitude 0° (horizon) to 90° (zenith)
        ax.set_rmax(90)
        ax.set_rticks([0, 15, 30, 45, 60, 75, 90])  # Altitude rings
        ax.set_rlabel_position(225)  # Move radial labels # type: ignore
        ax.set_rgrids(self._rticks, labels=self._rtick_labels, fontsize=12)
        ax.grid(True)
        ax.set_thetalim([0, 2 * np.pi])
        ax.set_thetagrids(
            np.rad2deg(np.linspace(0, 2 * np.pi, 9)[1:]),
            labels=("45°", "E", "135°", "S", "225°", "W", "315°", "N"),
            fontsize=12,
        )

        # update plot
        self.draw()
