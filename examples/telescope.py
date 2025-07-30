import sys

from PyQt5.QtCore import pyqtSlot, QTimer
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QHBoxLayout,
    QVBoxLayout,
    QWidget,
    QPushButton,
    QLineEdit,
    QLabel,
    QFormLayout,
    QGroupBox,
)
from astropy.coordinates import EarthLocation, AltAz, SkyCoord, ICRS
import astropy.units as u
from astropy.time import Time

from qtelescopepositionplot import QTelescopePositionPlot


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        QMainWindow.__init__(self)

        # the plot
        self.location = EarthLocation(lat=9 * u.deg, lon=52 * u.deg, height=200 * u.m)
        self.plot = QTelescopePositionPlot(location=self.location)
        self.plot.setMinimumSize(400, 400)

        # status
        group_status = QGroupBox("Status")
        layout = QFormLayout()
        group_status.setLayout(layout)
        self.textCurrentRightAscension = QLabel()
        self.textCurrentDeclination = QLabel()
        self.textCurrentAzimuth = QLabel()
        self.textCurrentAltitude = QLabel()
        layout.addRow("Az", self.textCurrentAzimuth)
        layout.addRow("Alt", self.textCurrentAltitude)
        layout.addRow("RA", self.textCurrentRightAscension)
        layout.addRow("Dec", self.textCurrentDeclination)

        # controls
        group_controls = QGroupBox("Controls")
        layout = QFormLayout()
        group_controls.setLayout(layout)
        self.buttonStop = QPushButton("Stop")
        self.textAzimuth = QLineEdit()
        self.textAltitude = QLineEdit()
        self.buttonMove = QPushButton("Move")
        self.textRightAscension = QLineEdit()
        self.textDeclination = QLineEdit()
        self.buttonTrack = QPushButton("Track")
        layout.addRow("", self.buttonStop)
        layout.addRow("Az", self.textAzimuth)
        layout.addRow("Alt", self.textAltitude)
        layout.addRow("", self.buttonMove)
        layout.addRow("RA", self.textRightAscension)
        layout.addRow("Dec", self.textDeclination)
        layout.addRow("", self.buttonTrack)

        # status/controls
        status_controls_layout = QVBoxLayout()
        status_controls_layout.addWidget(group_status)
        status_controls_layout.addWidget(group_controls)

        # put all together
        main_layout = QHBoxLayout()
        main_layout.addWidget(self.plot)
        main_layout.addLayout(status_controls_layout)
        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)

        # coordinates
        self.target = SkyCoord(0.0 * u.deg, 90.0 * u.deg, frame="altaz")
        self.current = SkyCoord(0.0 * u.deg, 90.0 * u.deg, frame="altaz")

        # buttons
        self.buttonStop.clicked.connect(self.stop)
        self.buttonMove.clicked.connect(self.move_altaz)
        self.buttonTrack.clicked.connect(self.track_radec)

        # updater
        self.update_timer = QTimer()
        self.update_timer.setInterval(100)
        self.update_timer.timeout.connect(self.update_position)
        self.update_timer.start()

    def stop(self) -> None:
        self.target = self.current

    def move_altaz(self) -> None:
        self.target = SkyCoord(
            float(self.textAzimuth.text()) * u.deg, float(self.textAltitude.text()) * u.deg, frame="altaz"
        )
        self.plot.set_target_position(self.target)

    def track_radec(self) -> None:
        self.target = SkyCoord(
            float(self.textRightAscension.text()) * u.deg, float(self.textDeclination.text()) * u.deg, frame="icrs"
        )
        self.plot.set_target_position(self.target)

    @pyqtSlot()
    def update_position(self) -> None:
        # convert target to AltAz
        if hasattr(self.target, "ra"):
            target = self.target.transform_to(AltAz(location=self.location, obstime=Time.now()))
        else:
            target = self.target

        # print(target)

        # distance small enough?
        alt, az = self.current.alt, self.current.az
        if abs(az - target.az) < 1.0 * u.deg:
            az = target.az
        elif az < target.az:
            az += 0.5 * u.deg
        else:
            az -= 0.5 * u.deg
        if abs(alt - target.alt) < 1.0 * u.deg:
            alt = target.alt
        elif alt < target.alt:
            alt += 0.5 * u.deg
        else:
            alt -= 0.5 * u.deg
        self.current = SkyCoord(az, alt, frame="altaz", obstime=Time.now(), location=self.location)

        # to RaDec
        radec = self.current.transform_to(ICRS())
        ra, dec = radec.to_string("hmsdms", sep=":", precision=1).split(" ")

        # set
        self.textCurrentAzimuth.setText(f"{self.current.az.degree:.2f}")
        self.textCurrentAltitude.setText(f"{self.current.alt.degree:.2f}")

        self.textCurrentRightAscension.setText(ra)
        self.textCurrentDeclination.setText(dec)

        # plot
        self.plot.set_current_position(self.current)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
