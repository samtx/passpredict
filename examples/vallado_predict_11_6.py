import datetime

from rich.console import Console
from rich.table import Table
from rich.align import Align

from passpredict import Location, Observer, SatellitePredictor, TLE
from passpredict.tle import jd_to_epoch_string
from passpredict._time import datetime2mjd


def vallado_predict_11_6():
    """
    Using orbital elements for MIR space station, p. 106 Vallado

    Note that this example propagates the satellite using SGP4.
    The example in the book uses pkepler. The results will be slightly different.
    """
    epoch_jd = 2450540.4
    tle = TLE(
        16609,
        (
            # "1 16609U 86017A   93352.53502934  .00007889  00000-0  10529-3 0   34",
            f"1 16609U 86017A   {jd_to_epoch_string(epoch_jd):14s}  .00007889  00000-0  10529-3 0   34",
            "2 16609  51.6190  13.3340 0005770 102.5680 257.5950 15.5911407044786"
        )
    )
    satellite = SatellitePredictor.from_tle(tle)
    location = Location('MIT', 42.38, -71.13, 24)
    observer = Observer(location, satellite)

    table = Table(title='Prediction Values for Example 11-6')
    table.add_column('Visibility')
    table.add_column('Range km', justify="right")
    table.add_column('\u03B2\u00B0', justify="right")
    table.add_column('el\u00B0', justify="right")
    table.add_column('Date', justify="center")
    table.add_column('UTC (h/min/s)', justify="center")

    time_step = datetime.timedelta(minutes=2)
    d = datetime.datetime(1997, 4, 1, 23, 28, 0, tzinfo=datetime.timezone.utc)
    for i in range(5):
        row = _get_row_data(observer, d)
        end_section = True if i == 4 else False
        table.add_row(*row, end_section=end_section)
        d = d + time_step

    d = datetime.datetime(1997, 4, 2, 1, 6, 0, tzinfo=datetime.timezone.utc)
    for i in range(5):
        row = _get_row_data(observer, d)
        end_section = True if i == 4 else False
        table.add_row(*row, end_section=end_section)
        d = d + time_step
        # print(pt)
    console = Console()
    console.print(table)


def _get_row_data(observer: Observer, d: datetime.datetime):
    #  aos_at_dg=0, tolerance_s=0.75
    pt = observer.point(d, visibility=True, aos_at_dg=0)
    vis = pt.type.capitalize() if pt.type else 'None'
    row = [
        vis,
        f"{pt.range:.4f}",
        f"{pt.azimuth:.4f}",
        f"{pt.elevation:.4f}",
        _format_date(pt.dt),
        _format_time(pt.dt)
    ]
    return row


def _format_date(d: datetime.datetime):
    s = f"{d.strftime('%b')} {d.day}, {d.strftime('%y')}"
    return s


def _format_time(d: datetime.datetime):
    s = f"{d.hour:02d}:{d.minute:02d}:{d.second + d.microsecond/1e6:.2f}"
    return s

if __name__ == "__main__":
    vallado_predict_11_6()
