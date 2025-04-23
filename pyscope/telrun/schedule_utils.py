import configparser
import datetime
import json
import logging
import os
import zoneinfo

import astroplan
import click
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import timezonefinder
import tqdm
from astroplan import plots as astroplan_plots
from astropy import coordinates as coord
from astropy import table
from astropy import time as astrotime
from astropy import units as u
from astroquery import mpc
from cmcrameri import cm as ccm
from matplotlib import ticker

from .. import utils
from ..observatory import Observatory, reconfig
from . import sch, schedtab

"""
This module contains functions for plotting the schedule of observations
"""

@click.command(
    epilog="""Check out the documentation at
               https://pyscope.readthedocs.io/en/latest/
               for more information."""
)
@click.argument(
    "schedule_table",
    type=click.Path(exists=True, resolve_path=True, dir_okay=False, readable=True),
)
@click.argument(
    "observatory",
    type=click.Path(exists=True, resolve_path=True, dir_okay=False, readable=True),
)
@click.version_option()
def plot_schedule_gantt_cli(schedule_table, observatory):
    """
    Plot the schedule of observations as a Gantt chart.
    
    Parameters
    ----------
    schedule_table : ???
        TBD
    observatory : ???
        TBD
            
    """
    logger = logging.getLogger(__name__)
    if type(schedule_table) is not table.Table:
        schedule_table = table.Table.read(schedule_table, format="ascii.ecsv")

    if type(observatory) is str:
        obs_cfg = configparser.ConfigParser()
        obs_cfg.read(observatory)
        observatory = astroplan.Observer(
            location=coord.EarthLocation(
                lon=obs_cfg.get("site", "longitude"),
                lat=obs_cfg.get("site", "latitude"),
            )
        )
        obs_lon = observatory.location.lon
        obs_lat = observatory.location.lat
    elif type(observatory) is Observatory:
        obs_lon = observatory.observatory_location.lon
        obs_lat = observatory.observatory_location.lat
    elif type(observatory) is astroplan.Observer:
        obs_lon = observatory.location.lon
        obs_lat = observatory.location.lat
    else:
        logger.error(
            "Observatory must be, a string, Observatory object, or astroplan.Observer object."
        )
        return
    location = coord.EarthLocation(lon=obs_lon, lat=obs_lat)

    tz = timezonefinder.TimezoneFinder().timezone_at(lng=obs_lon.deg, lat=obs_lat.deg)
    tz = zoneinfo.ZoneInfo(tz)
    date = str(np.min(schedule_table["start_time"]-1).isot)
    date = datetime.datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f")
    t0 = astrotime.Time(
        datetime.datetime(date.year, date.month, date.day, 12, 0, 0, tzinfo=tz)
    )
    t1 = t0 + 1 * u.day

    # Only keep scheduled blocks
    schedule_table = schedule_table[schedule_table["status"] == "S"]

    obscodes = list(np.unique(schedule_table["code"]))

    fig, ax = plt.subplots(1, 1, figsize=(12, len(obscodes) * 0.75))
    mdates.set_epoch(t0.strftime("%Y-%m-%dT%H:%M:%S"))

    # Create list for all y-axis labels
    y_labels = obscodes.copy()
    y_labels.append("All")

    for i in range(len(obscodes)):
        print(f"Plotting observer {obscodes[i]}")
        plot_blocks = [
            block for block in schedule_table if block["code"] == obscodes[i]
        ]

        for block in plot_blocks:
            start_time = astrotime.Time(np.float64(block["start_time"].jd), format="jd")
            end_time = astrotime.Time(np.float64(block["end_time"].jd), format="jd")
            length_min = int((end_time - start_time).sec / 60 + 1)
            times = (
                start_time
                + np.linspace(0, length_min, length_min, endpoint=True) * u.minute
            )
            airmass = []
            for t in times:
                obj = coord.SkyCoord(
                    block["target"], obstime=t, location=location
                ).transform_to("altaz")
                airmass.append(utils.airmass(obj.alt.rad))

            ax.scatter(
                times.datetime,
                i * np.ones(len(times)),
                lw=0,
                marker="s",
                c=airmass,
                cmap=ccm.batlow,
                vmin=1,
                vmax=2.3,
            )

            scatter = ax.scatter(
                times.datetime,
                len(obscodes) * np.ones(len(times)),
                lw=0,
                marker="s",
                c=airmass,
                cmap=ccm.batlow,
                vmin=1,
                vmax=2.3,
            )

    # obscodes.append("All")

    twilight_times = [
        t0,
        astrotime.Time(observatory.sun_set_time(t0, which="next"), scale="utc"),
        astrotime.Time(
            observatory.twilight_evening_civil(t0, which="next"), scale="utc"
        ),
        astrotime.Time(
            observatory.twilight_evening_nautical(t0, which="next"), scale="utc"
        ),
        astrotime.Time(
            observatory.twilight_evening_astronomical(t0, which="next"), scale="utc"
        ),
        astrotime.Time(
            observatory.twilight_morning_astronomical(t0, which="next"), scale="utc"
        ),
        astrotime.Time(
            observatory.twilight_morning_nautical(t0, which="next"), scale="utc"
        ),
        astrotime.Time(
            observatory.twilight_morning_civil(t0, which="next"), scale="utc"
        ),
        astrotime.Time(observatory.sun_rise_time(t0, which="next"), scale="utc"),
        t1,
    ]
    opacities = [0.8, 0.6, 0.4, 0.2, 0, 0.2, 0.4, 0.6, 0.8]
    for i in range(len(twilight_times) - 1):
        ax.axvspan(
            twilight_times[i].datetime,
            twilight_times[i + 1].datetime,
            color="grey",
            alpha=opacities[i],
        )

    ax.set_xlabel(
        "Time beginning %s [UTC]"
        % (twilight_times[1] - 0.5 * u.hour).strftime("%Y-%m-%d")
    )
    """ax.set_xlim(
        [
            (twilight_times[1] - 0.5 * u.hour).datetime,
            (twilight_times[-2] + 0.5 * u.hour).datetime,
        ]
    )"""
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.xaxis.set_tick_params(rotation=45)

    ax1 = ax.twiny()
    ax1.set_xlim(ax.get_xlim())
    ax1.xaxis.set_major_locator(mdates.HourLocator(interval=2, tz=tz))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=tz))
    ax1.xaxis.set_minor_locator(mdates.HourLocator(interval=1, tz=tz))
    ax1.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax1.xaxis.set_tick_params(rotation=45)
    ax1.set_xlabel("Observatory Local Time (%s)" % tz)

    # Original y-axis labels
    # ax.set_ylabel("Observer Code")
    # ax.set_ylim([len(obscodes) - 0.5, 0.5])
    # ax.yaxis.set_major_locator(ticker.FixedLocator(np.arange(len(obscodes))))
    # ax.yaxis.set_major_formatter(ticker.FixedFormatter(obscodes))
    # ax.yaxis.set_minor_locator(ticker.NullLocator())
    # ax.yaxis.set_minor_formatter(ticker.NullFormatter())

    # Use the y_labels list for setting the y-axis labels
    ax.set_ylabel("Observer Code")
    ax.set_ylim([len(y_labels) - 0.5, -1 + 0.5])
    ax.yaxis.set_major_locator(ticker.FixedLocator(np.arange(len(y_labels))))
    ax.yaxis.set_major_formatter(ticker.FixedFormatter(y_labels))
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_ticks([1, 1.5, 2, 2.3])
    cbar.set_label("Airmass", rotation=270, labelpad=20)

    """ax.set_title(
        "Observing Schedule: %s"
        % t0.strftime("%Y-%m-%d"), fontsize=14
    )"""
    ax.grid(linestyle=":", color="black")

    fig.set_facecolor("white")
    fig.set_dpi(300)

    return fig, ax

@click.command(
    epilog="""Check out the documentation at
               https://pyscope.readthedocs.io/en/latest/
               for more information."""
)
@click.argument(
    "schedule_table",
    type=click.Path(exists=True, resolve_path=True, dir_okay=False, readable=True),
)
@click.argument(
    "observatory",
    type=click.Path(exists=True, resolve_path=True, dir_okay=False, readable=True),
)
@click.version_option()
def plot_schedule_sky_cli(schedule_table, observatory):
    """
    Plot the schedule of observations on the sky.
    
    Parameters
    ----------
    schedule_table : ???
        TBD
    observatory : ???
        TBD
        
    """
    logger = logging.getLogger(__name__)
    if type(schedule_table) is not table.Table:
        schedule_table = table.Table.read(schedule_table, format="ascii.ecsv")

    # Only keep scheduled blocks
    schedule_table = schedule_table[schedule_table["status"] == "S"]

    if type(observatory) is str:
        obs_cfg = configparser.ConfigParser()
        obs_cfg.read(observatory)
        observatory = astroplan.Observer(
            location=coord.EarthLocation(
                lon=obs_cfg.get("site", "longitude"),
                lat=obs_cfg.get("site", "latitude"),
            )
        )
        obs_lon = observatory.location.lon
        obs_lat = observatory.location.lat
    elif type(observatory) is Observatory:
        obs_lon = observatory.observatory_location.lon
        obs_lat = observatory.observatory_location.lat
    elif type(observatory) is astroplan.Observer:
        obs_lon = observatory.location.lon
        obs_lat = observatory.location.lat
    else:
        logger.error(
            "Observatory must be, a string, Observatory object, or astroplan.Observer object."
        )
        return

    # Get unique targets in the schedule
    target_times = {}

    for row in schedule_table:
        if row["name"] == "TransitionBlock" or row["name"] == "EmptyBlock":
            continue
        target_string = row["target"].to_string("hmsdms")
        target_name = row["name"]
        if target_string not in target_times:
            target_times[target_string] = {
                "name": target_name,
                "times": [row["start_time"]],
            }
        else:
            target_times[target_string]["times"].append(row["start_time"])

    # targets = [t.to_string("hmsdms") for t in schedule_table["target"]]

    fig, ax = plt.subplots(1, 1, figsize=(7, 7), subplot_kw={"projection": "polar"})
    for target, target_dict in target_times.items():
        times = target_dict["times"]
        try:
            label = target_dict["name"].strip()
        except:
            label = target
        target = coord.SkyCoord(target, unit=(u.hourangle, u.deg))
        ax = astroplan_plots.plot_sky(
            astroplan.FixedTarget(target),
            observatory,
            times,
            ax=ax,
            style_kwargs={"label": label},
        )

    handles, labels = ax.get_legend_handles_labels()
    # unique = [
    #     (h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]
    # ]
    # ax.legend(*zip(*unique), loc=(1.1, 0))

    # Add title to plot based on date
    t0 = np.min(
        schedule_table["start_time"]
    )  # -1 corrects for UTC to local time, should be cleaned up
    t0 = astrotime.Time(t0, format="mjd")
    t0.format = "iso"

    ax.set_title(
        f"Observing Schedule: Night of {t0.to_string().split(' ')[0]} UTC", fontsize=14
    )

    ax.legend(labels, loc=(1.1, 0))

    fig.set_facecolor("white")
    fig.set_dpi(300)

    return fig, ax


def format_exptime(exptime):
    return (
        f"{exptime:.0f}" if exptime.is_integer() else f"{exptime:.2g}".replace(".", "-")
    )

plot_schedule_gantt = plot_schedule_gantt_cli.callback
plot_schedule_sky = plot_schedule_sky_cli.callback