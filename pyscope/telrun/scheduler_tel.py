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
from pyscope.telrun.schedule_utils import format_exptime, plot_schedule_gantt, plot_schedule_sky

from .. import utils
from ..observatory import Observatory, reconfig
from . import sch, schedtab

"""
This module contains the Scheduler Script
"""

# TODO - Add click options for the remaining parameters
@click.command(
    epilog="""Check out the documentation at
               https://pyscope.readthedocs.io/en/latest/
               for more information."""
)
@click.version_option
def scheduler_cli(    
    catalog_path=None,
    block_read_method=None, 
    date=None,
    observatory=None,
    max_altitude=-12,
    elevation=30,
    airmass=3,
    moon_separation=30,
    scheduler_type=("", ""),
    gap_time=60,
    resolution=5,
    name_format="{code}_{target}_{filter}_{exposure}s_{start_time}",
    filename=None,
    telrun=False,
    plot=None,
    quiet=False,
    verbose=0,
    reconfig_file=None,
):
    """
    Scheduler CLI

    This is the command line interface for the scheduler.

    Parameters
    ----------
    catalog_path : str
        Path to a '.sch' file or '.cat' file containing observation blocks. 
        TODO: Support for JSON files is not implemented yet.
        If not provided, defaults to 'schedule.cat' in '$TELHOME/schedules/'.
    block_read_method : str
        Method to read the blocks. Default is None for debugging
        # TODO: Can get this information from catalog_path extension - Remove this option
    date : str
        Date to schedule for. If not provided, the current date is used.
    observatory : ???
        Observatory to use. If not provided, the default observatory is used.
    max_altitude : float
        Maximum altitude of the sun in degrees. Default is -12 degrees.
    elevation : float
        Minimum elevation of the target in degrees. Default is 30 degrees.
    airmass : float
        Maximum airmass of the target. Default is 3.
    schedule_length : int
        Length of the schedule in days. Defaults to 1.
    moon_separation : float
        Minimum separation of the moon in degrees. Default is 30 degrees.
    scheduler_type : tuple # TODO: Unused
        Scheduler to use. Default is ("", ""). TODO: Add more options.
    gap_time : int
        Gap time between observations in seconds. Default is 60 seconds.
    resolution : int
        Time resolution of the scheduler in seconds. Default is 5 seconds.
    name_format : str
        Format string for the name of the observation. Default is "{code}_{target}_{filter}_{exposure}s_{start_time}".
    filename : str
        Filename to save the schedule to. If not provided, the schedule is not saved.
    telrun : bool
        If True, the schedule is run in telrun. Default is False.
    plot : str # TODO: to be implemented
        Plot the schedule. If "gantt", a Gantt plot is created. If "sky", a sky plot is created. Default is None.
    quiet : bool  # TODO: Unused
        If True, suppresses all output. Default is False.
    verbose : int # TODO: Unused
        Verbosity level.
    reconfig_file : str
        Path to the reconfiguration file. If not provided, the default reconfiguration file is used.
        
    Returns
    -------
    schedule_table : astropy.table.Table
        The schedule table.

    """
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    logger.info("Starting Scheduler CLI")
    # Load configuration
    #TODO Parse supplied arguments / Error Checking
    #TODO: Error handling of arguments
    # Set the schedule time [The current time when running the script]
    sched_time = astrotime.Time.now()
    sched_time.format = "mjd"
    # Create Global Constraints
    logger.info("Defining global constraints")
    global_constraints = [
        astroplan.AtNightConstraint(max_solar_altitude=max_altitude * u.deg),
        astroplan.AltitudeConstraint(min=elevation * u.deg),
        astroplan.AirmassConstraint(max=airmass, boolean_constraint=False),
        astroplan.MoonSeparationConstraint(min=moon_separation * u.deg),
    ]
    # Create Observatory Object
    observer = _load_observatory_config(observatory)
    # Get the observatory location
    location = observer.location #TODO: Unused currently
    # Load hardware configuration
    # TODO: Load hardware configuration from the config file
    instrument_reconfig_times = {"filter": {"default": 1 * u.second}} 
    # Create Transitioner Object
    transitioner = astroplan.Transitioner(
        observer.slew_rate, instrument_reconfig_times=instrument_reconfig_times
    )
    # Initialize the priority scheduler
    # TODO: Add more options to select the scheduler :: scheduler_type
    scheduler = astroplan.PriorityScheduler(
        constraints=global_constraints,
        observer=observer,
        transitioner=transitioner,
        gap_time=gap_time * u.second,
        time_resolution=resolution * u.second,
    )
    # Get the schedule time range
    start_time, end_time = _get_schedule_times(
        observer, date=date, schedule_length=1 * u.day
    )
    schedule = astroplan.Schedule(start_time, end_time)
    # Create blocks [Hardcoded for testing]
    # TODO - Choose import method [JSON file, .sch file]
    blocks = []
    if block_read_method == "json":
        # TODO - Create Blocks from JSON file
        blocks.append(_parse_json(catalog_path))
    elif block_read_method == "sch":
        # TODO - Create Blocks from .sch file
        blocks.append(_parse_sch(catalog_path))
    else:
        logging.error(
            "Block read method must be 'json' or 'sch'. Entering debug mode"
        )
        # Default Blocks
        deneb = astroplan.FixedTarget.from_name('Deneb')
        blocks = [astroplan.ObservingBlock(deneb, 20*u.minute, 0)]
        m13 = astroplan.FixedTarget.from_name('M13')
        blocks.append(astroplan.ObservingBlock(m13, 20*u.minute, 0))
    # Call the schedule with the observing blocks and schedule to schedule the blocks
    scheduler(blocks, schedule)
    # Create the schedule table
    logger.info("Creating schedule table")
    schedule_table = schedule.to_table()
    # Construct the path for the schedule file and output
    logger.debug("Creating schedule path")
    if filename is None or telrun:
        first_time = schedule[0]["start_time"].strftime("%Y-%m-%dT%H-%M-%S")
        filename = "telrun_" + first_time + ".ecsv"
    if not telrun:
        write_fname = filename
    else:
        path = os.environ.get("TELRUN_EXECUTE")
        if path is None:
            path = os.environ.get("TELHOME")
        if path is None:
            path = os.getcwd() + "/schedules/execute/"
        else:
            path += "/schedules/execute/"
        if not os.path.exists(path):
            os.makedirs(path)
        logger.info("Creating directory %s" % path)
        write_fname = path + filename
    # Write the schedule to a file
    logger.debug("Writing schedule to file")
    schedule_table.write(write_fname, overwrite=True, format="ascii.ecsv")

def _load_observatory_config(observatory=None):
    """
    Load the observatory configuration from the config file.

    Parameters
    ----------
    observatory : str or Observatory or astroplan.Observer
        Observatory to use. 
        If a string, it is assumed to be the path to the configuration file.
        This configuration file is used to create an astroplan.Observer object.
        If an Observatory object, it is converted to an astroplan.Observer object.
        If an astroplan.Observer object, it is used as is.

    Returns
    -------
    observatory : astroplan.Observer
        The astroplan observer object. 

    """

    logger = logging.getLogger(__name__)
    logger.info("Loading the observatory config")
    if type(observatory) is str:
        obs_cfg = configparser.ConfigParser()
        obs_cfg.read(observatory)
        slew_rate = obs_cfg["scheduling"].getfloat("slew_rate") * u.deg / u.second
        instrument_reconfig_times = json.loads(
            obs_cfg["scheduling"].get("instrument_reconfig_times")
        )
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
        slew_rate = observatory.slew_rate * u.deg / u.second
        instrument_reconfig_times = observatory.instrument_reconfig_times
        observatory = astroplan.Observer(
            location=coord.EarthLocation(lon=obs_lon, lat=obs_lat)
        )
    elif type(observatory) is astroplan.Observer:
        obs_lon = observatory.location.lon
        obs_lat = observatory.location.lat
        slew_rate = observatory.slew_rate * u.deg / u.second
        instrument_reconfig_times = observatory.instrument_reconfig_times
    else:
        logger.error(
            "Observatory must be, a string, Observatory object, or astroplan.Observer object."
        )
        return

    return observatory

def _get_schedule_times(observer,date=None,schedule_length=1):
    """
    Get the start and end times for the schedule.

    Parameters
    ----------
    observer : astroplan.Observer
        The observer object.
    date : str
        The date to schedule for. If not provided, the current date is used.
        The date should be in the format YYYY-MM-DD.
    schedule_length : int
        The length of the schedule in days. Default is 1 day.

    Returns
    -------
    start_time : astropy.time.Time
        Start time of the schedule.
    end_time : astropy.time.Time
        End time of the schedule.

    """

    logger = logging.getLogger(__name__)
    # Get Timezone of the observer
    tz = timezonefinder.TimezoneFinder().timezone_at(lng=observer.location.lon.deg, lat=observer.location.lat.deg)
    tz = zoneinfo.ZoneInfo(tz)
    logger.debug(f"tz = {tz}")
    # Check the supplied date
    if date is None:
        logger.debug("Using current date at observatory location")
        date = datetime.datetime.now()
    else:
        date = datetime.datetime.strptime(date, "%Y-%m-%d")
    date = datetime.datetime(date.year, date.month, date.day, 12, 0, 0, tzinfo=tz)
    start_time = astrotime.Time(
        datetime.datetime(date.year, date.month, date.day, 12, 0, 0, tzinfo=tz),
        format="datetime",
    )
    end_time = start_time + schedule_length + u.day
    logger.info("Schedule time range: %s to %s (UTC)" % (start_time.iso, end_time.iso))
    return start_time, end_time


def _parse_json(catalog_path):
    """
    Parse the JSON file and create the blocks. 
    TODO: Not implemented yet

    Parameters
    ----------
    catalog_path : str
        Path to the JSON file.

    Returns
    -------
    blocks : list [astroplan.ObservingBlock]
        List of blocks created from the JSON file.

    Example
    -------
        deneb = astroplan.FixedTarget.from_name('Deneb')
        blocks = [astroplan.ObservingBlock(deneb, 20*u.minute, 0)]
        m13 = astroplan.FixedTarget.from_name('M13')
        blocks.append(astroplan.ObservingBlock(m13, 20*u.minute, 0))

    """

    logger = logging.getLogger(__name__)
    logger.info("Parsing JSON file")
    with open(catalog_path, "r") as f:
        data = json.load(f)
        blocks = []
        for block in data["blocks"]:
            block["target"] = astroplan.FixedTarget.from_name(block["target"])
            blocks.append(block)
    return blocks

def _parse_sch(catalog_path):
    """
    Parse the .sch file and create the blocks. 
    TODO: Not implemented yet

    Parameters
    ----------
    catalog_path : str
        Path to the .sch file.

    Returns
    -------
    blocks : list [astroplan.ObservingBlock]
        List of blocks created from the .sch file.

    Example
    -------
        deneb = astroplan.FixedTarget.from_name('Deneb')
        blocks = [astroplan.ObservingBlock(deneb, 20*u.minute, 0)]
        m13 = astroplan.FixedTarget.from_name('M13')
        blocks.append(astroplan.ObservingBlock(m13, 20*u.minute, 0))

    """
    
    logger = logging.getLogger(__name__)
    logger.info("Parsing .sch file")
    with open(catalog_path, "r") as f:
        data = f.read()
        blocks = []
        for block in data.split("\n"):
            block = block.split(",")
            block[0] = astroplan.FixedTarget.from_name(block[0])
            blocks.append(block)
    return blocks