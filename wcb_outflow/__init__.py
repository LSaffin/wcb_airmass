import pathlib
import datetime

from pylagranto import trajectory

data_path = pathlib.Path(__file__).parent / "data"


class CaseStudy:

    timestep = datetime.timedelta(hours=6)
    forecast_end_lead_time = datetime.timedelta(hours=120)

    def __init__(self, name, start_time, outflow_lead_time, outflow_theta):
        self.name = name
        self.start_time = start_time
        self.outflow_lead_time = outflow_lead_time
        self.outflow_theta = outflow_theta

    @property
    def data_path(self):
        return data_path / self.name

    @property
    def datestr(self):
        return self.start_time.strftime("%Y%m%d_%H")

    @property
    def outflow_time(self):
        return self.start_time + self.outflow_lead_time

    @property
    def forecast_end_time(self):
        return self.start_time + self.forecast_end_lead_time

    def filename(self, time):
        # Files output every 12 hours
        lead_time = time - self.start_time
        file_lead_time = lead_time - (lead_time % datetime.timedelta(hours=12))
        lead_time_str = int(file_lead_time.total_seconds()) // 3600

        return "{}/prodm_op_gl-mn_{}_d{:03d}.pp".format(
            self.datestr, self.datestr, lead_time_str)

    def time_to_filename_mapping(self):
        mapping = dict()
        filename_t0 = self.filename(self.start_time)

        time = self.start_time
        while time <= self.forecast_end_time:
            mapping[time] = [filename_t0, self.filename(time)]

            time += self.timestep

        return mapping

    def filename_winds(self, time):
        return str(self.data_path / time.strftime("%Y%m%d_%H")) + "_winds.nc"

    def time_to_filename_winds_mapping(self, start=None, end=None):
        if start is None:
            start = self.start_time
        if end is None:
            end = self.outflow_time

        mapping = dict()

        time = start
        while time <= end:
            mapping[time] = self.filename_winds(time)

            time += self.timestep

        return mapping

    def filename_theta(self, time):
        # Tracer files ("b" and "c") output every 6 hours
        lead_time = time - self.start_time
        lead_time_str = int(lead_time.total_seconds()) // 3600

        filenames = [
            "{}/prodm_op_gl-mn_{}_{}{:03d}_thsfcs_5K.nc".format(
                self.name, self.datestr, letter, lead_time_str)
            for letter in ["c"]
        ]

        # Meteorological files ("d") output every 12 hours
        lead_time = lead_time - (lead_time % datetime.timedelta(hours=12))
        lead_time_str = int(lead_time.total_seconds()) // 3600

        filenames.append(
            "{}/prodm_op_gl-mn_{}_d{:03d}_thsfcs_5K.nc".format(
                self.name, self.datestr, lead_time_str)
        )

        return filenames


case_studies = dict(
    IOP3=CaseStudy(
        name="IOP3",
        start_time=datetime.datetime(2016, 9, 22, 12),
        outflow_lead_time=datetime.timedelta(hours=42),
        outflow_theta=[320, 325, 330],
    ),
    IOP5=CaseStudy(
        name="IOP5",
        start_time=datetime.datetime(2016, 9, 26, 12),
        outflow_lead_time=datetime.timedelta(hours=36),
        outflow_theta=[325, 330, 335],
    ),
    IOP6=CaseStudy(
        name="IOP6",
        start_time=datetime.datetime(2016, 9, 30, 12),
        outflow_lead_time=datetime.timedelta(hours=42),
        outflow_theta=[310, 315, 320],
    ),
    IOP7=CaseStudy(
        name="IOP7",
        start_time=datetime.datetime(2016, 10, 3, 12),
        outflow_lead_time=datetime.timedelta(hours=24),
        outflow_theta=[310, 315, 320],
    ),
)
