import pathlib
import datetime

data_path = pathlib.Path(__file__).parent / "data"


class CaseStudy:

    timestep = datetime.timedelta(hours=6)

    def __init__(self, name, start_time, outflow_lead_time):
        self.name = name
        self.start_time = start_time
        self.outflow_lead_time = outflow_lead_time

    @property
    def data_path(self):
        return data_path / self.name

    @property
    def datestr(self):
        return self.start_time.strftime("%Y%m%d_%H")

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
        outflow_time = self.start_time + self.outflow_lead_time
        while time <= outflow_time:
            mapping[time] = [filename_t0, self.filename(time)]

            time += self.timestep

        return mapping


case_studies = dict(
    IOP3=CaseStudy(
        name="IOP3",
        start_time=datetime.datetime(2016, 9, 22, 12),
        outflow_lead_time=datetime.timedelta(hours=42),
    ),
    IOP5=CaseStudy(
        name="IOP5",
        start_time=datetime.datetime(2016, 9, 26, 12),
        outflow_lead_time=datetime.timedelta(hours=36),
    ),
    IOP6=CaseStudy(
        name="IOP6",
        start_time=datetime.datetime(2016, 9, 30, 12),
        outflow_lead_time=datetime.timedelta(hours=42),
    ),
    IOP7=CaseStudy(
        name="IOP7",
        start_time=datetime.datetime(2016, 10, 3, 12),
        outflow_lead_time=datetime.timedelta(hours=24),
    ),
)
