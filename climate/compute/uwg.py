from pathlib import Path
from climate.urban_weather_generator.uwg import UWG

def run_uwg(self, parameters=None, output=None):

    input_epw = self.file_path

    if parameters is None:
        parameters = Path(__file__).resolve().parents[1] / "common" / "uwg_parameters.uwg"
    self.uwg_parameters = parameters

    if output is None:
        output = self.file_path.with_suffix('.epw_uwg')
    self.uwg_weatherfile = output

    uwg_ = UWG(
        epwFileName=str(input_epw),
        uwgParamFileName=str(parameters),
        destinationFileName=str(output)
    )
    uwg_.run()

    print("Urban weather generator successful: {}".format(output))

    return output
