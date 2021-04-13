from datetime import datetime
import hydrosim.delft3d as d3d
from hydrosim.functions import log, error
import json


# Define parameters
log("------------STARTING NEW RUN------------", start=True)
meteolakes_simulation_folder = "./example-sim"
output_simulation_folder = "./output-sim"
cosmo_folder = "./demo-data/Meteodata"
hydro_folder = "./demo-data/Hydrodata"
start_date = datetime(2021, 1, 1)
end_date = datetime(2021, 1, 5)
with open("parameters.json", encoding="utf-8") as json_file:
    parameters = json.load(json_file)

# Copy previous simulation files
d3d.copy_simulation_input_files(meteolakes_simulation_folder, output_simulation_folder)
start_date = d3d.copy_simulation_restart_file(meteolakes_simulation_folder, output_simulation_folder, start_date)
d3d.update_control_file(output_simulation_folder, start_date, end_date)

# Get a list of cosmo input files
cosmo_files = d3d.list_cosmo_files(cosmo_folder, start_date, end_date)

# Create river data
d3d.create_river_files(parameters, hydro_folder, start_date, end_date, cosmo_files, output_simulation_folder)

# Create meteo files
#d3d.create_meteo_files(output_simulation_folder, cosmo_files, parameters["grid"])



