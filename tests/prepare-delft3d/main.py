from datetime import datetime
import hydrosim.delft3d as d3d
import json


# Define parameters
meteolakes_simulation_folder = "./example-sim"
output_simulation_folder = "./output-sim"
cosmo_folder = "./demo-data/Meteodata"
start_date = datetime(2021, 1, 1)
end_date = datetime(2021, 1, 4)
grid = {"minx": 490000, "miny": 110000, "maxx": 570000, "maxy": 160000, "dx": 1000, "dy": 1000}
with open("parameters.json") as json_file:
    parameters = json.load(json_file)

# Copy previous simulation files
d3d.copy_simulation_input_files(meteolakes_simulation_folder, output_simulation_folder)
start_date = d3d.copy_simulation_restart_file(meteolakes_simulation_folder, output_simulation_folder, start_date)

# Get a list of cosmo input files
cosmo_files = d3d.list_cosmo_files(cosmo_folder, start_date, end_date)

# Create river data
coordinates = list(map(lambda x: x["coordinates"], parameters["inflows"]))
air_temperature = d3d.cosmo_point_timeseries(coordinates, "T_2M", cosmo_files)



# Create meteo files
#d3d.create_meteo_files(output_simulation_folder, cosmo_files, grid)



