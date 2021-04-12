import os
import math
import numpy as np
from scipy.interpolate import griddata, interp1d
import pandas as pd
from shutil import copyfile
from datetime import datetime, timedelta
from hydrosim.functions import log, error
import hydrosim.functions as hf
import netCDF4
import matplotlib.pyplot as plt


def copy_simulation_input_files(input_dir, output_dir):
    log("Copying simulation input files.")
    required_files = ['Simulation_Web.mdf',
                      'Bathymetry.dep',
                      'Grid.grd',
                      'GridEnclosure.enc',
                      'ObservationPoints.obs',
                      "RiversOperations.src"]
    for file in required_files:
        if os.path.isfile(os.path.join(input_dir, file)):
            copyfile(os.path.join(input_dir, file), os.path.join(output_dir, file))
        else:
            log('Failed to copy. Required file: "'+file+'" is not in the input directory.', 1)
    log("Completed copying simulation restart file.")


def copy_simulation_restart_file(input_dir, output_dir, start_date):
    log("Copying simulation restart file.")
    input_files = os.listdir(input_dir)
    restart_dates = []
    for i in range(len(input_files)):
        if "tri-rst.Simulation_Web_rst" in input_files[i]:
            date_str = input_files[i].split(".")[2]
            date = datetime.strptime(date_str, '%Y%m%d')
            restart_dates.append(date)
    closest_date = hf.closest_previous(restart_dates, start_date)
    if closest_date != start_date:
        log("No restart file for " + start_date.strftime(
              "%Y%m%d") + " using restart file from " + closest_date.strftime("%Y%m%d"), 1)
    restart_file = 'tri-rst.Simulation_Web_rst.'+closest_date.strftime("%Y%m%d")+'.000000'
    if os.path.isfile(os.path.join(input_dir, restart_file)):
        copyfile(os.path.join(input_dir, restart_file), os.path.join(output_dir, restart_file))
    else:
        log('Failed to copy. Restart file: "' + restart_file + '" is not in the input directory.', 1)
    log("Completed copying simulation restart file.")
    return closest_date


def update_control_file(folder, start_date, end_date, origin=datetime(2008, 3, 1), period=180):
    with open(os.path.join(folder, "Simulation_Web.mdf"), 'r') as f:
        lines = f.readlines()
    start = "{:.7e}".format((start_date - origin).total_seconds() / 60)
    end = "{:.7e}".format(((end_date - origin).total_seconds() / 60)-period)

    for i in range(len(lines)):
        if "Tstart" in lines[i]:
            lines[i] = "Tstart = " + start + "\n"
        if "Tstop" in lines[i]:
            lines[i] = "Tstop = " + end + "\n"
        if "Restid" in lines[i]:
            lines[i] = "Restid =  #Simulation_Web_rst.{}.000000#\n".format(start_date.strftime("%Y%m%d"))
        if lines[i].split(" ")[0] in ["Flmap", "Flhis", "Flwq"]:
            lines[i] = "{} = {} {} {}\n".format(lines[i].split(" ")[0], start, str(period), end)

    with open(os.path.join(folder, "Simulation_Web.mdf"), 'w') as f:
        f.writelines(lines)


def list_cosmo_files(dir, start_date, end_date, template="cosmo2_epfl_lakes_"):
    log("Listing cosmo files.")
    cosmo_files = []
    cosmo_dates = []
    cosmo_dates_str = []
    for path, subdirs, files in os.walk(dir):
        for name in files:
            if template in name:
                date = datetime.strptime(name.split(".")[0].split("_")[-1], '%Y%m%d')
                if start_date <= date <= end_date:
                    cosmo_files.append(os.path.join(path, name))
                    cosmo_dates.append(date)
                    cosmo_dates_str.append(date.strftime("%Y%m%d"))
    # Check coverage
    for day in np.arange(start_date, end_date, timedelta(days=1)).astype(datetime):
        if day.strftime("%Y%m%d") not in cosmo_dates_str:
            error("COSMO data does not cover full simulation period.")

    # Sort files in order
    df = pd.DataFrame(list(zip(cosmo_dates, cosmo_files)), columns=['dates', 'files'])
    df = df.sort_values(by=['dates'])
    log("Completed listing cosmo files.")
    return list(df["files"])


def create_timeseries(start_date, end_date, step):
    return np.arange(start_date, end_date, timedelta(seconds=step)).astype(datetime)


def cosmo_point_timeseries(coordinates, parameter, files):
    log("Extracting " + parameter + " from COSMO files.")
    out = [[] for _ in range(len(coordinates))]
    time = []
    for file in files:
        log("Processing file " + file, 1)
        date = datetime.strptime(file.split(".")[-2].split("_")[-1], '%Y%m%d')
        time = time + list(np.arange(date, date + timedelta(days=1), timedelta(hours=1)).astype(datetime))
        nc = netCDF4.Dataset(file, mode='r', format='NETCDF4_CLASSIC')
        lat = nc.variables["lat_1"][:]
        lon = nc.variables["lon_1"][:]

        for i in range(len(coordinates)):
            dist = np.sqrt((lat - coordinates[i][0]) ** 2 + (lon - coordinates[i][1]) ** 2)
            min_dist = np.argwhere(dist == np.min(dist))[0]
            out[i] = out[i] + list(np.array(nc.variables[parameter][:, min_dist[0], min_dist[1]]) - 273.15)
        nc.close()

    pt = []
    for i in range(len(coordinates)):
        df = pd.DataFrame(list(zip(time, out[i])), columns=['datetime', parameter])
        df = df.sort_values(by=['datetime'])
        pt.append(df)
    log("Completed extracting " + parameter + " from COSMO files.")
    return pt


def latlng_to_ch1900(lat, lng):
    lat = lat * 3600
    lng = lng * 3600
    lat_aux = (lat - 169028.66) / 10000
    lng_aux = (lng - 26782.5) / 10000
    x = 2600072.37 + 211455.93 * lng_aux - 10938.51 * lng_aux * lat_aux - 0.36 * lng_aux * lat_aux ** 2 - 44.54 * lng_aux ** 3 - 2000000
    y = 1200147.07 + 308807.95 * lat_aux + 3745.25 * lng_aux ** 2 + 76.63 * lat_aux ** 2 - 194.56 * lng_aux ** 2 * lat_aux + 119.79 * lat_aux ** 3 - 1000000
    return x, y


def create_meteo_files(dir, files, grid, origin=datetime(2008, 3, 1)):
    log("Creating meteo files")
    dict = [
        {"filename": 'CloudCoverage.amc', "parameter": "CLCT", "quantity": "cloudiness", "unit": "%", "adjust": 0},
        {"filename": 'Pressure.amp', "parameter": "PMSL", "quantity": "air_pressure", "unit": "Pa", "adjust": 0},
        {"filename": 'RelativeHumidity.amr', "parameter": "RELHUM_2M", "quantity": "relative_humidity", "unit": "%", "adjust": 0},
        {"filename": 'ShortwaveFlux.ams', "parameter": "GLOB", "quantity": "sw_radiation_flux", "unit": "W/m2", "adjust": 0},
        {"filename": 'Temperature.amt', "parameter": "T_2M", "quantity": "air_temperature", "unit": "Celsius", "adjust": -273.15},
        {"filename": 'WindU.amu', "parameter": "U", "quantity": "x_wind", "unit": "m s-1", "adjust": 0},
        {"filename": 'WindV.amv', "parameter": "V", "quantity": "y_wind", "unit": "m s-1", "adjust": 0},
    ]

    log("Read latitude and longitude", 1)
    nc = netCDF4.Dataset(files[0], mode='r', format='NETCDF4_CLASSIC')
    lat = nc.variables["lat_1"][:]
    lon = nc.variables["lon_1"][:]
    nc.close()

    log("Reduce size of MeteoSwiss dataset", 1)
    mx, my = latlng_to_ch1900(lat, lon)
    mxx, myy = mx.flatten(), my.flatten()
    mask = np.logical_and(np.logical_and(mxx >= grid["minx"] - 3 * grid["dx"], mxx <= grid["maxx"] + 3 * grid["dx"]),
                          np.logical_and(myy >= grid["miny"] - 3 * grid["dy"], myy <= grid["maxy"] + 3 * grid["dy"]))
    ind = np.where(mask)
    mxxx, myyy = mxx[ind], myy[ind]

    log("Creating the grid", 1)
    gx = np.arange(grid["minx"], grid["maxx"] + grid["dx"], grid["dx"])
    gy = np.arange(grid["miny"], grid["maxy"] + grid["dy"], grid["dy"])
    gxx, gyy = np.meshgrid(gx, gy)

    log("Create the output meteo files and write the header", 1)
    write_files = []
    for i in range(len(dict)):
        f = open(os.path.join(dir, dict[i]["filename"]), "w")
        f.write('FileVersion = 1.03')
        f.write('\nfiletype = meteo_on_equidistant_grid')
        f.write('\nNODATA_value = -999.00')
        f.write('\nn_cols = '+ str(len(gx)))
        f.write('\nn_rows = '+str(len(gy)))
        f.write('\ngrid_unit = m')
        f.write('\nx_llcenter = '+str(gx[0]))
        f.write('\ny_llcenter = '+str(gy[0]))
        f.write('\ndx = '+str(grid["dx"]))
        f.write('\ndy = '+str(grid["dy"]))
        f.write('\nn_quantity = 1')
        f.write('\nquantity1 = '+dict[i]["quantity"])
        f.write('\nunit1 = '+dict[i]["unit"]+'\n')
        write_files.append(f)

    log("Extracting data from files.", 1)
    for file in files:
        log("Processing file " + file, 2)
        nc = netCDF4.Dataset(file, mode='r', format='NETCDF4_CLASSIC')
        date = datetime.strptime(file.split(".")[-2].split("_")[-1], '%Y%m%d')
        for i in range(len(dict)):
            log("Processing parameter " + dict[i]["parameter"], 3)
            var = nc.variables[dict[i]["parameter"]][:]
            var = var + dict[i]["adjust"]
            nDims = len(nc.variables[dict[i]["parameter"]].dimensions)
            for j in range(var.shape[0]):
                dt = date + timedelta(hours=j)
                diff = dt - origin
                hours = diff.total_seconds() / 3600
                time_str = "TIME = "+str(hours)+"0 hours since "+origin.strftime("%Y-%m-%d %H:%M:%S") + " +00:00"
                write_files[i].write(time_str)
                if nDims == 3:
                    data = var[j, :, :]
                elif nDims == 4:
                    data = var[j, 0, :, :]
                else:
                    error("Incorrect number of dimensions for variable "+dict[i]["parameter"]+" in file "+file)
                grid_interp = griddata((mxxx, myyy), data.flatten()[ind], (gxx, gyy))
                grid_interp[np.isnan(grid_interp)] = -999.00
                write_files[i].write("\n")
                np.savetxt(write_files[i], grid_interp, fmt='%.2f')
    nc.close()

    for i in range(len(write_files)):
        write_files[i].close()
    log("Completed creating meteo files")


def create_river_files(parameters, folder, start_date, end_date, cosmo_files, output_simulation_folder):
    parameters = get_raw_river_data(parameters, folder, start_date, end_date)
    parameters = clean_raw_river_data(parameters, start_date, end_date)
    parameters = flow_balance_ungauged_rivers(parameters)
    parameters = estimate_flow_temperature(parameters, cosmo_files, start_date)
    # verify_against_meteolakes_all(parameters)
    # verify_against_meteolakes(parameters)
    write_river_data_to_file(parameters, output_simulation_folder)


def get_raw_river_data(parameters, folder, start_date, end_date, date_format="%d.%m.%Y%H:%M"):
    log("Getting raw river data.")
    for inflow in parameters["inflows"]:
        if "folder" in inflow:
            inflow["files"] = []
            datetime_arr = []
            temperature_arr = []
            flow_arr = []
            ideal_files = []
            for day in np.arange(start_date, end_date, timedelta(days=1)).astype(datetime):
                ideal_files.append(inflow["prefix"] + day.strftime('%Y%m%d') + ".txt")
            for path, subdirs, files in os.walk(folder):
                for name in files:
                    if name in ideal_files:
                        inflow["files"].append(os.path.join(path, name))
                        df = pd.read_csv(os.path.join(path, name), sep='\t')
                        df["datetime"] = pd.to_datetime(df['Date'] + df["Time"], format=date_format)
                        datetime_arr = datetime_arr + list(df["datetime"])
                        if "temperature" in inflow:
                            temperature_arr = temperature_arr + list(df[inflow["temperature"]])
                        else:
                            temperature_arr = temperature_arr + [np.nan] * len(df["datetime"])
                        if "flow" in inflow:
                            flow_arr = flow_arr + list(df[inflow["flow"]])
                        else:
                            flow_arr = flow_arr + [np.nan] * len(df["datetime"])
            inflow["data"] = pd.DataFrame(zip(datetime_arr, flow_arr, temperature_arr), columns=["datetime", "flow", "temperature"])

    if "outflow" in parameters:
        parameters["outflow"]["files"] = []
        datetime_arr = []
        temperature_arr = []
        flow_arr = []
        ideal_files = []
        for day in np.arange(start_date, end_date, timedelta(days=1)).astype(datetime):
            ideal_files.append(parameters["outflow"]["prefix"] + day.strftime('%Y%m%d') + ".txt")
        for path, subdirs, files in os.walk(folder):
            for name in files:
                if name in ideal_files:
                    parameters["outflow"]["files"].append(os.path.join(path, name))
                    df = pd.read_csv(os.path.join(path, name), sep='\t')
                    df["datetime"] = pd.to_datetime(df['Date'] + df["Time"], format=date_format)
                    datetime_arr = datetime_arr + list(df["datetime"])
                    if "temperature" in parameters["outflow"]:
                        temperature_arr = temperature_arr + list(df[parameters["outflow"]["temperature"]])
                    else:
                        temperature_arr = temperature_arr + [np.nan] * len(df["datetime"])
                    if "flow" in parameters["outflow"]:
                        flow_arr = flow_arr + list(df[parameters["outflow"]["flow"]])
                    else:
                        flow_arr = flow_arr + [np.nan] * len(df["datetime"])
        parameters["outflow"]["data"] = pd.DataFrame(zip(datetime_arr, flow_arr, temperature_arr), columns=["datetime", "flow", "temperature"])

    if "waterlevel" in parameters:
        parameters["waterlevel"]["files"] = []
        datetime_arr = []
        level_arr = []
        ideal_files = []
        for day in np.arange(start_date, end_date, timedelta(days=1)).astype(datetime):
            ideal_files.append(parameters["waterlevel"]["prefix"] + day.strftime('%Y%m%d') + ".txt")
        for path, subdirs, files in os.walk(folder):
            for name in files:
                if name in ideal_files:
                    parameters["waterlevel"]["files"].append(os.path.join(path, name))
                    df = pd.read_csv(os.path.join(path, name), sep='\t')
                    df["datetime"] = pd.to_datetime(df['Date'] + df["Time"], format=date_format)
                    datetime_arr = datetime_arr + list(df["datetime"])
                    level_arr = level_arr + list(df[parameters["waterlevel"]["level"]])
        parameters["waterlevel"]["data"] = pd.DataFrame(zip(datetime_arr, level_arr), columns=["datetime", "level"])
    log("Completed getting raw river data.")
    return parameters


def clean_raw_river_data(parameters, start_date, end_date, dt=timedelta(minutes=10), interpolate="linear"):
    log("Cleaning raw river data (interpolate missing values and smooth).")
    master_datetime = np.arange(start_date, end_date, dt)
    df = pd.DataFrame(master_datetime, columns=["datetime"])
    for inflow in parameters["inflows"]:
        if "folder" in inflow:
            data = df.merge(inflow["data"], on='datetime', how='left')
            data = data.interpolate(method=interpolate)
            inflow["data"] = data

    if "outflow" in parameters:
        data = df.merge(parameters["outflow"]["data"], on='datetime', how='left')
        data = data.interpolate(method=interpolate)
        data["raw_flow"] = data["flow"]
        data["flow"] = data["flow"].rolling(window=12, win_type='gaussian', center=True).mean(std=2)
        data = data.interpolate(method=interpolate, limit_direction='backward')
        data = data.interpolate(method=interpolate)
        parameters["outflow"]["data"] = data

    if "waterlevel" in parameters:
        data = df.merge(parameters["waterlevel"]["data"], on='datetime', how='left')
        data = data.interpolate(method=interpolate)
        data["level"] = data["level"].rolling(window=144, win_type='gaussian', center=True).mean(std=32)
        data = data.interpolate(method=interpolate, limit_direction='backward')
        data = data.interpolate(method=interpolate)
        parameters["waterlevel"]["data"] = data
    log("Completed cleaning raw river data.")
    return parameters


def get_water_level_change_flow_equivalent(df_wl, altitude, bathymetry):
    df_wl['level+'] = df_wl['level'].shift(-1)
    df_wl.iloc[-1, df_wl.columns.get_loc('level+')] = df_wl['level'].iloc[-1]
    depth_area = interp1d(bathymetry["depth"], bathymetry["area"])
    df_wl["area"] = depth_area(df_wl['level'] - altitude)
    df_wl["area+"] = depth_area(df_wl['level+'] - altitude)
    df_wl["area_mean"] = df_wl[['area', 'area+']].mean(axis=1)
    df_wl["dV"] = (df_wl['level+'] - df_wl['level']) * df_wl["area_mean"]
    dt = (df_wl["datetime"].iloc[1] - df_wl["datetime"].iloc[0]).total_seconds()
    df_wl["flow"] = df_wl["dV"] / dt
    return np.array(df_wl["flow"])


def verify_against_meteolakes(parameters):
    df = pd.read_csv("comp.csv")
    for inflow in parameters["inflows"]:
        ax = plt.subplot(211)
        ax.plot(inflow["data"]["flow"], label="james_" + inflow["name"])
        ax.plot(df[inflow["name"]], label="theo_" + inflow["name"])
        ax.legend()
        ax = plt.subplot(212)
        ax.plot((df[inflow["name"]] - inflow["data"]["flow"]) / df[inflow["name"]])
        plt.show()

    ax = plt.subplot(211)
    ax.plot(parameters["outflow"]["data"]["flow"], label="james outflow")
    ax.plot(df["Out"], label="theo")
    ax.legend()
    ax = plt.subplot(212)
    ax.plot((df["Out"] - parameters["outflow"]["data"]["flow"]) / df["Out"])
    plt.show()

    ax = plt.subplot(211)
    ax.plot(parameters["waterlevel"]["data"]["level"], label="james_level")
    ax.plot(df["Level"], label="theo")
    ax = plt.subplot(212)
    ax.plot((df["Level"] - parameters["waterlevel"]["data"]["level"]) / df["Level"])
    plt.show()


def verify_against_meteolakes_all(parameters):
    df = pd.read_csv("comp.csv")
    for inflow in parameters["inflows"]:
        plt.plot(inflow["data"]["flow"], label="james_" + inflow["name"])
        plt.plot(df[inflow["name"]], label="theo_" + inflow["name"])

    plt.plot(parameters["outflow"]["data"]["flow"], label="james_outflow")
    plt.plot(df["Out"], label="theo_outflow")

    plt.legend()
    plt.show()


def flow_balance_ungauged_rivers(parameters, min_flow=0.0001, dt=600):
    log("Calculating flow balance to assign flow rates to ungauged rivers")
    outflow = np.array(parameters["outflow"]["data"]["flow"])
    master_datetime = np.array(parameters["outflow"]["data"]["datetime"])
    temperature = [np.nan] * len(master_datetime)
    water_level_flow = get_water_level_change_flow_equivalent(parameters["waterlevel"]["data"], parameters["altitude"], parameters["bathymetry"])
    extra_flow = outflow + water_level_flow

    for inflow in parameters["inflows"]:
        if "folder" in inflow:
            extra_flow = extra_flow - np.array(inflow["data"]["flow"])
    for inflow in parameters["inflows"]:
        if "folder" not in inflow:
            flow = inflow["contribution"] * extra_flow
            flow[np.isnan(flow)] = 0
            neg_flow = np.copy(flow)
            neg_flow[neg_flow > min_flow] = 0
            added_vol = np.abs(np.sum(neg_flow * dt))
            flow[flow <= min_flow] = min_flow
            total_vol = np.sum(flow * dt)
            flow = flow * (1 - added_vol/total_vol)
            inflow["data"] = pd.DataFrame(zip(master_datetime, flow, temperature), columns=["datetime", "flow", "temperature"])

    log("Completed calculating flow balance.")
    return parameters


def toffolon_air_to_water_temperature(start_temperature, air_temperature, flow, a, dt=600, ty=365.2422):
    a = np.array(a) / (3600 * 24)
    water_temperature = np.array([start_temperature] * len(air_temperature))
    """for i in range(1, len(air_temperature)):
        theta = np.abs(flow[i] / np.nanmean(flow))
        delta = theta * a[3]
        t = dt / (3600 * 24)
        change_temperature = 1 / delta * (a[0] + a[1] * air_temperature[i] - a[2] * water_temperature[i - 1] + theta * (a[4] + a[5] * math.cos(2 * math.pi * (t / ty - a[6])) - a[7] * water_temperature[i - 1])) * dt
        water_temperature[i] = water_temperature[i] + change_temperature

    print(water_temperature)
    plt.plot(water_temperature)
    plt.plot(air_temperature)
    plt.show()
    exit()"""
    return water_temperature


def estimate_flow_temperature(parameters, files, start_date, temperature="T_2M"):
    log("Estimating river temperatures")
    coordinates = list(map(lambda x: x["coordinates"], parameters["inflows"]))
    air_temperature = cosmo_point_timeseries(coordinates, temperature, files)

    for i in range(len(parameters["inflows"])):
        if parameters["inflows"][i]["data"]["temperature"].isnull().values.any():
            df = parameters["inflows"][i]["data"].merge(air_temperature[i], on="datetime", how="left")
            df = df.interpolate(method="linear")
            start_temperature = parameters["inflows"][i]["temperature"][start_date.month - 1]
            df["temperature"] = toffolon_air_to_water_temperature(start_temperature, df[temperature],
                                                                  parameters["inflows"][i]["data"]["flow"],
                                                                  parameters["inflows"][i]["a"])
            parameters["inflows"][i]["data"] = df
    log("Completed estimating river temperatures")
    return parameters


def write_river_data_to_file(parameters, folder, filename="RiversOperationsQuantities.dis", origin=datetime(2008, 3, 1)):
    log("Writing river data to file.")
    f = open(os.path.join(folder, filename), "w")
    for inflow in parameters["inflows"]:
        f.write("table-name           'Discharge : {}'\n".format(str(inflow["id"])))
        f.write("contents             'momentum'\n")
        f.write("location             '{}'\n".format(inflow["name"]))
        f.write("time-function        'non-equidistant'\n")
        f.write("reference-time       {}\n".format(origin.strftime("%Y%m%d")))
        f.write("time-unit            'minutes'\n")
        f.write("interpolation        'linear'\n")
        f.write("parameter            'time                '                     unit '[min]'\n")
        f.write("parameter            'flux/discharge rate '                     unit '[m3/s]'\n")
        f.write("parameter            'Temperature         '                     unit '[°C]'\n")
        f.write("parameter            'flow magnitude      '                     unit '[m/s]'\n")
        f.write("parameter            'flow direction      '                     unit '[deg]'\n")
        f.write("records-in-table     {}\n ".format(str(len(inflow["data"]))))

        df = inflow["data"]
        df["flow_velocity"] = inflow["flow_velocity"]
        df["flow_direction"] = inflow["flow_direction"]
        df["minutes"] = (df["datetime"]-origin).astype('timedelta64[m]')
        out = np.column_stack((df["minutes"], df["flow"], df["temperature"], df["flow_velocity"], df["flow_direction"]))
        np.savetxt(f, out, fmt="%.7e", delimiter="   ", newline="\n ")

    f.write("table-name           'Discharge : {}'\n".format(str(parameters["outflow"]["id"])))
    f.write("contents             'regular  '\n")
    f.write("location             '{}           '\n".format(parameters["outflow"]["name"]))
    f.write("time-function        'non-equidistant'\n")
    f.write("reference-time       {}\n".format(origin.strftime("%Y%m%d")))
    f.write("time-unit            'minutes'\n")
    f.write("interpolation        'linear'\n")
    f.write("parameter            'time                '                     unit '[min]'\n")
    f.write("parameter            'flux/discharge rate '                     unit '[m3/s]'\n")
    f.write("parameter            'Temperature         '                     unit '[°C]'\n")
    f.write("records-in-table     {}\n ".format(str(len(parameters["outflow"]["data"]))))

    df = parameters["outflow"]["data"]
    df["temperature"] = 6
    df["minutes"] = (df["datetime"] - origin).astype('timedelta64[m]')
    out = np.column_stack((df["minutes"], -df["flow"], df["temperature"]))
    np.savetxt(f, out, fmt="%.7e", delimiter="   ", newline="\n ")

    f.close()
    log("Completed writing river data to file.")
