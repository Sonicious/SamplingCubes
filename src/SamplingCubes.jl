module SamplingCubes

#=
 Todo:
 * config.toml
 * correlation distance
=#

using Blosc, Zarr, YAXArrays
Blosc.set_num_threads(Threads.nthreads())
using DataFrames, Dates, CSV
using GeoDatasets
using Interpolations # for Constant Interpolation
using Images, GLMakie
using ProgressMeter

const global CUBEPATH = Ref{String}()
const global CSVPATH = Ref{String}()

const global EXTREMECUBES = Ref{Vector{YAXArray}}()
const global EXTREMEDATAFRAME = Ref{DataFrame}()

export InitExtremes, GetEvent, GetMaskedEvent, MakeT2mVideo

function __init__()
    CUBEPATH[] = "D:/ERA5Data.zarr"
    CSVPATH[] = "src/EventPart1_csv_Ltime.csv"
end

"""
    InitExtremes()

This function initializes all the data and provides the Data
"""
function InitExtremes()
    # first the csv data
    EXTREMEDATAFRAME[] = CSV.read(CSVPATH[], DataFrame, dateformat=dateformat"yyyy.mm.dd")
    # now the cubes
    fullcube = Cube(CUBEPATH[])
    numcubes = nrow(EXTREMEDATAFRAME[])
    events = Vector{YAXArray}(undef, numcubes)
    for index = 1:nrow(EXTREMEDATAFRAME[],)
        time = (EXTREMEDATAFRAME[][index, :when_from], EXTREMEDATAFRAME[][index, :when_until])
        longitude = Tuple(EXTREMEDATAFRAME[][index, [:where_W, :where_E]])
        latitude = Tuple(EXTREMEDATAFRAME[][index, [:where_S, :where_N]])
        events[index] = fullcube[longitude=longitude, latitude=latitude, time=time, Variable=["t2m"]]
    end
    EXTREMECUBES[] = events
end

# image(ev[:,:,1], axis = (aspect = DataAspect(), title = "Europe Drought",))

"""
    GetEvent(eventidx = 1, celsius = true)

This function returns the t2m-data of the event `eventidx` as pure data, not `YAXArray`. Conversion to degree Celsius can be done as well.
"""
function GetEvent(eventidx=1, celsius=true)
    if (celsius)
        return collect(EXTREMECUBES[][eventidx].data)[:, :, :] .- 273.15
    else
        return collect(EXTREMECUBES[][eventidx].data)[:, :, :]
    end
end

"""
    GetEventLandSeaMask(event_idx)

Returns a mask for the event of the same size as the cube. Sea and Lakes are considered as 0 and landmass as 1
"""
function GetEventLandSeaMask(event_idx)
    eventlons, eventlats = GetEventCoordinates(event_idx)
    lon, lat, data = GeoDatasets.landseamask(grid=1.25, resolution='f')

    data_lats = findall(x -> eventlats[1] <= x <= eventlats[2], lat)
    data_lons = findall(x -> eventlons[1] <= x <= eventlons[2], lon)
    mask = imresize(data[data_lons, data_lats], size(EXTREMECUBES[][event_idx])[1:2], method=Constant())
    mask .= reverse(mask, dims=2)
    mask .= mask .> 0
    return mask
end

"""
    GetEventCoordinates(event_idx)

Helper function for getting the coordinates of the event
"""
function GetEventCoordinates(event_idx)
    dfr = EXTREMEDATAFRAME[][event_idx, :]
    lons = (dfr.where_W, dfr.where_E)
    lats = (dfr.where_S, dfr.where_N)
    return lons, lats
end

"""
    MakeT2mVideo(event_idx)

This function creates a video for the extreme event `event_idx` and stores it as mp4 in the results directory
"""
function MakeT2mVideo(event_idx)
    cube = EXTREMECUBES[][event_idx]
    # get lat/lon coordinates as string for axes
    x = [string(idx) for idx in cube.longitude.values]
    y = [string(idx) for idx in cube.latitude.values]

    dfr = EXTREMEDATAFRAME[][event_idx, :]
    eventname = dfr."Name"
    eventtype = dfr."Event type"
    eventyear = dfr."Year"

    # collect data
    mask = GetEventLandSeaMask(event_idx)
    data = GetEvent(event_idx) .* mask

    # get min and max data for axis
    datamax = maximum(data)
    datamin = minimum(data)
    longitudevalues = collect(cube.longitude.values)
    latitudevalues = collect(cube.latitude.values)

    # construct initial heatmap at timepoint 1
    fig, ax, hm = heatmap(
        longitudevalues, latitudevalues, data[:, :, 1],
        colorrange=(datamin, datamax),
    )
    ax.title = eventname * "-" * eventtype * "-" * string(Date(cube.time.values[1]))
    Colorbar(fig[:, end+1], hm)

    filename = "./results/" * eventname * "-" * eventtype * "-" * string(eventyear) * ".mp4"

    # loop through all time-points in 2D time-series beginning at idx 2
    timepoints = range(2, size(data, 3))
    record(fig, filename, timepoints; framerate=10) do timepoint
        hm[3][] = data[:, :, timepoint]
        ax.title = eventname * "-" * eventtype * "-" * string(Date(cube.time.values[timepoint]))
    end

end

end # module
