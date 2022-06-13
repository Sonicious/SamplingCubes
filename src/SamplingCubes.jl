module SamplingCubes

using EarthDataLab
using YAXArrays
using GeoDatasets
using Plots
using Interpolations # for Constant Interpolation
using Images
using DataFrames
using Dates: DateTime
using CSV
using ProgressMeter

global const CUBEPATH = Ref{String}()
global const CSVPATH = Ref{String}()

global const EXTREMECUBES = Ref{Vector{YAXArray}}()
global const EXTREMEDATAFRAME = Ref{DataFrame}()

export InitExtremes, GetExtremeEventData, GetExtremeEventCubes, GetEventCoordinates, GetEventLandSeaMask, GetMaskedEvent, GetSummedValues

function __init__()
    gr()
    CUBEPATH[] = "D:/deep_extremes/ERA5Data.zarr"
    CSVPATH[] = "src/EventPart1_csv_Ltime.csv"
end

function InitExtremes()
    # first the csv data
    tempframe = CSV.read(CSVPATH[], DataFrame)
    str_to_date(date) = DateTime(date, "yyyy.mm.dd")
    tempframe[!, :when_from] = str_to_date.(tempframe[!, :when_from])
    tempframe[!, :when_until] = str_to_date.(tempframe[!, :when_until])
    EXTREMEDATAFRAME[] = tempframe
    # now the cubes
    fullcube = Cube(CUBEPATH[])
    numcubes = size(EXTREMEDATAFRAME[], 1)
    events = Vector{YAXArray}(undef, numcubes)
    for index = 1:size(tempframe, 1)
        time = (tempframe[index, :when_from], tempframe[index, :when_until])
        longitude = Tuple(tempframe[index, [:where_W, :where_E]])
        latitude = Tuple(tempframe[index, [:where_S, :where_N]])
        events[index] = subsetcube(fullcube, lon=longitude, lat=latitude, time=time, Variable=["t2m", "t2mmin", "t2mmax"])
    end
    EXTREMECUBES[] = events
    return nothing
end

function GetExtremeEventData(event_idx)
    return collect(EXTREMECUBES[][event_idx].data)
end

function GetExtremeEventCubes()
    return EXTREMECUBES[]    
end

function GetEventCoordinates(event_idx)
    dfr = EXTREMEDATAFRAME[][event_idx, :]
    lons =  (dfr.where_W, dfr.where_E)
    lats = (dfr.where_S, dfr.where_N)
    return lons, lats
end

function GetEventLandSeaMask(event_idx)
    eventlons, eventlats = GetEventCoordinates(event_idx)
    lon, lat, data = GeoDatasets.landseamask(grid=1.25, resolution='f')

    data_lats = findall(x -> eventlats[1] <= x <= eventlats[2], lat)
    data_lons = findall(x -> eventlons[1] <= x <= eventlons[2], lon)
    mask = imresize(data[data_lons, data_lats]', reverse(size(EXTREMECUBES[][event_idx])[1:2]), method=Constant())
    mask .= mask .> 0
    return mask
end

function GetMaskedEvent(event_idx::Int64; datapoint=1, normalize=true)
    mask = GetEventLandSeaMask(event_idx)
    data = collect(subsetcube(EXTREMECUBES[][event_idx],Variable="t2mmax").data)
    if normalize
        data .-= 273.15
    end
    data = reverse(data[:, :, datapoint], dims=2)'
    return mask .* data
end

function GetSummedValues(event_idx)
    mask = GetEventLandSeaMask(event_idx)
    data = collect(subsetcube(EXTREMECUBES[][event_idx], Variable="t2mmax").data)
    data = dropdims(sum(data, dims=3), dims=3)
    data = reverse(data, dims=2)'
    mask .* data
end

function ShowCubeHistogram
end

function ShowSamplingCube
end

function t2mmax_video(event_idx)
    mask = GetEventLandSeaMask(event_idx)
    tempcube = subsetcube(EXTREMECUBES[][event_idx], Variable="t2mmax")
    x = [string(idx) for idx in tempcube.longitude.values]
    y = [string(idx) for idx in tempcube.latitude.values]

    dfr = EXTREMEDATAFRAME[][event_idx,:]
    eventname = dfr."Name"
    eventtype = dfr."Event type"
    eventyear = dfr."Year"

    datamax = maximum(collect(tempcube.data)) - 273.15
    datamin = minimum(collect(tempcube.data)) - 273.15

    animation = Animation()

    @showprogress for timepoint in tempcube.time.values
        titlestring = eventname*"-"*eventtype*"-"*string(Date(timepoint))
        data = reverse(collect(tempcube[time=timepoint].data) .- 273.15, dims=2)'
        data .= mask .* data
        plot = heatmap(x, y, data,
            title=titlestring,
            clim=(datamin, datamax),
        )
        frame(animation)
    end

    gif(animation, "./results/"*eventname*"-"*eventtype*"-"*string(eventyear)*".gif", fps=5)
end

end # module
