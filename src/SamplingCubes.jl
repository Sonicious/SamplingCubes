module SamplingCubes

using Blosc, Zarr, YAXArrays
Blosc.set_num_threads(Threads.nthreads())
using Interpolations # for Constant Interpolation
using Images, GLMakie, Colors, ColorSchemes # everything plotting related
using Statistics, StatsBase

const global CUBEPATH = Ref{String}()
const global CUBE = Ref{YAXArray}()

"""
    __init__

initializes the cubes which are empty from the beginning
"""
function __init__()
    # Getting roughly Europe
    arealongitude = (-10, 35)
    arealatitude = (35, 65)
    CUBEPATH[] = "D:/ERA5Data.zarr"
    CUBE[] = Cube(CUBEPATH[])
    CUBE[] = CUBE[][longitude=arealongitude, latitude=arealatitude, Variable=["t2m"]]
    lon, lat, data = GeoDatasets.landseamask(grid=1.25, resolution='f')
    data_lons = findall(x -> arealongitude[1] <= x <= arealongitude[2], lon)
    data_lats = findall(x -> arealatitude[1] <= x <= arealatitude[2], lat)
    mask = imresize(data[data_lons, data_lats], size(CUBE[])[1:2], method=Constant())
    mask .= reverse(mask, dims=2)
    mask .= mask .> 0
    const global LSM = mask
end

"""
"""
function initialPoints(strategy::String, lsm::Bool, value)
    initialpositions = Point2f[]
    if strategy == "random"
        randompoints = value
        if lsm
            for i in 1:randompoints
                while true
                    newpoint = (sample(lonSize, 1)..., sample(latSize, 1)...)
                    if mask[newpoint...] > 0
                        push!(initialpositions, Point2f(longitudevalues[newpoint[1]], latitudevalues[newpoint[2]]))
                        break
                    end 
                end
            end
        else
            for i in 1:randompoints
                newpoint = (sample(lonSize, 1)..., sample(latSize, 1)...)
                if mask[newpoint...] > 0
                    push!(initialpositions, Point2f(longitudevalues[newpoint[1]], latitudevalues[newpoint[2]]))
                    break
                end 
            end
        end
    elseif strategy == "regular"
    else
    end
end

longitudevalues = collect(samplecube.longitude.values)
latitudevalues = collect(samplecube.latitude.values)

sampling = "regular"

lonSize = 1:size(longitudevalues, 1)
latSize = 1:size(latitudevalues, 1)

if sampling == "random"
    randompoints = 150
    for i in 1:randompoints
        while true
            newpoint = (sample(lonSize, 1)..., sample(latSize, 1)...)
            if mask[newpoint...] > 0
                push!(landpositions, Point2f(longitudevalues[newpoint[1]], latitudevalues[newpoint[2]]))
                break
            end
        end
    end
end

if sampling == "regular"
    griddistance = 10
    for x in 1:griddistance:lonSize, y in 1:griddistance:latSize
        newpoint = (x, y)
        if mask[newpoint...] > 0
            push!(landpositions, Point2f(longitudevalues[newpoint[1]], latitudevalues[newpoint[2]]))
        end
    end
end

data = mask;

fig = Figure()
ax = GLMakie.Axis(fig[1, 1])
ax.aspect = DataAspect()
hm = heatmap!(ax, longitudevalues, latitudevalues, mask, colorrange=(0, 1))
scatter!.(landpositions, color=:red, marker=:x, markersize=10)

for i in CartesianIndices((3, 3))
    
end

numClasses = size(landpositions, 1) + 1
cm = distinguishable_colors(numClasses + 1, lchoices=30:15:130)
cm[1] = RGB(0, 0, 0)

end