using Blosc, Zarr, YAXArrays
Blosc.set_num_threads(Threads.nthreads())
using GeoDatasets
using Interpolations # for Constant Interpolation
using Images, GLMakie, Colors, ColorSchemes # everything plotting related
using Statistics, StatsBase

fullcube = Cube("D:/ERA5Data.zarr")
arealongitude = (-10, 35)
arealatitude = (35, 65)
landpositions = Point2f[]
samplecube = fullcube[longitude=arealongitude, latitude=arealatitude, Variable=["t2m"]]

_, _, ls = GeoDatasets.landseamask(grid=1.25, resolution='f')

data_lats = findall(x -> arealatitude[1] <= x <= arealatitude[2], lat)
data_lons = findall(x -> arealongitude[1] <= x <= arealongitude[2], lon)
mask = imresize(ls[data_lons, data_lats], size(samplecube[])[1:2], method=Constant())
mask .= reverse(mask, dims=2)
mask .= mask .> 0

longitudevalues = collect(samplecube.longitude.values)
latitudevalues = collect(samplecube.latitude.values)

sampling = "regular"

lonSize = 1:size(longitudevalues,1)
latSize = 1:size(latitudevalues, 1)

if sampling=="random"
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
        newpoint = (x,y)
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

numClasses = size(landpositions,1)+1
cm = distinguishable_colors(numClasses+1, lchoices=30:15:130)
cm[1] = RGB(0,0,0)

