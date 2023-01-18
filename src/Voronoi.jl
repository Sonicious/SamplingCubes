# part of the SamplingCube module

export InitialPoints
export createDistance

export BruteVoronoi, JflVoronoi

export Brutetest, JFLTest, TestInitialPoints

global testpoints = Point2f[
    (22.625, 41.875),
    (7.625, 58.125),
    (-1.375, 49.125),
    (-2.875, 52.125),
    (20.125, 49.125),
    (-3.375, 41.875)
]

#= function VarianceHistogram(event_idx)
    cube = EXTREMECUBES[][event_idx]
    mask = GetEventLandSeaMask(event_idx)
    mm = Array{Union{Float32,Missing}}(undef, size(mask)...)
    mm[mask.>0] .= 1.0
    data = GetEvent(event_idx)

end =#

"""
"""
function BruteTest(numCenters = 200)
    data = GetFilteredEvent(7)
    dots = InitialPoints(7, numCenters)
    img  = zeros(RGB{N0f8}, size(data)[1:2])
    for idx in eachindex(data), pp = 1:numCenters
        distances = dcor(data[idx], dots[pp]) # distance
        nn = findmin(distances)[2]
        img[x,y]  = dots[nn,:][3]
    end
    return img
end

function BruteVoronoi(event_idx, numCenters = 200)
    dist = (point,vector) -> sqrt.((point[1].-vector[:,1]).^2 .+ (point[2].-vector[:,2]).^2)
    dots = [rand(1:height, numCenters) rand(1:width, numCenters) rand(RGB{N0f8}, numCenters)]
    img  = zeros(RGB{N0f8}, height, width)
    for x in 1:height, y in 1:width
        distances = dist([x,y],dots) # distance
        nn = findmin(distances)[2]
        img[x,y]  = dots[nn,:][3]
    end
    return img
end

"""
"""
function JFLTest(height=600, width=800, numCenters = 200)
    dist = (point,vector) -> sqrt.((point[1].-vector[:,1]).^2 .+ (point[2].-vector[:,2]).^2)
    dots = [rand(1:height, numCenters) rand(1:width, numCenters) rand(RGB{N0f8}, numCenters)]
    img  = zeros(RGB{N0f8}, height, width)
    for x in 1:height, y in 1:width
        distances = dist([x,y],dots) # distance
        nn = findmin(distances)[2]
        img[x,y]  = dots[nn,:][3]
    end
    return img
end

"""
"""
function JumpFlooding(initialPoints)
    
end

"""
    findnearest(sortedlist, value)

Helper function: find nearest in sorted list.
"""
function findnearest(sortedlist, value)
    return findmin(abs.(sortedlist .- value))
end

"""
    InitialPoints(event_idx, num_samples)

Calculate the initial points for the creation of the voronoi diagram
"""
function InitialPoints(event_idx, num_samples)
    cube = EXTREMECUBES[][event_idx]
    mask = GetEventLandSeaMask(event_idx)
    samplePositions = sample(findall(mask .== 1), num_samples)

    longitudevalues = collect(cube.longitude.values)
    latitudevalues = collect(cube.latitude.values)
    landpositions = Point2f[]

    for point ∈ samplePositions
        push!(landpositions, (longitudevalues[point[1]], latitudevalues[point[2]]))
    end
    return landpositions
end

"""
"""
function TestInitialPoints()
    cube = EXTREMECUBES[][9]
    mask = GetEventLandSeaMask(9)
    longitudevalues = collect(cube.longitude.values)
    latitudevalues = collect(cube.latitude.values)
    positions = InitialPoints(9, 5)
    mm = Array{Union{Float32,Missing}}(undef, size(mask)...)
    mm[mask.>0] .= 1.0

    fig = Figure()
    ax = GLMakie.Axis(fig[1, 1])
    ax.aspect = DataAspect()

    hm = heatmap!(ax,
        longitudevalues, latitudevalues, mm,
        colorrange=(0, 1)
    )

    scatter!.(positions, color=:red, marker=:x, markersize=30)

    return fig
end

"""
    createDistance(data, position; measure=dcor)

Calculate distances of one point in the data matrix to all other entries. The position is an index position inside the matrix. The measure for distance calculation is a function (a,b)->distance, default is the distance correlation.
"""
function createDistance(data, position; measure=dcor)
    distances = zeros(size(data)[1:2])
    testpoint = data[position..., :]
    for i ∈ axes(data, 1), j ∈ axes(data, 2)
        distances[i, j] = measure(data[i, j, :], testpoint)
    end
    return distances
end


########################################################################
# Figure functions
########################################################################
export FigureCircle, FigureDistance, FigurePoints, FigureLines, FigureBox, FigureDistanceArea

function FigureDistance(event_idx, position; masked=true, measure=dcor)

    cube = EXTREMECUBES[][event_idx]
    mask = GetEventLandSeaMask(event_idx)
    mm = Array{Union{Float32,Missing}}(undef, size(mask)...)
    mm[mask.>0] .= 1.0
    data = GetFilteredEvent(event_idx)

    longitudevalues = collect(cube.longitude.values)
    latitudevalues = collect(cube.latitude.values)
    positionlongitude = findnearest(longitudevalues, position[1])[2]
    positionlatitude = findnearest(latitudevalues, position[2])[2]

    data .= createDistance(data, (positionlongitude, positionlatitude), measure=measure)

    # get min and max data for axis
    datamax = maximum(data)
    datamin = minimum(data)

    fig = Figure()
    ax = GLMakie.Axis(fig[1, 1])
    ax.aspect = DataAspect()

    if !masked
        mm[mask.<=0] .= -1.0
        datamin = -datamax
    end

    hm = heatmap!(ax,
        longitudevalues, latitudevalues, data[:, :, 1] .* mm,
        colorrange=(datamin, datamax),
    )
    Colorbar(fig[:, end+1], hm)
    scatter!([position...]', color=:red, marker=:x, markersize=30)
    return fig
end

function FigureLines()
    cube = EXTREMECUBES[][9]
    longitudevalues = collect(cube.longitude.values)
    latitudevalues = collect(cube.latitude.values)

    position = (-8, 40)

    longitudevalues = collect(cube.longitude.values)
    latitudevalues = collect(cube.latitude.values)
    positionlongitude = findnearest(longitudevalues, position[1])[2]
    positionlatitude = findnearest(latitudevalues, position[2])[2]

    data1 = GetEvent(9)
    data1 = data1[positionlongitude, positionlatitude, :]
    data2 = GetEvent(9)
    data2 = data2[positionlongitude+5, positionlatitude+5, :]

    fig = Figure()
    lines(data1, color=:blue)
    lines!(data2, color=:red)

    return data1, data2
end

function FigureBox()
    mask = GetEventLandSeaMask(9)
    cube = EXTREMECUBES[][9]
    longitudevalues = collect(cube.longitude.values)
    latitudevalues = collect(cube.latitude.values)

    fig = Figure()
    ax = GLMakie.Axis(fig[1, 1])
    ax.aspect = DataAspect()

    hm = heatmap!(ax,
        longitudevalues, latitudevalues, mask,
        colorrange=(0, 1)
    )

    scatter!.(testpoints, color=:red, marker=:x, markersize=40)

    for (index, point) in enumerate(testpoints)
        polygon = Point2f[]
        tt = 1.5
        push!(polygon, point .+ (-tt, -tt))
        push!(polygon, point .+ (-tt, +tt))
        push!(polygon, point .+ (+tt, +tt))
        push!(polygon, point .+ (+tt, -tt))
        poly!(Polygon(polygon))
    end

    return fig
end

function FigurePoints()
    mask = GetEventLandSeaMask(9)
    cube = EXTREMECUBES[][9]
    longitudevalues = collect(cube.longitude.values)
    latitudevalues = collect(cube.latitude.values)

    fig = Figure()
    ax = GLMakie.Axis(fig[1, 1])
    ax.aspect = DataAspect()

    hm = heatmap!(ax,
        longitudevalues, latitudevalues, mask,
        colorrange=(0, 1)
    )

    scatter!.(testpoints, color=:red, marker=:x, markersize=40)

    return fig
end

function FigureCircle()
    mask = GetEventLandSeaMask(9)
    cube = EXTREMECUBES[][9]
    longitudevalues = collect(cube.longitude.values)
    latitudevalues = collect(cube.latitude.values)

    fig = Figure()
    ax = GLMakie.Axis(fig[1, 1])
    ax.aspect = DataAspect()

    hm = heatmap!(ax,
        longitudevalues, latitudevalues, mask,
        colorrange=(0, 1)
    )

    scatter!.(testpoints, color=:red, marker=:x, markersize=40)

    for (index, point) in enumerate(testpoints)
        tt = 1.5
        poly!(Circle(point, tt), Color=:red)
    end

    return fig
end

function FigureDistanceArea()
    mask = GetEventLandSeaMask(9)
    cube = EXTREMECUBES[][9]
    longitudevalues = collect(cube.longitude.values)
    latitudevalues = collect(cube.latitude.values)

    fig = Figure()
    ax = GLMakie.Axis(fig[1, 1])
    ax.aspect = DataAspect()

    hm = heatmap!(ax,
        longitudevalues, latitudevalues, mask,
        colorrange=(0, 1)
    )

    scatter!.(testpoints, color=:red, marker=:x, markersize=40)

    for (index, point) in enumerate(testpoints)
        tt = 1.5
        poly!(Circle(point, tt), Color=:red)
    end

    return fig
end