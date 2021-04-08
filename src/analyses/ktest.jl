using CSV
using DataFrames
using Distances: Euclidean
pointData = CSV.read("../Data/Stem_Map.csv", DataFrame)
my_win = ppWindow((-150, 150), (-150, 150))
#### Testing Data
my_array = Array{Union{Nothing, Int}}(nothing, 2, 5)
my_pp=CartesianPointPattern((-100,100),my_array)

abstract type PointPattern end

struct ppWindow
    xdim::Tuple{Float32, Float32}
    ydim::Tuple{Float32, Float32}
end

my_win = ppWindow((-150, 150), (-150, 150))
my_points = pointData[:,2:3]

function ppWin_area(win::ppWindow)
    x_range = abs(win.xdim[1]-win.xdim[2])
    y_range = abs(win.yrange[1]-win.yrange[2])
    area = x_range * y_range
    return area
end

struct CartesianPointPattern <: PointPattern
    #Use with basic x,y coordinate data
    pattern::Union{Array,DataFrame}
    window::ppWindow
    marks::Union{Vector,Nothing}
end

struct GeographicPointPattern <: PointPattern
    #Possible interfaces with geographic projections
myPPP.pattern[:,1]

myPPP = CartesianPointPattern(my_win, pointData[:,(2:3)], pointData[:,(4)])

function AsPointPattern(data::Any, win::Any)
    try
        #df = convert(data, DataFrame)
        if size(data,2) > 2 
            pp_pattern = df[:,(0:1)]
            pp_marx = df[:,Not(0:1)]
        elseif size(x,2) = 2
            pp_pattern = df[:,(0:1)]
            pp_marx = nothing
        else
            throw(DomainError(x, "Input must have at least 2 columns"))
        end
        pp_win = convert(win,Tuple{Float32,Float32})
        pp_out = CartesianPointPattern(pp_win, pp_pattern, pp_marx)
        return pp_out
    catch e
        if isa(e, DomainError)


    end
    df = convert(data, DataFrame)
    if size(data,2) > 2 
        pp_pattern = df[:,(0:1)]
        pp_marx = df[:,Not(0:1)]
    elseif size(x,2) = 2
        pp_pattern = df[:,(0:1)]
        pp_marx = nothing
    else
        throw(DomainError(x, "Input must have at least 2 columns"))
    end
    pp_win = convert(win,Array{Tuple{Float32,Float32}, Tuple{Float32,Float32}})
    pp_out = CartesianPointPattern(pp_win, pp_pattern, pp_marx)
    return pp_out
end

function distance_calc(pp::CartesianPointPattern)
    output_dim = length(pp.pattern[:,1]) #finds size of output Array
    dist_matrix = Array{Float64}(undef, output_dim, output_dim)
    for i in 1:lastindex(pp.pattern[:,1])
        for j in 1:lastindex(pp.pattern[:,1])            
            dist_matrix[i,j] = sqrt((pp.pattern[i,1]-pp.pattern[j,1])^2 + (pp.pattern[i,2]-pp.pattern[j,2])^2)
        end
    end
    return dist_matrix
end

function distance_calc(x::Array)
    output_dim = length(x[:,1]) #finds size of output Array
    dist_matrix = Array{Float64}(undef, output_dim, output_dim)
    for i in 1:lastindex(x[:,1])
        for j in 1:lastindex(x[:,1])            
            dist_matrix[i,j] = sqrt((x[i,1]-x[j,1])^2 + (x[i,2]-x[j,2])^2)
        end
    end
    return dist_matrix
end

function RipleyK(pp::CartesianPointPattern, r::Union{Float64,Int}, rmax::Union{Float64,Int}, breaks = Union{Any, Nothing}, correction::Union{Array{String},String} = ["border", "isotropic", "Ripley", "translate"], nlarge::Union{Int,Nothing} = 3000, domain = Union{Any, Nothing}, return_var::Bool=false, return_ratio::Bool=false )
    npts = length(pp)
    win = pp.window
    K_area = ppWin_area(win)
    λ = npts/K_area
    λ2 = (npts * (npts-1))/(area^2)


