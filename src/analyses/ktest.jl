using CSV
using DataFrames

pointData = CSV.read("Data/Stem_Map.csv", DataFrame)
#### Testing Data
my_array = Array{Union{Nothing, Int}}(nothing, 2, 5)
my_pp=CartesianPointPattern((-100,100),my_array)
abstract type PointPattern end

struct CartesianPointPattern <: PointPattern
    #Use with basic x,y coordinate data
    window::Tuple{Float32, Float32}
    pattern::Array
    marks::Union{Vector,Nothing}
end

struct GeographicPointPattern <: PointPattern
    #Possible interfaces with geographic projections


myPPP = pointpattern((-150,150), pointData[:,(2:3)], pointData[:,(4)])

function AsPointPattern(data::Any, win::Any)
    try
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
    pp_win = convert(win,Tuple{Float32,Float32})
    pp_out = CartesianPointPattern(pp_win, pp_pattern, pp_marx)
    return pp_out
end