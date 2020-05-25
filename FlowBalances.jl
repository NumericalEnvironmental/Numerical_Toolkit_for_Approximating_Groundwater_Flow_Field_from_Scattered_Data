###################################################################################################
#
# FlowBalances.jl
#
# by Walt McNab (May 2020)
#
# (1) Approximate potentiometric surface and hydraulic conductivity distributions by interpolation with kriging
# (2) Estimate groundwater flow vector field by Darcy's law; optional particle tracking
# (3) Numerical integration across polygon drawn on flow field to compute fluxes
#
# Current version assumes datum for head values is equal to an elevation of zero.
#
###################################################################################################


using DelimitedFiles
using Statistics
using Polynomials
using QuadGK


# miscellaneous fixed parameters
const epsilon = 0.01                	# factor used to prevent binning problems with variogram construction
const delta = 0.1 						# distance offset used to estimate hydraulic gradient


### defined types ###


mutable struct Polygon
    x::Array{Float64, 1}                # apex locations
    y::Array{Float64, 1}
    m::Array{Float64, 1}                # slope of polygon side
    n::Array{Any, 1}                    # normal vector to polygon side
    L::Array{Float64, 1}                # length of polygon side    
end


mutable struct TrackParams
	dStep::Float64 						# step size for particle tracking
	phi::Float64 						# porosity
end


mutable struct Sources              	# source/sink terms (i.e., wells)
	activeSources::Bool
    x::Array{Float64, 1}
    y::Array{Float64, 1}
    Q::Array{Float64, 1}    
end


mutable struct Particle
	x::Array{Float64, 1} 				# particle locations, starting with initial
	y::Array{Float64, 1}
	vx::Array{Float64, 1} 				# particle velocities (for plotting)
	vy::Array{Float64, 1}	
	t::Array{Float64, 1}	
	fwd::Float64 						# movement direction: 1 = forward, -1 = backward
	active::Bool
end


mutable struct Krig 					# variogram parameters
	varMatrix::Array{Float64, 2}
	m::Float64
	b::Float64
end


mutable struct Points
    x::Array{Float64, 1}
    y::Array{Float64, 1} 
    v::Array{Float64, 1} 				# head measurements, or log hydraulic conductivity
end


mutable struct Grid
	x0::Float64
	xf::Float64
	y0::Float64
	yf::Float64
    x::Array{Float64, 1}
    y::Array{Float64, 1}
end


### utility functions ###


function DistVector(xp::Float64, yp::Float64, xArray::Array{Float64, 1}, yArray::Array{Float64, 1})::Array{Float64, 1}
	# return a vector of distances
	d = Float64[]
    for i = 1:length(xArray)
		dist = sqrt((xp-xArray[i])^2 + (yp-yArray[i])^2)
        push!(d, dist)
    end	
	return d
end


function dot(a::Array{Float64, 1}, b::Array{Float64, 1})::Float64
    # vector dot product
    return sum(a.*b)
end


### hydro modeling functions ###


function Grad(x::Float64, y::Float64, headData::Points, krigHead::Krig)
	# estimate hydraulic gradient by finite differences
	h = Interpolate(x, y, headData, krigHead)
	hx = Interpolate(x+delta, y, headData, krigHead)
	hy = Interpolate(x, y+delta, headData, krigHead)
	dhdx = (hx-h)/delta
	dhdy = (hy-h)/delta
	return dhdx, dhdy
end


function Darcy(x::Float64, y::Float64, headData::Points, krigHead::Krig, KData::Points, krigLogK::Krig)
	# local groundwater Darcy velocity
	K = 10.0^Interpolate(x, y, KData, krigLogK)
	dhdx, dhdy = Grad(x, y, headData, krigHead)
	qx = -K * dhdx
	qy = -K * dhdy
	return qx, qy
end


function Displace(x0::Float64, y0::Float64, vx::Float64, vy::Float64, d::Float64, fwd::Float64)
	# particle movement
	vm = sqrt(vx^2 + vy^2) 	# velocity magnitude
	dt = d/vm 				# time step to move "d"
	x = x0 + fwd*vx*dt 		# implied new position
	y = y0 + fwd*vy*dt
	return x, y, dt
end


function CheckEndState(part::Particle, trackParams::TrackParams, grid::Grid, wells::Sources)::Particle
	# note if particle has exited model domain or has approached too close to a well
	part.active = (part.x[end]>=grid.x0) & (part.x[end]<=grid.xf) & (part.y[end]>=grid.y0) & (part.y[end]<=grid.yf)
	if wells.activeSources == true
		dWells = DistVector(part.x[end], part.y[end], wells.x, wells.y)
		if minimum(dWells) <= trackParams.dStep
			part.active = false
		end
	end
	return part
end


function Track(part::Particle, trackParams::TrackParams, headData::Points, krigHead::Krig, KData::Points, krigLogK::Krig, wells::Sources, grid::Grid)::Particle
	# follow the progress of a particle (initial conditions already defined)
	while part.active == true
		# single-pass particle position estimator
		xf, yf, dt = Displace(part.x[end], part.y[end], part.vx[end], part.vy[end], trackParams.dStep, part.fwd)
		vxf, vyf = Darcy(xf, yf, headData, krigHead, KData, krigLogK)./trackParams.phi
		vx = 0.5*part.vx[end] + 0.5*vxf
		vy = 0.5*part.vy[end] + 0.5*vyf
		x, y, dt = Displace(part.x[end], part.y[end], vx, vy, trackParams.dStep, part.fwd)
		push!(part.x, x)
		push!(part.y, y)		
		push!(part.vx, vx)
		push!(part.vy, vy)	
		push!(part.t, part.t[end]+dt)
		part = CheckEndState(part, trackParams, grid, wells)        # check to see if exit condition is met
	end
	return part
end


function FluxNorm(x::Float64, headData::Points, krigHead::Krig, KData::Points, krigLogK::Krig, polygon::Polygon, indx::Int64)
    # vertically-integrated, normal groundwater flux at a point along a transect
    y = polygon.y[indx] + polygon.m[indx]*(x-polygon.x[indx])
    qx, qy = Darcy(x, y, headData, krigHead, KData, krigLogK).*Interpolate(x, y, headData, krigHead) 
    return -dot([qx, qy], polygon.n[indx]) 
end
        

function SumPolyFluxes(polygon::Polygon, headData::Points, krigHead::Krig, KData::Points, krigLogK::Krig)::Float64
    # integrate groundwater fluxes across polygon/polyline sides
    QTot = 0.0
    for i = 1:length(polygon.x)-1
        x0 = minimum([polygon.x[i], polygon.x[i+1]])
        xf = maximum([polygon.x[i], polygon.x[i+1]]) 
        QTot += (polygon.L[i]/abs(xf-x0)) * quadgk(x -> FluxNorm(x, headData, krigHead, KData, krigLogK, polygon, i), x0, xf)[1]
    end
    return QTot
end


### kriging functions ###


function MakeKrig(data, setName)
	# parameters to inform kriging
    S = Squareform(data.x, data.y)      			# distance matrix	
    m, b = Variogram(data, S, setName)       		# linear variogram	
    varMatrix = (S.*m).+b	
	return Krig(varMatrix, m, b)
end


function Squareform(xArray::Array{Float64, 1}, yArray::Array{Float64, 1})::Array{Float64, 2}
	# squareform distance matrix
	N = length(xArray)
	S = zeros(Float64, N, N)
    for i = 1:N
		S[i, :] = DistVector(xArray[i], yArray[i], xArray, yArray)
    end		
	return S
end


function Variogram(points::Points, S::Array{Float64, 2}, setName::AbstractString)
    # create variogram; write variogram summary to file and return linear variogram model
    dist = Float64[]
    gamma = Float64[]
    lagGroup = []
    lagMeanDist = Float64[]
    xVar = Float64[]                      
    yVar = Float64[]
    N = length(points.x)        # number of lag bins
    # difference in value versus distance offset
    for i = 1:N
        for j = i+1:N 
            push!(dist, S[i,j])
            push!(gamma, (points.v[i]-points.v[j])^2)
        end
    end    
    # populate bins
    lagStepSize = maximum(dist)/N
    for i = 1:N
        push!(lagGroup, [])                             # set up list of lag sets
        push!(lagMeanDist, (i-0.5)*lagStepSize)         # mean distance offset for each of the lag groups
    end
    # assign lags to groups
    for i = 1:length(dist)
        lagIndex = floor(Int, (dist[i]-epsilon)/lagStepSize) + 1
        push!(lagGroup[lagIndex], gamma[i])    
    end
    # calculate variogram
    for i = 1:N
        if length(lagGroup[i]) > 1     # binned lag group is large enough to calculate a variance
            push!(xVar, lagMeanDist[i])
            push!(yVar, mean(lagGroup[i])/2)    
        end
    end
    WriteVariogram(xVar, yVar, setName)					# write to output file
    # linear regression (i.e. linear variogram model)
    r = fit(xVar, yVar, 1)
    return r[1], r[0]               # = return m, b
end


function Interpolate(x::Float64, y::Float64, points::Points, krig::Krig)::Float64
    # interpolate field value at (x, y)
    dist = DistVector(x, y, points.x, points.y)          # calculate covariance vector for distances between (x,y) and data points
    k = (dist.*krig.m).+krig.b
    wt = \(krig.varMatrix, k)                            # compute weights
    z = sum(wt.*points.v)
    return z
end


### I/O functions ###


function ReadParticles(trackParams::TrackParams, headData::Points, krigHead::Krig, KData::Points, krigLogK::Krig)::Array{Particle, 1}
	# initialize particles
	part = []
	data = readdlm("particleStarts.txt", '\t', header=true)
    for i = 1:size(data[1], 1)
        x = Float64(data[1][i, 1])
        y = Float64(data[1][i, 2])
        fwd = Float64(data[1][i, 3])		
		vx, vy = Darcy(x, y, headData, krigHead, KData, krigLogK)./trackParams.phi
		t = 0.0
		active = true
        push!(part, Particle([x], [y], [vx], [vy], [t], fwd, active))
    end
    println("Read particle starting positions and initialized.")
	return part
end


function ReadPolygon()::Polygon
    # clockwise list convention; repeat 1st point as last point
    x = Float64[]
    y = Float64[]
    m = Float64[]    
    n = []
    L = Float64[]
    data = readdlm("polygon.txt", '\t', header=true)
    for i = 1:size(data[1], 1)
        push!(x, Float64(data[1][i, 1]))
        push!(y, Float64(data[1][i, 2]))   
    end
    for i = 1:length(x)-1
        dx = x[i+1] - x[i]            
        dy = y[i+1] - y[i]
        push!(m, dy/dx)
        d = sqrt(dx^2 + dy^2)
        push!(L, d)
        push!(n, [-dy/d, dx/d])
    end
    println("Read and processed bounding polygon.")
    return Polygon(x, y, m, n, L)
end


function ReadPoints(dataSet::AbstractString)::Points
    # read point data set and return array of point structs
    x = Float64[]
    y = Float64[]
    v = Float64[]	
    data = readdlm(dataSet * ".txt", '\t', header=true)
    for i = 1:size(data[1], 1)
        push!(x, Float64(data[1][i, 1]))
        push!(y, Float64(data[1][i, 2]))
        push!(v, Float64(data[1][i, 3]))
    end
    println("Read " * dataSet * ".")
    return Points(x, y, v)
end


function ReadGrid()
    # set up grid for kriging and inverse model
    x = Float64[]
    y = Float64[]
    data = readdlm("grid.txt", '\t', header=true)
    x0 = Float64(data[1][1, 2])
    y0 = Float64(data[1][1, 3])
    xf = Float64(data[1][2, 2])
    yf = Float64(data[1][2, 3])
    nx = Int64(data[1][3, 2])
    ny = Int64(data[1][3, 3])    
    dx = (xf-x0)/nx
    dy = (yf-y0)/ny
    for j = 1:ny
        for i = 1:nx
            push!(x, x0 + (i-0.5)*dx)
            push!(y, y0 + (j-0.5)*dy)
        end
    end
    println("Read grid constraints.")
    return Grid(x0, xf, y0, yf, x, y)
end


function ReadSources(dataSet::AbstractString)::Sources
    # read source/sink terms
    x = Float64[]
    y = Float64[]
    Q = Float64[]
	activeSources = true
    data = readdlm(dataSet * ".txt", '\t', header=true)
    for i = 1:size(data[1], 1)
        push!(x, Float64(data[1][i, 1]))
        push!(y, Float64(data[1][i, 2]))
        push!(Q, Float64(data[1][i, 3]))
    end
    println("Read " * dataSet * ".")
    return Sources(activeSources, x, y, Q)
end


function ReadTrack()
	# miscellaneous particle tracking parameters
    data = readdlm("trackParams.txt", '\t', header=false)
    dStep = Float64(data[1, 2])
	phi = Float64(data[2, 2])
    println("Read tracking parameters.")	
	return TrackParams(dStep, phi)
end


function WriteVariogram(xVar::Array{Float64, 1}, yVar::Array{Float64, 1}, setName::AbstractString)
	# variogram output file
    csvfile = open("variogram_" * setName * ".csv","w")
    line_out = "h" * "," * "gamma"
    println(csvfile, line_out)    
    for i = 1:length(xVar)
        line_out = string(xVar[i]) * "," * string(yVar[i])
        println(csvfile, line_out)    
    end
    close(csvfile)
end


function WriteGridOut(grid::Grid, head::Array{Float64, 1}, logK::Array{Float64, 1})
	# write grid output, with interpolated heads and log K estimates, to file
    csvfile = open("gridResults.csv","w")
    line_out = "x" * "," * "y" * "," * "head" * "," * "logK"
    println(csvfile, line_out)    
    for i = 1:length(grid.x)
        line_out = string(grid.x[i]) * "," * string(grid.y[i]) * "," * string(head[i]) * "," * string(logK[i])
        println(csvfile, line_out)    
    end
    close(csvfile)
	println("Wrote grid output file.")
end


function WriteParticles(part::Array{Particle, 1})
	# write grid output, with interpolated heads and log K estimates, to file
    csvfile = open("particleTracks.csv","w")
    line_out = "particle" * "," * "time" * "," * "x" * "," * "y" * "," * "vx" * "," * "vy" * "," * "v-magn"
    println(csvfile, line_out)    
    for (i, p) in enumerate(part)
		for j = 1:length(p.x)
			vMagn = sqrt(p.vx[j]^2 + p.vy[j]^2)
			line_out = string(i) * "," * string(p.t[j]) * "," * string(p.x[j]) * "," * string(p.y[j]) * "," * string(p.vx[j]) * "," * string(p.vy[j]) * "," * string(vMagn)
			println(csvfile, line_out)
		end
    end
    close(csvfile)
	println("Wrote particle tracks output file.")
end


### main ###


function FlowBalances(gridFlag, trackFlag, polyFlag, wellsFlag)

    # problem setup
	dataSet = ["heads", "logKs"]
    grid = ReadGrid()                       	# read grid
	headData = ReadPoints(dataSet[1])          	# head measurements
	KData = ReadPoints(dataSet[2])          	# log K measurements
	if wellsFlag == true 						# well locations and pumping rates, if applicable
		wells = ReadSources("wells")
	else
		wells = Sources(false, Float64[], Float64[], Float64[]) 		# empty struct, with activeSources = false
	end

	# set up kriging models
	krigHead = MakeKrig(headData, dataSet[1])
	krigLogK = MakeKrig(KData, dataSet[2])
	
	# read polygon, if used
    if polyFlag == true
        polygon = ReadPolygon()                                                     # read input file
        QIntegrate = SumPolyFluxes(polygon, headData, krigHead, KData, krigLogK)    # integrate
        println("Flux across polygon/polyline = ", "\t", QIntegrate)
    end
	
	# conduct particle tracking, if implemented
	if trackFlag == true
		trackParams = ReadTrack()													# read parameters file
		part = ReadParticles(trackParams, headData, krigHead, KData, krigLogK)		# initialize particles
		for i = 1:length(part)														# track particles
			if part[i].active == true
				part[i] = Track(part[i], trackParams, headData, krigHead, KData, krigLogK, wells, grid)
			end
		end
		WriteParticles(part)							# write to output file
	end
	
	# write to output grid for plotting
    if gridFlag == true
        head = Float64[]
        logK = Float64[]
        println("Interpolating heads across grid ...")	
        for i = 1:length(grid.x)
            push!(head, Interpolate(grid.x[i], grid.y[i], headData, krigHead))
        end
        println("Interpolating log K estimates across grid ...")	
        for i = 1:length(grid.x)
            push!(logK, Interpolate(grid.x[i], grid.y[i], KData, krigLogK))
        end	
        WriteGridOut(grid, head, logK)        # write to output file
	end
    
    println("Finished.")

end

### run model ###

gridFlag = false
trackFlag = false
polyFlag = true
wellsFlag = true
FlowBalances(gridFlag, trackFlag, polyFlag, wellsFlag)          


