#=
Julia program to simulate spreading of COVID-19
Julia translation of https://github.com/MuAuan/collective_particles
=#

"""
Module for global parameters and variables
"""
module ParamVar
    struct Parameters
        num_particles::Int64  # Number of particles
        max_iteration::Int64  # Total number of time iteration
        #
        x_range::Float64  # Size of computational domain [0, x_range]
        y_range::Float64
        #
        vel_std::Float64  # Range of particle moving speed
        #
        recovery_time::Int64  # Recovery time
        infection_chance::Float64  # Chance of infection per unit time
        #
        radius_infection::Float64  # Radius of infection-zone
    end

    mutable struct Variables
        num_g::Int64  # Number of never-infected particles
        num_r::Int64  # Number of infected particles
        num_o::Int64  # Number of had-infected particles
    end

    mutable struct Particle
        pos_x::Float64
        pos_y::Float64
        status::Char  # 'g':never been infected, 'r':infected, 'o':had been infected
        t_ifcn::Int64
        flag_ifcn::Bool
        # Constructor
        Particle() = new()
    end
end


"""
Module for simulation
"""
module TimeMarch
using Distributions
    """
    Compute (squared) relative distance of two particles
    """
    function compute_relative_distance(x1, y1, x2, y2)
        return (x1-x2)^2 + (y1-y2)^2
    end


    """
    Compute new coordinates of particles ensuring periodic boundary condition
    point ∈ [0, range]
    """
    function ensure_periodic_bc(point, range)
        if point > range
            new_point = point - range
        elseif point < 0.0
            new_point = point + range
        else
            new_point = point
        end

        return new_point
    end


    """
    Update particle properties
    """
    function update_particles(param, ptcl)
        x_ifcn = Array{Float64}(undef, 0)
        y_ifcn = Array{Float64}(undef, 0)

        # Update properties of infected particles
        for itr_ptcl = 1:param.num_particles
            # Recovery from infection
            s = ptcl[itr_ptcl].status
            t = ptcl[itr_ptcl].t_ifcn
            if s == 'r' && t >= param.recovery_time
                ptcl[itr_ptcl].status = 'o'  # Hold infection history
                ptcl[itr_ptcl].flag_ifcn = true  # Hold infection history
            end

            # Store position of infected particles
            if s == 'r'
                append!(x_ifcn, ptcl[itr_ptcl].pos_x)
                append!(y_ifcn, ptcl[itr_ptcl].pos_y)
                ptcl[itr_ptcl].t_ifcn += 1
            end
        end

        # Update properties of never-infected particles
        for itr_ptcl = 1:param.num_particles
            x = ptcl[itr_ptcl].pos_x
            y = ptcl[itr_ptcl].pos_y
            s = ptcl[itr_ptcl].status
            t = ptcl[itr_ptcl].t_ifcn
            f = ptcl[itr_ptcl].flag_ifcn

            # Loop of infected particles
            for itr_ifcn = 1:length(x_ifcn)
                r2 = compute_relative_distance(x, y, x_ifcn[itr_ifcn], y_ifcn[itr_ifcn])
                if param.radius_infection^2 <= r2 && rand(Uniform(0.0, 1.0)) <= param.infection_chance
                    if s == 'g'  # If the particle has never been infected, get infected
                        s = 'r'
                        t = 0
                        f = true
                    end
                end
            end

            ptcl[itr_ptcl].pos_x = x
            ptcl[itr_ptcl].pos_y = y
            ptcl[itr_ptcl].status = s
            ptcl[itr_ptcl].t_ifcn = t
            ptcl[itr_ptcl].flag_ifcn = f
        end

        # Update position of all particles
        for itr_ptcl = 1:param.num_particles
            x = ptcl[itr_ptcl].pos_x
            y = ptcl[itr_ptcl].pos_y

            vx = param.vel_std * rand(Uniform(-1.0, 1.0))
            vy = param.vel_std * rand(Uniform(-1.0, 1.0))

            x_new = x + vx
            y_new = y + vy

            # Ensure periodic boundary condition
            x_new = ensure_periodic_bc(x_new, param.x_range)
            y_new = ensure_periodic_bc(y_new, param.y_range)

            ptcl[itr_ptcl].pos_x = x_new
            ptcl[itr_ptcl].pos_y = y_new
        end
    end


    """
    Count number of
    - Number of never-infected particles
    - Number of infected particles
    - Number of had-infected particles
    """
    function count_gro(param, ptcl)
        num_g, num_r, num_o = 0, 0, 0

        for itr_ptcl = 1:param.num_particles
            # Count number of particles according to their status
            if ptcl[itr_ptcl].status == 'g'  # never-infected
                num_g += 1
            elseif ptcl[itr_ptcl].status == 'r'  # infected
                num_r += 1
            elseif ptcl[itr_ptcl].status == 'o'  # had-infected
                num_o += 1
            else
                throw(DomainError(ptcl[itr_ptcl].status, "status must be 'g', 'r' or 'o'"))
            end
        end

        return num_g, num_r, num_o
    end
end


"""
Module for plot
"""
module Output
using Plots
gr()
font = Plots.font("Times New Roman", 20)
    """
    Scatter plot of particles
    """
    function plot_particles(itr, param, ptcl)
        x_list = Array{Float64}(undef, param.num_particles)
        y_list = Array{Float64}(undef, param.num_particles)

        # Extract necessary information
        for itr_ptcl = 1:param.num_particles
            x_list[itr_ptcl] = ptcl[itr_ptcl].pos_x
            y_list[itr_ptcl] = ptcl[itr_ptcl].pos_y
        end

        itr_str = lpad(itr, 4, "0")
        filename = string("fig/itr_", itr_str, ".png")
        p = scatter(
            x_list, y_list,
            aspect_ratio = 1,
            xlims = (0.0, param.x_range),
            ylims = (0.0, param.y_range),
            xaxis = nothing,
            yaxis = nothing,
            title = string("itr = ", itr))
        savefig(p, filename)
    end
end


# ========================================
# Main function
# ========================================

## Declare modules
using Distributions
using Printf
using .ParamVar
using .TimeMarch:
update_particles,
count_gro
using .Output

# ----------------------------------------
## Set parameters & variables
# ----------------------------------------
num_particles = 10  # 100
max_iteration = 10  # 100

x_range = 10.0
y_range = 10.0

vel_std = 1.0

recovery_time = 3  # 30
infection_chance = 1.0  # 0.03

radius_infection = 5

### Declare parameters
param = ParamVar.Parameters(
    num_particles,max_iteration,
    x_range,y_range,
    vel_std,
    recovery_time,infection_chance,
    radius_infection)

num_g = 0
num_r = 0
num_o = 0

### Declare parameters
var = ParamVar.Variables(
    num_g,num_r,num_o)

### Define array of particle properties
particles = Array{ParamVar.Particle}(undef, param.num_particles)
for itr_ptcl = 1:param.num_particles
    particles[itr_ptcl] = ParamVar.Particle()
end
# particles .= ParamVar.Particle()  # ERROR: LoadError: MethodError: no method matching length(::Main.ParamVar.Particle)


# ----------------------------------------
## Set initial condition of particles
# ----------------------------------------
for itr_ptcl = 1:param.num_particles
    particles[itr_ptcl].pos_x = rand(Uniform(0.0, param.x_range))
    particles[itr_ptcl].pos_y = rand(Uniform(0.0, param.y_range))
    particles[itr_ptcl].status = 'g'  # Initially not infected
    if itr_ptcl == 1  # One particle is initially infected
        particles[itr_ptcl].status = 'r'
    end
    particles[itr_ptcl].t_ifcn = 0
    particles[itr_ptcl].flag_ifcn = false
end
# particles.pos_x .= rand(Uniform(0.0, param.x_range))  # ERROR: LoadError: type Array has no field pos_x

#=
println("x = ", getfield.(particles, :pos_x))
println("y = ", getfield.(particles, :pos_y))
println("status = ", getfield.(particles, :status))
println("t_ifcn = ", getfield.(particles, :t_ifcn))
println("flag_ifcn = ", getfield.(particles, :flag_ifcn))
=#

# ----------------------------------------
## Time iteration of infection simulation
# ----------------------------------------
for itr_time = 1:param.max_iteration
    update_particles(param, particles)
    var.num_g, var.num_r, var.num_o = count_gro(param, particles)
    Output.plot_particles(itr_time, param, particles)
    # tmp_string = @sprintf "itr_time = %i x[1] = %6.3f y[1] = %6.3f" itr_time particles[1].pos_x particles[1].pos_y
    # println(tmp_string)
    println("itr_time = ", itr_time, " status = ", getfield.(particles, :status))
end
