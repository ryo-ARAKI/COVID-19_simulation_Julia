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
        vel_mean::Float64  # Average of velocity component
        vel_σ::Float64  # Standard deviation
        vel_flc::Float64  # Range of fluctuation
        #
        ratio_infection_init::Float64  # Initial percentage of infected particles
        recovery_time::Int64  # Recovery time
        infection_chance::Float64  # Chance of infection per unit time
        #
        radius_infection::Float64  # Radius of infection-zone
    end

    mutable struct Variables
        num_not_infected::Int64  # Number of never-infected particles
        num_infected::Int64  # Number of infected particles
        num_recovered::Int64  # Number of recovered particles
    end

    mutable struct Particle
        pos_x::Float64  # Position x-component
        pos_y::Float64
        vel_x::Float64  # Velocity x-component
        vel_y::Float64
        status::Char  # 'g':never been infected, 'r':infected, 'o':had been infected
        t_ifcn::Int64  # Infected time
        flag_ifcn::Bool  # Infection
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
    Compute initial condition
    """
    function set_initial_condition(param, ptcl)
        num_infected_init = Int64(param.ratio_infection_init * param.num_particles)
        for itr_ptcl = 1:param.num_particles
            ptcl[itr_ptcl].pos_x = rand(Uniform(0.0, param.x_range))  # Uniform distribution
            ptcl[itr_ptcl].pos_y = rand(Uniform(0.0, param.y_range))
            ptcl[itr_ptcl].vel_x = rand(Normal(param.vel_mean, param.vel_σ))  # Gaussian distribution
            ptcl[itr_ptcl].vel_y = rand(Normal(param.vel_mean, param.vel_σ))
            ptcl[itr_ptcl].status = 'g'  # Initially not infected
            if itr_ptcl <= num_infected_init  # Some particles are initially infected
                ptcl[itr_ptcl].status = 'r'
            end
            ptcl[itr_ptcl].t_ifcn = 0
            ptcl[itr_ptcl].flag_ifcn = false
        end
        # ptcl.pos_x .= rand(Uniform(0.0, param.x_range))  # ERROR: LoadError: type Array has no field pos_x
    end


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
                if r2 < param.radius_infection^2 && rand(Uniform(0.0, 1.0)) <= param.infection_chance
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
            vx = ptcl[itr_ptcl].vel_x
            vy = ptcl[itr_ptcl].vel_y

            # Update velocity
            vx_flc = param.vel_flc * rand(Uniform(-1.0, 1.0))
            vy_flc = param.vel_flc * rand(Uniform(-1.0, 1.0))

            vx_new = vx + vx_flc
            vy_new = vy + vy_flc

            # Update position
            x_new = x + vx
            y_new = y + vy

            # Ensure periodic boundary condition
            x_new = ensure_periodic_bc(x_new, param.x_range)
            y_new = ensure_periodic_bc(y_new, param.y_range)

            ptcl[itr_ptcl].pos_x = x_new
            ptcl[itr_ptcl].pos_y = y_new
            ptcl[itr_ptcl].vel_x = vx_new
            ptcl[itr_ptcl].vel_y = vy_new
        end
    end


    """
    Count number of
    - Number of never-infected particles
    - Number of infected particles
    - Number of recovered particles
    """
    function count_status(param, ptcl)
        num_not_infected, num_infected, num_recovered = 0, 0, 0

        for itr_ptcl = 1:param.num_particles
            # Count number of particles according to their status
            if ptcl[itr_ptcl].status == 'g'  # never-infected
                num_not_infected += 1
            elseif ptcl[itr_ptcl].status == 'r'  # infected
                num_infected += 1
            elseif ptcl[itr_ptcl].status == 'o'  # recovered
                num_recovered += 1
            else
                throw(DomainError(ptcl[itr_ptcl].status, "status must be 'g', 'r' or 'o'"))
            end
        end

        return num_not_infected, num_infected, num_recovered
    end
end


"""
Module for plot
"""
module Output
    using Plots
    font = Plots.font("Times New Roman", 20)
    """
    Scatter plot of particles
    """
    function plot_particles(itr, param, var, ptcl)
        x_g = Array{Float64}(undef, var.num_not_infected)
        y_g = Array{Float64}(undef, var.num_not_infected)
        x_r = Array{Float64}(undef, var.num_infected)
        y_r = Array{Float64}(undef, var.num_infected)
        x_o = Array{Float64}(undef, var.num_recovered)
        y_o = Array{Float64}(undef, var.num_recovered)

        # Extract necessary information
        count_g, count_r, count_o = 1, 1, 1
        for itr_ptcl = 1:param.num_particles
            if ptcl[itr_ptcl].status == 'g'  # never-infected
                x_g[count_g] = ptcl[itr_ptcl].pos_x
                y_g[count_g] = ptcl[itr_ptcl].pos_y
                count_g += 1
            elseif ptcl[itr_ptcl].status == 'r'  # infected
                x_r[count_r] = ptcl[itr_ptcl].pos_x
                y_r[count_r] = ptcl[itr_ptcl].pos_y
                count_r += 1
            elseif ptcl[itr_ptcl].status == 'o'  # recovered
                x_o[count_o] = ptcl[itr_ptcl].pos_x
                y_o[count_o] = ptcl[itr_ptcl].pos_y
                count_o += 1
            end
        end

        itr_str = lpad(itr, 4, "0")
        filename = string("fig/itr_", itr_str, ".png")
        p = scatter(  # never-infected
            x_g, y_g,
            markercolor = :green,
            label = "Never infected")
        p! = scatter!(  # infected
            x_r, y_r,
            markercolor = :red,
            label = "Infected")
        p! = scatter!(  # recovered
            x_o, y_o,
            aspect_ratio = 1,
            markercolor = :orange,
            label = "recovered",
            xlims = (0.0, param.x_range),
            ylims = (0.0, param.y_range),
            axis = nothing,
            size=(640, 640),
            title = string("itr = ", itr))
        plot(p)
        # savefig(p, filename)
    end

    """
    Make gif video
    """
    function make_gif(param,anim)
        gif(
            anim,
            "fig/particles.gif",
            fps=5)
    end
end


# ========================================
# Main function
# ========================================

## Declare modules
using ProgressMeter
using Distributions
using Printf
using Plots
gr(
    markerstrokewidth = 0,
    markersize = 10
)
using .ParamVar
using .TimeMarch:
set_initial_condition,
update_particles,
count_status
using .Output:
plot_particles,
make_gif

# ----------------------------------------
## Set parameters & variables
# ----------------------------------------
num_particles = 100
max_iteration = 100

x_range = 1.0
y_range = 1.0

vel_mean = 0.01
vel_σ = 0.01
vel_flc = 0.3 * vel_mean

ratio_infection_init = 0.03
recovery_time = 5
infection_chance = 0.3

radius_infection = 0.1

### Declare parameters
param = ParamVar.Parameters(
    num_particles,max_iteration,
    x_range,y_range,
    vel_mean,vel_σ,vel_flc,
    ratio_infection_init,recovery_time,infection_chance,
    radius_infection)

num_not_infected, num_infected, num_recovered = 0, 0, 0

### Declare variables
var = ParamVar.Variables(
    num_not_infected,num_infected,num_recovered)

### Define array of particle properties
particles = Array{ParamVar.Particle}(undef, param.num_particles)
for itr_ptcl = 1:param.num_particles
    particles[itr_ptcl] = ParamVar.Particle()
end
# particles .= ParamVar.Particle()  # ERROR: LoadError: MethodError: no method matching length(::Main.ParamVar.Particle)


# ----------------------------------------
## Set initial condition of particles
# ----------------------------------------
set_initial_condition(param, particles)


# ----------------------------------------
## Time iteration of infection simulation
# ----------------------------------------
progress = Progress(param.max_iteration)
anim = @animate for itr_time = 1:param.max_iteration
    # Update particle properties
    update_particles(param, particles)

    # Count not-infected, infected & recovered number of particles
    var.num_not_infected, var.num_infected, var.num_recovered = count_status(param, particles)

    # Plot particles for gif video
    plot_particles(itr_time, param, var, particles)

    # tmp_string = @sprintf "itr_time = %i x[1] = %6.3f y[1] = %6.3f" itr_time particles[1].pos_x particles[1].pos_y
    # println(tmp_string)
    # println("itr_time = ", itr_time, " g = ", var.num_not_infected, " r = ", var.num_infected, " o = ", var.num_recovered)
    next!(progress)

    # Finish if there are no infected particles any more
    if var.num_infected == 0
        println("\n No patients at itr_time = ", itr_time, " :Exit")
        break
    end
end


# ----------------------------------------
## Make gif video of simulation
# ----------------------------------------
make_gif(param, anim)
