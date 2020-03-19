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

    struct Flags
        flag_multiple_infection::Bool  # Multiple infection
        flag_infected_isolation::Bool  # Isolation (stop) of infected particle
    end

    mutable struct Particle
        pos_x::Float64  # Position x-component
        pos_y::Float64
        vel_x::Float64  # Velocity x-component
        vel_y::Float64
        status::AbstractString  # "not_infected" "infected" or "recovered"
        t_ifcn::Int64  # Infected time
        past_ifcn::Int64  # Past infection history
        # Constructor
        Particle() = new()
    end

    mutable struct NumSnapshot
        not_infected::Int64  # Number of never-infected particles
        infected::Int64  # Number of infected particles
        recovered::Int64  # Number of recovered particles
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
            ptcl[itr_ptcl].t_ifcn = 0
            ptcl[itr_ptcl].past_ifcn = 0
            ptcl[itr_ptcl].status = "not_infected"  # Initially not infected
            # Some particles are initially infected
            if itr_ptcl <= num_infected_init
                ptcl[itr_ptcl].status = "infected"
                ptcl[itr_ptcl].past_ifcn = 1
            end
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
    Check whether the particle will be infected or not
    """
    function compute_is_infected(
        infection_chance, flag_multiple_infection,
        status, past_ifcn
    )
        isinfected = false

        # Infection chance & past infection criteria
        if rand(Uniform(0.0, 1.0)) <= infection_chance && status == "not_infected"
            isinfected = true
        end

        if flag_multiple_infection  # Multiple infection is allowed
            if status == "recovered"  # Infected in the past
                # Infection chance criteria is modified by past number of infection
                if rand(Uniform(0.0, 1.0)) <= infection_chance * (0.5^past_ifcn)
                    isinfected = true
                end
            end
        end

        return isinfected
    end


    """
    Compute updated velocity components
    """
    function compute_new_velocity(
        flag_infected_isolation, vel_flc, status,
        vx, vy
    )

        vx_flc = vel_flc * rand(Uniform(-1.0, 1.0))
        vy_flc = vel_flc * rand(Uniform(-1.0, 1.0))

        vx_new = vx + vx_flc
        vy_new = vy + vy_flc

        if flag_infected_isolation && status == "infected"
            vx_new, vy_new = 0.0, 0.0
        end

        return vx_new, vy_new
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
    function update_particles(param, flag, ptcl)
        x_ifcn = Array{Float64}(undef, 0)
        y_ifcn = Array{Float64}(undef, 0)

        # Update properties of infected particles
        for itr_ptcl = 1:param.num_particles
            # Recovery from infection
            s = ptcl[itr_ptcl].status
            t = ptcl[itr_ptcl].t_ifcn
            if s == "infected" && t >= param.recovery_time
                ptcl[itr_ptcl].status = "recovered"  # Hold infection history
                ptcl[itr_ptcl].past_ifcn += 1  # Hold infection history
            end

            # Store position of infected particles
            if s == "infected"
                append!(x_ifcn, ptcl[itr_ptcl].pos_x)
                append!(y_ifcn, ptcl[itr_ptcl].pos_y)
                ptcl[itr_ptcl].t_ifcn += 1
            end
        end

        # Update properties of not_infected and recovered particles
        for itr_ptcl = 1:param.num_particles
            x = ptcl[itr_ptcl].pos_x
            y = ptcl[itr_ptcl].pos_y
            s = ptcl[itr_ptcl].status
            t = ptcl[itr_ptcl].t_ifcn
            p = ptcl[itr_ptcl].past_ifcn

            # Loop of infected particles
            for itr_ifcn = 1:length(x_ifcn)
                r2 = compute_relative_distance(x, y, x_ifcn[itr_ifcn], y_ifcn[itr_ifcn])
                if r2 < param.radius_infection^2  # Relative distance criteria
                    isinfected = compute_is_infected(
                        param.infection_chance, flag.flag_multiple_infection,
                        s, p)
                    if isinfected
                        s = "infected"
                        t = 0  # Time since infection
                        p += 1  # Past number of infection
                    end
                end
            end

            ptcl[itr_ptcl].status = s
            ptcl[itr_ptcl].t_ifcn = t
            ptcl[itr_ptcl].past_ifcn = p
        end

        # Update position and velocity of all particles
        for itr_ptcl = 1:param.num_particles
            x = ptcl[itr_ptcl].pos_x
            y = ptcl[itr_ptcl].pos_y
            vx = ptcl[itr_ptcl].vel_x
            vy = ptcl[itr_ptcl].vel_y
            s = ptcl[itr_ptcl].status

            # Update velocity
            vx_new, vy_new = compute_new_velocity(
                flag.flag_infected_isolation, param.vel_flc, s,
                vx, vy)

            # Update position
            x_new = x + vx_new
            y_new = y + vy_new

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
    - not infected particles
    - infected particles
    - recovered particles
    """
    function count_status(param, ptcl, num_snapshot, num_timeseries)
        num_not_infected, num_infected, num_recovered = 0, 0, 0

        for itr_ptcl = 1:param.num_particles
            # Count number of particles according to their status
            if ptcl[itr_ptcl].status == "not_infected"  # not infected
                num_not_infected += 1
            elseif ptcl[itr_ptcl].status == "infected"  # infected
                num_infected += 1
            elseif ptcl[itr_ptcl].status == "recovered"  # recovered
                num_recovered += 1
            else
                throw(DomainError(ptcl[itr_ptcl].status, "status must be 'not_infected', 'infected' or 'recovered'"))
            end
        end

        # Define mutable construct num_snapshot
        num_snapshot.not_infected = num_not_infected
        num_snapshot.infected = num_infected
        num_snapshot.recovered = num_recovered

        # Add timeseries data
        push!(num_timeseries, num_snapshot)
    end
end


"""
Module for plot
"""
module Output
    using Printf
    using Base.Filesystem
    using Plots
    using StatsPlots
    font = Plots.font("Times New Roman", 20)

    """
    Stdout simulation condition
    """
    function stdout_condition(param, flag)
        println("#========== Coronavirus spreading simulation ====================")
        println(@sprintf "num_particles = %i,  max_iteration = %i" param.num_particles param.max_iteration)
        println(@sprintf "(x, y) = [0.00, %.2f] × [0.00, %.2f]" param.x_range param.y_range)
        println(@sprintf "vel_mean = %.3f,  vel_σ = %.3f,  vel_flc = %.3f" param.vel_mean param.vel_σ param.vel_flc)
        println(@sprintf "ratio_infection_init = %.3f" param.ratio_infection_init)
        println(@sprintf "recovery_time = %i,  infection_chance = %.2f,  radius_infection = %.2f" param.recovery_time param.infection_chance param.radius_infection)
        println("flag_multiple_infection = ", flag.flag_multiple_infection ? "true" : "false")
        println("flag_infected_isolation = ", flag.flag_infected_isolation ? "true" : "false")
        println("#================================================================\n")

        out_dir = @sprintf "dat/np=%i_vm=%.3f_vσ=%.3f_vf=%.3f_rii=%i_rt=%i_ic=%.2f_ri=%.2f_fmi=%d_fii=%d/" param.num_particles param.vel_mean param.vel_σ param.vel_flc param.ratio_infection_init param.recovery_time param.infection_chance param.radius_infection flag.flag_multiple_infection flag.flag_infected_isolation
        mkpath(out_dir)
        return out_dir
    end


    """
    Scatter plot of particles
    """
    function plot_particles(itr, param, num_snapshot, out_dir, ptcl)
        x_not_infected = Array{Float64}(undef, num_snapshot.not_infected)
        y_not_infected = Array{Float64}(undef, num_snapshot.not_infected)
        x_infected = Array{Float64}(undef, num_snapshot.infected)
        y_infected = Array{Float64}(undef, num_snapshot.infected)
        x_recovered = Array{Float64}(undef, num_snapshot.recovered)
        y_recovered = Array{Float64}(undef, num_snapshot.recovered)

        # Extract necessary information
        count_not_infeced, count_infected, count_recovered = 1, 1, 1
        for itr_ptcl = 1:param.num_particles
            if ptcl[itr_ptcl].status == "not_infected"
                x_not_infected[count_not_infeced] = ptcl[itr_ptcl].pos_x
                y_not_infected[count_not_infeced] = ptcl[itr_ptcl].pos_y
                count_not_infeced += 1
            elseif ptcl[itr_ptcl].status == "infected"
                x_infected[count_infected] = ptcl[itr_ptcl].pos_x
                y_infected[count_infected] = ptcl[itr_ptcl].pos_y
                count_infected += 1
            elseif ptcl[itr_ptcl].status == "recovered"
                x_recovered[count_recovered] = ptcl[itr_ptcl].pos_x
                y_recovered[count_recovered] = ptcl[itr_ptcl].pos_y
                count_recovered += 1
            end
        end

        itr_str = lpad(itr, 4, "0")
        filename = string(out_dir, "itr_", itr_str, ".png")
        p = scatter(
            x_not_infected, y_not_infected,
            markercolor = :deepskyblue,
            label = "Not infected")
        p! = scatter!(
            x_infected, y_infected,
            markercolor = :orangered,
            label = "Infected")
        p! = scatter!(
            x_recovered, y_recovered,
            aspect_ratio = 1,
            markercolor = :gold,
            label = "Recovered",
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
    function make_gif(param, anim, out_dir)
        gif(
            anim,
            string(out_dir, "particles.gif"),
            fps=5)
    end

    """
    Plot number of
    - Not infected
    - Infected
    - Recovered
    particle as timeseries
    """
    function plot_num_timeseries(param, num_timeseries, out_dir)
        filename = string(out_dir, "timeseries.png")

        # Array of struct -> struct  #####POSSIBLE BETTER SOLUTION?#####
        tseries = Array{Int64}(undef, length(num_timeseries), 3)
        for itr_time = 1:length(num_timeseries)
            tseries[itr_time, 1] = num_timeseries[itr_time].recovered  # Order is altered for visualisation
            tseries[itr_time, 2] = num_timeseries[itr_time].infected
            tseries[itr_time, 3] = num_timeseries[itr_time].not_infected
        end

        p = groupedbar(
            tseries,
            bar_position=:stack,
            bar_width=0.7,
            color = [:gold :orangered :deepskyblue],
            label = ["Recovered" "Infected" "Not_infected"],
            xaxis = ("Time step"),
            yaxis = ("Number of particles"),
            size=(960, 640),
            )
        savefig(p, filename)
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
stdout_condition,
plot_particles,
make_gif,
plot_num_timeseries

# ----------------------------------------
## Set parameters & variables
# ----------------------------------------

### Declare parameters
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

param = ParamVar.Parameters(
    num_particles,max_iteration,
    x_range,y_range,
    vel_mean,vel_σ,vel_flc,
    ratio_infection_init,recovery_time,infection_chance,
    radius_infection)

### Declare flags
flag_multiple_infection = true
flag_infected_isolation = false

flag = ParamVar.Flags(
    flag_multiple_infection,
    flag_infected_isolation
)

### Stdout simulation condition & define output directory name
out_dir = stdout_condition(param, flag)

### Define array of particle properties
particles = Array{ParamVar.Particle}(undef, param.num_particles)
for itr_ptcl = 1:param.num_particles
    particles[itr_ptcl] = ParamVar.Particle()
end
# particles .= ParamVar.Particle()  # ERROR: LoadError: MethodError: no method matching length(::Main.ParamVar.Particle)

### Define array of number of not_infected/infected/recovered particle in a snapshot
not_infected, infected, recovered = 0, 0, 0

num_snapshot = ParamVar.NumSnapshot(
    not_infected, infected, recovered)

### Define array of number of not_infected/infected/recovered particle in a timeseries
num_timeseries = Array{ParamVar.NumSnapshot}(undef, 0)  # push! in temporal iteration

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
    update_particles(param, flag, particles)

    # Count not-infected, infected & recovered number of particles
    count_status(param, particles, num_snapshot, num_timeseries)

    # Plot particles for gif video
    plot_particles(itr_time, param, num_snapshot, out_dir, particles)

    # tmp_string = @sprintf "itr_time = %i x[1] = %6.3f y[1] = %6.3f" itr_time particles[1].pos_x particles[1].pos_y
    # println(tmp_string)
    # println("itr_time = ", itr_time, " not infected = ", num_snapshot.not_infected, " infected = ", num_snapshot.infected, " recovered = ", num_snapshot.recovered)
    next!(progress)

    # Finish if there are no infected particles any more
    if num_snapshot.infected == 0
        println("\n No patients at itr_time = ", itr_time, " :Exit")
        break
    end
end


# ----------------------------------------
## Make gif video of simulation
# ----------------------------------------
make_gif(param, anim, out_dir)

# ----------------------------------------
## Make figure of timeseries
# ----------------------------------------
plot_num_timeseries(param, num_timeseries, out_dir)
