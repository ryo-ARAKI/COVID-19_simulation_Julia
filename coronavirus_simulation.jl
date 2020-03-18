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
        x_range::Float64  # Size of computational domain [-x_range/2, x_range/2]
        y_range::Float64
        #
        recovery_time::Int64  # Recovery time
        infection_chance::Float64  # Chance of infection per unit time
    end

    mutable struct Particle
        pos_x::Float64
        pos_y::Float64
        state::Char  # 'g':never been infected, 'r':infected, 'o':had been infected
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
    """
    Update particle properties
    """
    function update_particles(param, ptcl)
        x_ifcn = Array{Float64}(undef, 0)
        y_ifcn = Array{Float64}(undef, 0)

        # Update properties of infected particles
        for itr_ptcl = 1:param.num_particles
            # Recovery from infection
            s = ptcl[itr_ptcl].state
            t = ptcl[itr_ptcl].t_ifcn
            if s == 'r' && t >= param.recovery_time
                ptcl[itr_ptcl].state = 'o'  # Hold infection history
                ptcl[itr_ptcl].flag_ifcn = true  # Hold infection history
            end
            # Store position of infected particles
            if s == 'r'
                append!(x_ifcn, ptcl[itr_ptcl].pos_x)
                append!(y_ifcn, ptcl[itr_ptcl].pos_y)
            end
        end

        println("x_ifcn = ", x_ifcn, " y_ifcn = ", y_ifcn)
    end
end


# ========================================
# Main function
# ========================================

## Declare modules
using Distributions
using Plots
using .ParamVar
using .TimeMarch:
update_particles

font = Plots.font("Times New Roman", 20)

# ----------------------------------------
## Set parameters & variables
# ----------------------------------------
num_particles = 10  # 100
max_iteration = 10  # 100

x_range = 10.0
y_range = 10.0

recovery_time = 3  # 30
infection_chance = 0.03

### Declare parameters
param = ParamVar.Parameters(
    num_particles,max_iteration,
    x_range,y_range,
    recovery_time,infection_chance)

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
    particles[itr_ptcl].pos_x = rand(Uniform(-0.5*param.x_range, 0.5*param.x_range))
    particles[itr_ptcl].pos_y = rand(Uniform(-0.5*param.y_range, 0.5*param.y_range))
    particles[itr_ptcl].state = 'g'  # Initially not infected
    if itr_ptcl == 1  # One particle is initially infected
        particles[itr_ptcl].state = 'r'
    end
    particles[itr_ptcl].t_ifcn = 0
    particles[itr_ptcl].flag_ifcn = false
end
# particles.pos_x .= rand(Uniform(-0.5*param.x_range, 0.5*param.x_range))  # ERROR: LoadError: type Array has no field pos_x

#=
println("x = ", getfield.(particles, :pos_x))
println("y = ", getfield.(particles, :pos_y))
println("state = ", getfield.(particles, :state))
println("t_ifcn = ", getfield.(particles, :t_ifcn))
println("flag_ifcn = ", getfield.(particles, :flag_ifcn))
=#

# ----------------------------------------
## Time iteration of infection simulation
# ----------------------------------------
println("x = ", particles[1].pos_x)
for itr_time = 1:param.max_iteration
    println("itr_time = ", itr_time)
    println("x = ", particles[1].pos_x)
    update_particles(param, particles)
end
