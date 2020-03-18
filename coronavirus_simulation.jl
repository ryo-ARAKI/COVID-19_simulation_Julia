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
        recovery_chance::Float64  # Chance of recovery per unit time
        infection_chance::Float64  # Chance of infection per unit time
    end

    mutable struct Particle
        pos_x::Float64
        pos_y::Float64
        state::Char  # 'g':never been infected, 'r':infected, 'o':had been infected
        t_inf::Int64
        flag_inf::Bool
        # Constructor
        Particle() = new()
    end
end


# ========================================
# Main function
# ========================================

## Declare modules
using Distributions
using Plots
using .ParamVar

font = Plots.font("Times New Roman", 20)

# ----------------------------------------
## Set parameters & variables
# ----------------------------------------
num_particles = 10  # 100
max_iteration = 100

x_range = 10.0
y_range = 10.0

recovery_chance = 0.30
infection_chance = 0.03

### Declare parameters
param = ParamVar.Parameters(
    num_particles,max_iteration,
    x_range,y_range,
    recovery_chance,infection_chance)

### Define array of particle properties
particles = Array{ParamVar.Particle}(undef, param.num_particles)
for itr_par = 1:param.num_particles
    particles[itr_par] = ParamVar.Particle()
end
# particles .= ParamVar.Particle()  # ERROR: LoadError: MethodError: no method matching length(::Main.ParamVar.Particle)

# ----------------------------------------
## Set initial condition of particles
# ----------------------------------------
for itr_par = 1:param.num_particles
    particles[itr_par].pos_x = rand(Uniform(-0.5*param.x_range, 0.5*param.x_range))
    particles[itr_par].pos_y = rand(Uniform(-0.5*param.y_range, 0.5*param.y_range))
    particles[itr_par].state = 'g'  # Initially not infected
    if itr_par == 1  # One particle is initially infected
        particles[itr_par].state = 'r'
    end
    particles[itr_par].t_inf = 0
    particles[itr_par].flag_inf = false
end
# particles.pos_x .= rand(Uniform(-0.5*param.x_range, 0.5*param.x_range))  # ERROR: LoadError: type Array has no field pos_x

#=
println("x = ", getfield.(particles, :pos_x))
println("y = ", getfield.(particles, :pos_y))
println("state = ", getfield.(particles, :state))
println("t_inf = ", getfield.(particles, :t_inf))
println("flag_inf = ", getfield.(particles, :flag_inf))
=#
