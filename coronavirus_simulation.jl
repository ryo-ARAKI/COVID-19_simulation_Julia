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
end


# ========================================
# Main function
# ========================================

## Declare modules
using .ParamVar
using Plots

font = Plots.font("Times New Roman", 20)

# ----------------------------------------
## Set parameters & variables
# ----------------------------------------
num_particles = 100
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
