# COVID-19_simulation_Julia

COVID-19 spreading simulation, Julia translation of https://github.com/MuAuan/collective_particles with several changes and new functionalities.

<img src="https://raw.githubusercontent.com/ryo-ARAKI/COVID-19_simulation_Julia/demo/gif/particles.gif" width="400"> <img src="https://raw.githubusercontent.com/ryo-ARAKI/COVID-19_simulation_Julia/demo/gif/timeseries.png" width="450">

## Model description

### Parameter

- `num_particles` particles are simulated over `max_iteration` temporal iteration.
- Computational domain is `[0, x_range]` × `[0, y_range]`.
- Initial particle velocity `(vel_x, vel_y)` obeys Gaussian distribution, with average `vel_mean` and standard deviation `vel_σ`.
- `ratio_infection_init` percentage of particles are initially infected by the virus.
- Particles in the nearby (below relative distance `radius_infection`) of the infected particle can be infected by the chance rate `infection_chance` per time step.
- Infected particles (deterministically) recover after `recovery_time` steps and never infected again.
- During temporal development, particle velocity alters by up to `vel_flc`.

### Flag

- `flag_multiple_infection` allows recovered particles to re-infect the virus.
  - In the case of multiple infection, `infection_chance` is modified as follows.

```math
  \textrm{multiple infection chance} = \textrm{infection chance} \times 0.5^\textrm{number of past infection}
```

- `flag_infected_isolation` forces infected particles to stop their movement during infection.

## Execution

- The code is tested on Julia version 1.3.1.
- `./fig/` subdirectory is needed for gif video.
- Module dependencies are [ProgressMeter](https://github.com/timholy/ProgressMeter.jl), [Distributions](https://github.com/JuliaStats/Distributions.jl), [Printf](https://github.com/JuliaLang/julia/blob/master/stdlib/Printf/src/Printf.jl) and [Plots](http://docs.juliaplots.org/latest/) with [GR](https://gr-framework.org/julia.html) backend.
