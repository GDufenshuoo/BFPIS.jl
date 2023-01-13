# a new program to simulate femi and bose
# started at 2022 2.11 by zlwang/clockwang
# Todo / finished
#
# |-
#
# This part offers the System structure 


"""
    System(; <keyword arguments>)

A physical system.

`atoms`, `atoms_data`, `coords` and `velocities` should have the same length.
This is a sub-type of `AbstractSystem` from AtomsBase.jl and implements the
interface described there.

# Arguments
- `particle_position`
- `particle_type`
- `particle_name`
- `particle_state`
- `path`
- `velocities`
- `box`
- `loggers`: record properties 
"""
mutable struct System{Pp, PT, Pn, Ps, P, V, B, L} <: AbstractSystem
    particle_position::Pp
    particle_type::PT
    particle_name::Pn
    particle_state::Ps
    path::P
    velocities::V
    box::B
    loggers::L
end


