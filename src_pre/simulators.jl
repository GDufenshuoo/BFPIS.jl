# a new program to simulate femi and bose
# started at 2022 2.11 by zlwang/clockwang
# Todo / finished
#
#

"""
    System(; <keyword arguments>)

A physical system to be simulated.
Properties unused in the simulation or in analysis can be left with their
default values.
`atoms`, `atoms_data`, `coords` and `velocities` should have the same length.
This is a sub-type of `AbstractSystem` from AtomsBase.jl and implements the
interface described there.

# Arguments
- `atoms::A`: the atoms, or atom equivalents, in the system. Can be
    of any type but should be a bits type if the GPU is used.
- `atoms_data::AD`: other data associated with the atoms, allowing the atoms to
    be bits types and hence work on the GPU.
- `pairwise_inters::PI=()`: the pairwise interactions in the system, i.e.
    interactions between all or most atom pairs such as electrostatics.
    Typically a `Tuple`.
- `specific_inter_lists::SI=()`: the specific interactions in the system,
    i.e. interactions between specific atoms such as bonds or angles. Typically
    a `Tuple`.
- `general_inters::GI=()`: the general interactions in the system,
    i.e. interactions involving all atoms such as implicit solvent. Typically
    a `Tuple`.
- `coords::C`: the coordinates of the atoms in the system. Typically a
    vector of `SVector`s of 2 or 3 dimensions.
- `velocities::V=zero(coords)`: the velocities of the atoms in the system.
- `box_size::B`: the size of the box in which the simulation takes place.
    Typically a `SVector` of 2 or 3 dimensions.
- `neighbor_finder::NF=NoNeighborFinder()`: the neighbor finder used to find
    close atoms and save on computation.
- `loggers::L=Dict()`: the loggers that record properties of interest during a
    simulation.
- `force_units::F=u"kJ * mol^-1 * nm^-1"`: the units of force of the system.
- `energy_units::E=u"kJ * mol^-1"`: the units of energy of the system.
- `gpu_diff_safe::Bool`: whether to use the code path suitable for the
    GPU and taking gradients. Defaults to `isa(coords, CuArray)`.
"""
mutable struct Simulator{SM, ST, PD, PC, V, B, NF, L, F, E}
    simulator_method::SM
    simulator_setting::ST
    particle_data::PD
    path_coords::PC
    velocities::V
    box_size::B
    neighbor_finder::NF
    force_units::F
    energy_units::E
    loggers::L
end

function MD!(system::System, universe::Universe)
    0
end



"""
# basic MC MC_move

 `MC_move!(universe.particle_list, distribution, E_0, MConf_setting)`
   | - `MC_distribution(E_conf, E)`
   | - `MC_generate_conf(MC_distribution, MC_movelength)`
 `density!(system.loggers.density, universe.conf)`

## change the  
- universe.conf
- universe.E_0

## give out `pass` and `cont`

"""
function MC!(system::System, universe::Universe)

    pass = 0.0
    cont = 0.0
    
    universe.conf = system.particle_position

    while pass < universe.block_step
        cont += MC_move!(universe.particle_list, universe.conf, universe.E_0, universe.MConf_setting)
        pass += 1
        E = MC_step[2]
        density!(system.loggers.density, universe.conf)
    end

    return pass, cunt
end


"""
# basic MC MC_move
-  `MC_move!(universe.particle_list, distribution, E_0, MConf_setting)`
-   | - `MC_distribution(E_conf, E)`
-   | - `MC_generate_conf(MC_distribution, MC_movelength)`

## change the  
- distribution 
- E_0

## give out `cont`

"""
function MC_move!(universe.particle_list, distribution, E_0, MConf_setting)
    i = 0
    cont = 0
    MC_step = false, E_0
    while i < length(distribution)
        MC_step[1] = false
        while MC_step[1] = false
            universe.fur_conf[i] = MC_generate_conf(distribution[i], MConf_setting)
            MC_step = MC_distribution(MC_energy(universe.fur_conf, i, universe.particle_list), E_0[i])
            cont += 1
        end
        E_0[i] = MC_step[2] 
        i += 1
    end
    distribution = universe.fur_conf
    return cont
end

"""
the part which decide which to be accpeted
according to Blotzmann distribution
"""
function MC_distribution(E_conf, E)
    if E_conf <= E || exp(E_conf - E) <= rand()
        return true, E_conf
    else
        return false, E
    end
end

"""
MC generate configuration 
"""
function MC_generate_conf(MC_distribution, MC_movelength)
    i = 0
    while i < length(distribution)
        MC_distribution[i] += (rand()-0.5)*MC_movelength
        i += 1
    end
    return MC_distribution
end




