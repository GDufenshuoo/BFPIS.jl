# a new program to simulate femi and bose
# started at 2022 2.11 by zlwang/clockwang
# Todo / finished
#
# quantum mente carlo   |- base programming
#                       |- path
#                       |- worm
#
# molecular dynamics    |- base
#                       |- path
#                       |- bose
#
# Other things 
# DFT support
#
# This part offers the final API to simulate

"""
    simulate!(system, simulator, method(MC))        --- programming
    simulate!(system, simulator, method(MD))        --- TODO

Run simulations using simulator with different method.
- system contains all the physical information
- simulator contains MD/MC DFT/PI 
- method contains n_step neighbors_settings and some things verious from methods to methods 
which is Artificially and deliberately added

"""
function simulate!(system::System, simulator::Simulator, method::Method)

    universe = setting_universe(simulator, method)
    run!(system, universe)
    observe(system, method.obs)

end

