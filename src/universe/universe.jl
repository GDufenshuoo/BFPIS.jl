




function setting_universe(simulator::Simulator, method::Method)

    if simulator.simulator_method == "PI"
        PI_universe!(simulator, method, universe)
    end

    if simulator.simulator_tool == "MC"
         MC_universe!(simulator, method, universe)
    end

    return universe
end

function PI_universe!(simulator::Simulator, method::Method, universe::Universe)
    
end

function MC_universe!(simulator::Simulator, method::Method, universe::Universe)
    
end
