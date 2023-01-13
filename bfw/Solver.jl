
T = 0.001
Num = 9
Beads = 10
include("Physics.jl")
include("Fermions.jl")

function V_r(P) #势能 只有一个正电中心
    N,dim,B = size(P)
    V = 0.0
    for b in 1:B
        for i in 1:N    #+
            V += -(Num)/sqrt(sum(abs2,P[i,:,b]))
        end
        for i in 2:N    #-
            for j in 1:i-1
                V += 1/sqrt(sum(abs2,P[i,:,b]-P[j,:,b]))
            end
        end
    end
    for b=1:B
        for j=1:N
            for i=1:N
                V += sum(abs2,P[i,:,b].-P[j,:,o(b,B)]).*(k*B)
            end
        end
    end
    return V/B
end

function Ep(p)  #总的势能
    N,dim,B = size(p)
    nc = size(p,2)
    Ep = zeros(nc)
    for i in 1:nc
        P = reshape(p[:,i],N,3,B)
        W_f =W_bf(P)   #交换虚拟势能
        V = V_r(P)/B
        Ep[i] = -1 *((W_f <0 ? -log(-W_f) : -log(W_f)) + V)
    end
    return Ep
end

function Ep_debug(p)
    nc = size(p,2)
    lp = zeros(nc,2)
    for i in 1:nc
        P = reshape(p[:,i],N,3,B)
        W_f = W_bf(P)
        V = V_r(P)/B
        lp[i,1] = -1 *((W_f <0 ? -10e2*log(-W_f) : -log(W_f)))
        println(W_f," ",lp[i,1]," ",V)
    end
end

using AdvancedHMC
using ForwardDiff

function Enging_HMC(P,E_P)
    initial_θ = P[:]
    metric = DiagEuclideanMetric(length(initial_θ))
    hamiltonian = Hamiltonian(metric, E_P, logDensity_deriv)
    # Define a leapfrog solver, with initial step size chosen heuristically
    initial_ϵ = find_good_stepsize(hamiltonian, initial_θ)
    integrator = Leapfrog(initial_ϵ)
    proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
    adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))
    samples, stats = sample(hamiltonian, proposal, initial_θ, n_samples, adaptor; progress=true)
end

function ln_E(P)
    Pos = reshape(P,Num,3,Beads)
    return(-β*V_r(Pos)
            #+log(abs(W_bf(Pos[1:fld(Num,2),:,:])))
            #+log(abs(W_bf(Pos[fld(Num,2)+1:Num,:,:])))
            )
end


# ForwardDiff.gradient(pl -> W_bf(reshape(pl,N,3,B)), pl)

# samples, stats = Enging_HMC(rand(N,3,B), x -> lnE(x))

using LinearAlgebra
using AdvancedMH
using Distributions

function AMC(P,steps::Real)

    density_model = AdvancedMH.DensityModel(x -> ln_E(x))
    # proposal = Array{Normal{Float64}}(undef, N,dim,beads)

    # for i in 1:N
    #     proposal[i,:,:] .= Normal(0, 1)
    # end
    proposal = RandomWalkProposal(MvNormal(zeros(Num,3,Beads)[:], I))
    sampler = AdvancedMH.MetropolisHastings(proposal)

    chain = AdvancedMH.sample(density_model, sampler, convert(Int, steps);
                            init_params=P[:])

    return chain
end

sign = 0
chain = AMC(rand(Num,3,Beads),5000)
configs = [reshape(config.params,Num,3,Beads) for config in chain]
E = [config.lp for config in chain]./(-β)

using Plots
plot(E)
