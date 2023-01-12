

# oS = 0.0
# oSn = 0.0
# l = 0
# for i in 1:1000
#     a = W_bf(rand(N,3,B))
#     oS += a/1000
#     for k in a
#         if k < 0
#             println(k," ",i)
#             l += 1
#         end
#     end
# end

# show("oS:$oS l:$l")


T = 0.000000000000000000001
N = 9
B = 10
include("Physics.jl")
include("Fermions.jl")

function V_r(P) #势能 只有一个正电中心
    V = 0.0
    for b in 1:B
        for i in 1:N    #+
            V += -1/sqrt(sum(abs2,P[i,:,b])) *1
        end
        # for i in 2:N    #-
        #     for j in 1:i-1
        #         V += 1/sqrt(sum(abs2,P[i,:,b]-P[j,:,b])) *E_epr
        #     end
        # end
    end
    V
end

function Ep(p)  #总的势能
    nc = size(p,2)
    Ep = zeros(nc)
    for i in 1:nc
        P = reshape(p[:,i],N,3,B)
        W_f =W_bf(P)   #交换虚拟势能
        V = V_r(P)/B
        Ep[i] = -1 *((W_f <0 ? -β*log(-W_f) : -log(W_f)) + V)
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



import AbstractMCMC as AMC
using AdvancedPS



# Obtain the initial sample and state.
sample, state = AMC.step(rng, model, sampler)

# Save the sample.
samples = AMC.samples(sample, model, sampler, N)
samples = AMC.save!!(samples, sample, 1, model, sampler, N)

# Step through the sampler.
for i in 2:Sn
    # Obtain the next sample and state.
    sample, state = AMC.step(rng, model, sampler, state)

    # Save the sample.
    samples = AMC.save!!(samples, sample, i, model, sampler, N)
end

AMC.bundle_samples(samples, model, sampler, state, chain_type)

Wbf = AMC.LogDensityModel(P -> W_bf(P))
AMC.Sample(Wbf, AdvancedPS.PG(N))
