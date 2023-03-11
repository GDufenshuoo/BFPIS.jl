# using("Physics.jl")
include("Space.jl")


"""
According to 
    “Path Integral Molecular Dynamics for Fermions...”

    -It's not necessary to consider all permutation.
"""
function W_bf(P)
    N = size(P)[1]
    B = size(P)[3]
    W_bf=zeros(BigFloat,N+1,N)
    W_bf[1,:].=1.0
    En = E_N(P)
    Rp = sortperm(En[:,1];rev = false)
    Eβ = big(-k*β).*En
    for ln in 1:N
        for n=1:N
            for k=1:n
                W_bf[n+1,ln] += (-1)^(k-1)*(
                    exp(Eβ[Rp[ro(N-k+ln,N)],1]+Eβ[Rp[N-k+1],Rp[ln]+1]) 
                    * W_bf[n+1-k,ln]
                    )/n
                # println("$n $ln $k exp($(Float32(Eβ[Rp[ro(N-k+ln,N)],1])) + $(Float32(Eβ[Rp[N-k+1],Rp[ln]+1]))) ")
                # println("ip:$(ip) ln:$(ln) n:$n k:$k")
            #debug 
            # println(" $(-Eβ[Rp[ro(N-k+ln,N)],1])+$(Eβ[Rp[N-k+1],Rp[ln]+1]) ")
            end
        end
    end
    println("\n \n sum(W_bf[N+1,:]) $(sum(W_bf[N+1,:]))")
    return W_bf#sum(W_bf[N+1,:]/N)#,Rp,Eβ
end


function D_W_bf(P)
    N = size(P)[1]
    B = size(P)[3]
    W_bf=zeros(BigFloat,N+1,N)
    D_W_bf=zeros(BigFloat,N+1,N)
    W_bf[1,:].=1.0
    Eβ = k*β.*E_N(P)
    for ln in 1:N
        for n=1:N
            for k=1:n
                W_bf[n+1,ln] += (-1)^(k-1)*(exp(-Eβ[ro(N-k+ln,N),1]+Eβ[N-k+1,ln+1])* W_bf[n+1-k,ln])/n
            #debug 
            #println(W_bf[n+1]," ",-exp(-big(β)*_E_N[k])," ",W_bf[n+1-k],"\n")
            end
        end
    end
    Eβ *= big(1/β*(β+1e-12β))
    for ln in 1:N
        for n=1:N
            for k=1:n
                D_W_bf[n+1,ln] += (-1)^(k-1)*(exp(-Eβ[ro(N-k+ln,N),1]+Eβ[N-k+1,ln+1])* W_bf[n+1-k,ln])/n
            #debug 
            #println(W_bf[n+1]," ",-exp(-big(β)*_E_N[k])," ",W_bf[n+1-k],"\n")
            end
        end
    end
    return β.*(-β⁻¹*log.(D_W_bf).+β⁻¹log.(W_bf))./(1e-12β) +(-β⁻¹log.(W_bf))
end
