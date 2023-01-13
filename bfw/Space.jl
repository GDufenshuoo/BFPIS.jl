using Hungarian
include("Physics.jl")

"""
Create a circle of 1-N and move backrward
"""
function o(i, N)
	if i == 1
		return N
	else
		return i-1
	end
end

"""
Create a circle of 1-N and not move
"""
function ro(i, N)
	if i > N
		return i%N
    elseif i < 1
		return i%N+N
    else
        return i
	end
end

"""
Create a circle of 1-N and move backrward
"""
function fo(i, N)
	if i == N
		return 1
	else
		return i+1
	end
end

"""

"""
function bosons()
end

"""
Forward fermions sampling
"""
function fermions(P)
    d2,d2_s = d2(P)
    db2,db2_s = db2(P)
end

"""

Give in Position N×3×B and 
out the d² between N particles N×N×B:
"""
function d2(P)
    N = size(P)[1]
    B = size(P)[3]

    d2 = zeros(N,N,B)
    for i=1:N
        d2[i,i,:].=1e33
    end
    for b=1:B
        for i=2:N
            for j=1:i-1
                d = sum(abs2,P[i,:,b].-P[j,:,b])
                d2[i,j,b] = d
                d2[j,i,b] = d
            end
        end
    end
    return d2
end

"""
Reshape and calculate the Path
"""
function E_N(P)
    N = size(P)[1]
    B = size(P)[3]

    E_d2 = zeros(N,N+1)
    db² = db2(P)

    permute = zeros(Int,N,B)
    for b in 1:B
        permute[:,b] .= hungarian(db²[:,:,b])[1]
    end

    # 1 <- B
    lp = zeros(Int,N,B)
    for i in 1:N
        n = i
        for j in 1:B
            n = permute[n,j]
            lp[i,j] = n
        end
    end

    # limitition: B > 1
    for b in 2:B
        for i in 1:N
            E_d2[i,1] += db²[lp[i,o(b,B)],lp[i,b],b]
            # println("$(db²[lp[i,o(b,B)],lp[i,b],b]) $i $(lp[i,o(b,B)]) $(lp[i,b]) $b")
        end
    end

    # <- 1 ~ B <-
    for i in 1:N, j in 1:N
        E_d2[N-i+1,1+j] += db²[lp[ro(j+i-1,N),B],lp[j,1],1]
        # println("N-i+1 = $(N-i+1) j+i-1 = $(j+i-1) j:$(j) i:$(i)    ",
        # "db²[lp[$(ro(j+i-1,N))]=$(lp[ro(j+i-1,N),B]),$(lp[j,1])]")
    end

    return E_d2
end

sort

"""
Reshape and calculate the exchange Path

`i,j -> -j-i- (B-1 -> B)`'
"""
function db2(P)
    N = size(P)[1]
    B = size(P)[3]
    db2 = zeros(N,N,B)

    for b=1:B
        for j=1:N
            for i=1:N
                d = sum(abs2,P[i,:,b].-P[j,:,o(b,B)])
                db2[i,j,b] = d
            end
        end
    end
    return db2
end
