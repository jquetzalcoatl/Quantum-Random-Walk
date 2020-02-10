using Plots, LinearAlgebra, DelimitedFiles, Arpack
plotlyjs()
# using Plotly
mutable struct Site
    pos
    occ
    prob
    spin
end
Ldim=20
sites=Array{Site}(undef,Ldim)

function initsys()
    for i=1:Ldim
        if Ldim/2 == i  # || Ldim/2 -2 == i
            # global sites[i] = Site(-Ldim/2+i, [1], [1], [1/sqrt(2),1/sqrt(2)])
            global sites[i] = Site(-Ldim/2+i, [1], [1.0], [1.,0.])
        else
            global sites[i] = Site(-Ldim/2+i, [0], [0.0], [0.,0.])
        end
    end
    println("init sys...")
end

function findnonzero(m)
    aRight = similar(sites, Int)
    aLeft = similar(sites, Int)
    countR, countL = 1, 1
    @inbounds for i in eachindex(sites)
        aRight[countR] = i
        aLeft[countL] = i
        countR += (sites[i].occ[end] != 0. && sum(m[1] * sites[i].spin) != 0)
        countL += (sites[i].occ[end] != 0. && sum(m[2] * sites[i].spin) != 0)
    end
    return resize!(aRight, countR-1), resize!(aLeft, countL-1)
end

function translateFunc(m)
    vR, vL = findnonzero(m)
    println(vR, vL)
    right = vR .+ 1
    left = vL .- 1
    maketrans(vR, vL, right,left)
    println("translation...")
end

function maketrans(vR, vL, right,left)
    @inbounds for i in eachindex(right)
        append!(sites[right[i]].occ, 1)
        # unNorm = m[1] * sites[vR[i]].spin
        # sites[right[i]].spin = unNorm/sqrt(unNorm' * unNorm)
        sites[right[i]].spin = m[1] * sites[vR[i]].spin
        append!(sites[right[i]].prob, (sites[right[i]].spin' * sites[right[i]].spin) )
    end
    @inbounds for i in eachindex(left)
        if sum(left[i] .== right) == 0
            append!(sites[left[i]].occ, 1)
            # unNorm = m[2] * sites[vL[i]].spin
            # sites[left[i]].spin = unNorm/sqrt(unNorm' * unNorm)
            sites[left[i]].spin = m[2] * sites[vL[i]].spin
            append!(sites[left[i]].prob, (sites[left[i]].spin' * sites[left[i]].spin) )
        else
            # append!(sites[left[i]].occ, 1)
            # unNorm = sites[left[i]].spin + m[2] * sites[vL[i]].spin
            # sites[left[i]].spin = unNorm/sqrt(unNorm' * unNorm)
            sites[left[i]].spin += m[2] * sites[vL[i]].spin
            sites[left[i]].prob[end] = (sites[left[i]].spin' * sites[left[i]].spin)
        end
    end
    s = setdiff(collect(1:Ldim),right)
    s = setdiff(s,left)
    @inbounds for i in s            # THIS MAKES SENSE WHEN THE PROB TO STAY IN STATE i IS ZERO!!!
        append!(sites[i].occ, 0)
        append!(sites[i].prob, 0 )
        sites[i].spin = [0.0, 0.0]
    end
end

function cointoss(h)
    @inbounds for i in eachindex(sites)
        if sites[i].occ[end] > 0
            # unNorm = h * sites[i].spin
            # sites[i].spin = unNorm/sqrt(unNorm' * unNorm)
            sites[i].spin = h * sites[i].spin
        end
    end
    println("coin toss...")
end

m=[[1 0;0 0],[0 0;0 1]]
h=1/sqrt(2) * [1 1; 1 -1]
initsys()
println([sites[i].spin for i=8:12])
println([sites[i].prob for i=8:12])
cointoss(h)
println([sites[i].spin for i=8:12])
println([sites[i].prob for i=8:12])
translateFunc(m)
cointoss(h)
