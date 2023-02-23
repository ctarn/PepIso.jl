module PepIso

import GLPK
import JuMP
import MesCore

include("IPV.jl")
using .IPV: build_ipv, calc_ipv, ipv_m, ipv_w, ipv_mz

prefilter(ion, spec, ε, V, max_mode=false) = begin
    f = x -> (x > 0) && !isempty(MesCore.query_ε(spec, IPV.ipv_mz(ion, x, V), ε))
    if max_mode
        i = argmax(IPV.ipv_w(ion, V))
        return f(i) && (f(i + 1) || f(i - 1))
    else
        return f(1) && f(2)
    end
end

exclude(ions, xs, ms, τ, ε, V) = begin
    map(zip(ions, xs, ms)) do (i, x, m)
        (x <= 0.0 || m <= 0.0) && return false
        for (i_, x_) in zip(ions, xs)
            if i_.z % i.z == 0 && x < x_ * τ
                MesCore.in_moe(i.mz, i_.mz, ε) && i_ !== i && return false
                MesCore.in_moe(i.mz, ipv_mz(i_, 2, V), ε) && return false
                MesCore.in_moe(i.mz, ipv_mz(i_, 3, V), ε) && return false
            end
        end
        return true
    end
end

build_constraints_lp(ions, spec, ε, V) = begin
    items = [(; i, j, mz=ipv_mz(ion, j, V)) for (i, ion) in enumerate(ions) for j in eachindex(ipv_w(ion, V))]
    cs = map(p -> (; p.mz, p.inten, slots=empty(items)), spec)
    for i in items
        l, r = searchsortedfirst(spec, (1 - ε) * i.mz), searchsortedlast(spec, (1 + ε) * i.mz)
        εs = map(p -> abs(i.mz - p.mz), spec[l:r])
        l <= r && push!(cs[argmin(εs)+l-1].slots, i)
        l > r && push!(cs, (; i.mz, inten=0.0, slots=[i]))
    end
    return filter(c -> !isempty(c.slots), cs)
end

solve_lp(ions, cs, V) = begin
    model = JuMP.Model(GLPK.Optimizer)
    JuMP.set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF)
    JuMP.@variable(model, x[eachindex(ions)] >= 0)
    JuMP.@variable(model, δ[eachindex(cs)] >= 0) # the constraint is redundant but faster
    JuMP.@objective(model, Min, sum(δ))
    s = map(c -> sum(s -> x[s.i] * ipv_w(ions[s.i], V)[s.j], c.slots), cs)
    i = map(c -> c.inten, cs)
    JuMP.@constraint(model, i .- s .<= δ)
    JuMP.@constraint(model, s .- i .<= δ)
    JuMP.optimize!(model)
    return (xs=JuMP.value.(x), δs=i .- JuMP.value.(s))
end

calc_match_lp(xs, cs, δs) = begin
    map(eachindex(xs)) do i
        s = map(c -> any(s -> s.i == i, c.slots), cs)
        return 1 - sum(abs, δs[s], init=0.0) / sum(c -> c.inten, cs[s], init=1.0e-16)
    end
end

evaluate_lp(ions, spec, ε, V) = begin
    cs = build_constraints_lp(ions, spec, ε, V)
    xs, δs = solve_lp(ions, cs, V)
    ms = calc_match_lp(xs, cs, δs)
    return xs, ms
end

build_constraints_rlp(ions, spec, ε, V) = begin
    items = [(; i, j, mz=ipv_mz(ion, j, V)) for (i, ion) in enumerate(ions) for j in eachindex(ipv_w(ion, V))]
    cs = map(p -> (; p.mz, p.inten, slots=[]), spec)
    tab = [[] for _ in ions]
    for item in items
        l, r = searchsortedfirst(spec, (1 - ε) * item.mz), searchsortedlast(spec, (1 + ε) * item.mz)
        for (k, c) in enumerate(cs[l:r])
            push!(c.slots, (; item..., k))
        end
        l > r && push!(cs, (; item.mz, inten=0.0, slots=[(; item..., k=1),]))
        push!(tab[item.i], l <= r ? Vector(l:r) : [length(cs)])
    end
    return cs, tab
end

solve_rlp(ions, cs, tab, V) = begin
    model = JuMP.Model(GLPK.Optimizer)
    JuMP.@variable(model, x[eachindex(ions)] >= 0)
    JuMP.@variable(model, u[i=eachindex(tab), j=eachindex(tab[i]), k=eachindex(tab[i][j])] >= 0)
    JuMP.@variable(model, δ[eachindex(cs)] >= 0) # the constraint is redundant but faster
    JuMP.@objective(model, Min, sum(δ))
    s = map(c -> sum(s -> u[s.i, s.j, s.k] * ipv_w(ions[s.i], V)[s.j], c.slots, init=0.0), cs)
    i = map(c -> c.inten, cs)
    JuMP.@constraint(model, i .- s .<= δ)
    JuMP.@constraint(model, s .- i .<= δ)
    for i in eachindex(tab)
        for j in eachindex(tab[i])
            JuMP.@constraint(model, sum(u[i, j, k] for k in eachindex(tab[i][j])) == x[i])
        end
    end
    JuMP.optimize!(model)
    return (xs=JuMP.value.(x), δs=i .- JuMP.value.(s), us=JuMP.value.(u))
end

calc_match_rlp(xs, cs, δs, us, tab) = begin
    map(enumerate(xs)) do (i, x)
        inten = 0.0
        delta = 0.0
        for (j, ps) in enumerate(tab[i])
            ws = [us[i, j, k] / max(x, 1.0e-16)  for k in eachindex(ps)]
            inten += sum(map(p -> cs[p].inten, ps) .* ws)
            delta += sum(map(p -> abs(δs[p]), ps) .* ws)
        end
        return 1 - delta / max(inten, 1.0e-16)
    end
end

evaluate_rlp(ions, spec, ε, V) = begin
    cs, tab = build_constraints_rlp(ions, spec, ε, V)
    xs, δs, us = solve_rlp(ions, cs, tab, V)
    ms = calc_match_rlp(xs, cs, δs, us, tab)
    return xs, ms
end

evaluators = Dict(
    :LP => evaluate_lp,
    :RLP => evaluate_rlp,
)

deisotope(ions, spec, τ_max, ε, V, evaluator=:LP) = begin
    evaluator = evaluators[evaluator]
    xs = ms = nothing
    τs = [τ_max / 4, τ_max / 2, τ_max]
    τ = popfirst!(τs)
    while true
        xs, ms = evaluator(ions, spec, ε, V)
        s = exclude(ions, xs, ms, τ, ε, V)
        while all(s)
            isempty(τs) && @goto done
            τ = popfirst!(τs)
            s = exclude(ions, xs, ms, τ, ε, V)
        end
        ions = ions[s]
    end
    @label done
    return [(; i.mz, i.z, x, m) for (i, x, m) in zip(ions, xs, ms)]
end

end
