using OrdinaryDiffEq
using Plots
plotly()

function getReplicatorTrajectory()
    ntypes = 5
    A = rand(ntypes, ntypes)
    x = rand(ntypes)
    x ./= sum(x)
    g(x, A) = x .* (A * x .- x' * A * x)
    Δt = 0.1
    T = 100.00
    trange = 0.0:Δt:T
    nt = length(trange)
    Xt = zeros(nt, ntypes)
    Xt[1, :] = x
    for i in 2:nt
        x .+= g(x, A) .* Δt
        Xt[i, :] = x
    end
    return trange, Xt
end

function getReplicatorTrajectoryWithSumxDots()
    ntypes = 5
    A = rand(ntypes, ntypes)
    x = rand(ntypes)
    x ./= sum(x)
    g(x, A) = x .* (A * x .- x' * A * x)
    Δt = 0.1
    T = 100.00
    trange = 0.0:Δt:T
    nt = length(trange)
    Xt = zeros(nt, ntypes)
    sumxdots = zeros(nt)
    Xt[1, :] = x
    for i in 2:nt
        ẋ = g(x, A)
        x .+= ẋ .* Δt
        Xt[i, :] = x
        sumxdots[i] = sum(ẋ)
    end
    return trange, Xt, sumxdots
end
trange, Xt, sumxdots = getReplicatorTrajectoryWithSumxDots()
# plot(trange, Xt)
plot!(trange, sumxdots)

function getReplicatorTrajectory()
    ntypes = 5
    A = rand(ntypes, ntypes)
    x = rand(ntypes); x ./= sum(x)
    g(x, A, t) = x .* (A * x .- x' * A * x)
    T = 100.00
    tspan = (0.0, T)
    prob = ODEProblem(g, x, tspan, A)
    sol = solve(prob, Tsit5())
    return sol 
end

Xt = getReplicatorTrajectory()
plot(Xt)