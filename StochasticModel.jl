begin
    using Random
    using Distributions
    using DifferentialEquations
    using DataFrames
    using Plots
end

struct ZeroInflatedDistribution{T<:Distributions.VariateForm, 
    S<:Distributions.ValueSupport} <:Distributions.Sampleable{T,S}
    p::Real
    d::Distributions.Sampleable{T,S}
end

Base.rand(d::ZeroInflatedDistribution) = Base.rand(Bernoulli(d.p)) == 1 ? 0.0 : Base.rand(d.d)

ZeroInflatedGamma(p, α, β) = ZeroInflatedDistribution(p, Gamma(α, β))
rand(ZeroInflatedGamma(0.1, 2.0, 1.0))


# Parameters 
β = [2.35 2.3; 2.0 2.5]  # Transmission rate
γ = [0.75, 0.5]   # Recovery rate


cost_infection_parms = [0.1, 1.3, 2.3]
cost_vaccination_parms = [0.1, 0.3, 1.3]


# Initial conditions
S_1 = 3999
S_2 = 6999
I_1 = 1
I_2 = 1
R_1 = 0
R_2 = 0

InitialSizes = [S_1 S_2 I_1 I_2 R_1 R_2]
N = sum(InitialSizes) 

# Gillespie algorithm for the model
function GillespieSimulation(β, γ, cost_infection_parms, cost_vaccination_parms, InitialSizes)
    times = Float64[]
    S1_values = Int[]
    S2_values = Int[]
    I1_values = Int[]
    I2_values = Int[]
    InfectionCosts = Float64[]
    # VaccinationCosts = Float64[]
    R1_values = Int[]
    R2_values = Int[]
    event_types = Int[]

    # Initial conditions
    t = 0.0
    S_1, S_2, I_1, I_2, R_1, R_2 = InitialSizes

    # Store results
    push!(times, t)
    push!(S1_values, S_1)
    push!(S2_values, S_2)
    push!(I1_values, I_1)
    push!(I2_values, I_2)
    push!(R1_values, R_1)
    push!(R2_values, R_2)
    push!(InfectionCosts, 0.0)
    # push!(VaccinationCosts, 0.0)
    push!(event_types, 0)


    while I_1 > 0 || I_2 > 0
        # Calculate rates
        rate_infection_1_1 = β[1, 1] * S_1 * I_1 / N
        rate_infection_1_2 = β[1, 2] * S_1 * I_2 / N
        rate_infection_2_1 = β[2, 1] * S_2 * I_1 / N
        rate_infection_2_2 = β[2, 2] * S_2 * I_2 / N
        rate_recovery_1 = γ[1] * I_1 
        rate_recovery_2 = γ[2] * I_2
        rate_total = rate_infection_1_1 + rate_infection_1_2 + rate_infection_2_1 + rate_infection_2_2 + rate_recovery_1 + rate_recovery_2

        # Time to next event
        τ = -log(rand()) / rate_total

        # Update time
        t += τ

        # Update species
        r = rand()
        if (r < rate_infection_1_1 / rate_total) & (S_1 > 0) & (I_1 > 0)
            S_1 -= 1
            I_1 += 1
            cost = rand(ZeroInflatedGamma(cost_infection_parms[1], cost_infection_parms[2], cost_infection_parms[3]))
            event_type = 1
        elseif (r < (rate_infection_1_1 + rate_infection_1_2) / rate_total) & (S_1 > 0) & (I_2 > 0)
            S_1 -= 1
            I_2 += 1
            cost = rand(ZeroInflatedGamma(cost_infection_parms[1], cost_infection_parms[2], cost_infection_parms[3]))
            event_type = 2
        elseif (r < (rate_infection_1_1 + rate_infection_1_2 + rate_infection_2_1) / rate_total) & (S_2 > 0) & (I_1 > 0)
            S_2 -= 1
            I_1 += 1
            cost = rand(ZeroInflatedGamma(cost_infection_parms[1], cost_infection_parms[2], cost_infection_parms[3]))
            event_type = 3
        elseif (r < (rate_infection_1_1 + rate_infection_1_2 + rate_infection_2_1 + rate_infection_2_2) / rate_total) & (S_2 > 0) & (I_2 > 0)
            S_2 -= 1
            I_2 += 1
            cost = rand(ZeroInflatedGamma(cost_infection_parms[1], cost_infection_parms[2], cost_infection_parms[3]))
            event_type = 4
        elseif (r < (rate_infection_1_1 + rate_infection_1_2 + rate_infection_2_1 + rate_infection_2_2 + rate_recovery_1) / rate_total) & (I_1 > 0)
            I_1 -= 1
            R_1 += 1
            cost = 0
            event_type = 5
        else (r < (rate_infection_1_1 + rate_infection_1_2 + rate_infection_2_1 + rate_infection_2_2 + rate_recovery_1) / rate_total) & (I_1 > 0)
            I_2 -= 1
            R_2 += 1
            cost = 0
            event_type = 6
        end

        # Store results
        push!(times, t)
        push!(S1_values, S_1)
        push!(S2_values, S_2)
        push!(I1_values, I_1)
        push!(I2_values, I_2)
        push!(R1_values, R_1)
        push!(R2_values, R_2)
        push!(InfectionCosts, cost)
        push!(event_types, event_type)
    end

    # Return results
    return DataFrame(times=times, S1=S1_values, S2=S2_values, I1=I1_values, I2=I2_values, R1=R1_values, R2=R2_values, InfectionCosts=InfectionCosts, event_types=event_types)

end

GillespieTrajectory = GillespieSimulation(β, γ, cost_infection_parms, cost_vaccination_parms, InitialSizes)


fig = plot(GillespieTrajectory[!, :times], GillespieTrajectory[!, :S1], color="lightgrey", linewidth=2, primary=false)
plot!(GillespieTrajectory[!, :times], GillespieTrajectory[!, :S2], color="grey", linewidth=2, primary=false)
plot!(GillespieTrajectory[!, :times], GillespieTrajectory[!, :I1], color="pink", linewidth=2, primary=false)
plot!(GillespieTrajectory[!, :times], GillespieTrajectory[!, :I2], color="red", linewidth=2, primary=false)
plot!(GillespieTrajectory[!, :times], GillespieTrajectory[!, :R1], color="lightgreen", linewidth=2, primary=false)
plot!(GillespieTrajectory[!, :times], GillespieTrajectory[!, :R2], color="green", linewidth=2, primary=false)


# Plot infection costs
histogram(GillespieTrajectory[!, :InfectionCosts], color="blue")



plot(GillespieTrajectory[!, :times], GillespieTrajectory[!, :I1], color="pink", linewidth=2, primary=false)
plot!(GillespieTrajectory[!, :times], GillespieTrajectory[!, :I2], color="red", linewidth=2, primary=false)
