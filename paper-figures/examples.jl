


## OU Model with Division Noise and failing condition I
function M0()
    λ(x) = 1.0 + x[1]

    α = 0.5
    σY = 0.05
    D = 0.025 * log(2)  # D = 0.1 * γ² * ln(2) / 2 = 0.1 * 1² * ln(2) / 2
    γ = 1. 
    σX = sqrt(0.0001)
    # 

    function β(x, z, t)
        y0 = z[1]- z[2]
        μ = log(2) - α * y0
        A = exp(-((z[1] - μ)^2) / (2 * σY^2))/(2 * π * σY^2)^(1/2)
        B = (1 - erf((z[1] - μ) / (sqrt(2) * σY)))
        return  A/B
    end

    # symmetric division with Gamma noise (same mean=0, variance=σX²)
    function gamma_noise(σ)
        # Gamma(1, σ) has mean=σ, variance=σ²
        # Shift by -σ to get mean=0, variance=σ²
        return rand(Gamma(1, σ)) - σ *0.9
    end
    h(z, x) = ([z[1] - log(2) + rand(Normal(0,0.05)),0],x)

    # OU process dynamics of x (exact distribution)
    L(x, t, dt) = x .* (exp(-γ*dt) .-1) .+ sqrt(D * (1 - exp(-2*γ*dt))/(2*γ)) .* rand(Normal(0,1),3)

    return λ, β, h, L
end

## OU Model with no Division Noise
function M1()
    λ(x) = 1.0 + x[1]
    α = 0.5
    σY = 0.05
    γ = 1.
    D = 0.025 * log(2)  # D = 0.1 * γ² * ln(2) / 2 = 0.1 * 1² * ln(2) / 2

    function β(x, z, t)
        y0 = z[1]- z[2]
        μ = log(2) - α * y0
        A = exp(-((z[1] - μ)^2) / (2 * σY^2))/(2 * π * σY^2)^(1/2)
        B = (1 - erf((z[1] - μ) / (sqrt(2) * σY)))
        return λ(x)*A/B
    end

    # symmetric division 
    h(z, x) = ([z[1] - log(2) + rand(Normal(0,0.05)),0],x)

    # OU process dynamics of x (exact distribution)
    L(x, t, dt) = x .* (exp(-γ*dt) -1)  .+ sqrt(D * (1 - exp(-2*γ*dt))/(2*γ)) .* rand(Normal(0,1),3)

    return λ, β, h, L
end


function M2()
    λ(x) = 1.0 + x[1]
    α = 0.5
    σY = 0.05
    σX = sqrt(0.005)  # σX = sqrt(v_DN) = sqrt(0.1)


    function β(x, z, t)
        y0 = z[1]- z[2]
        μ = log(2) - α * y0
        A = exp(-((z[1] - μ)^2) / (2 * σY^2))/(2 * π * σY^2)^(1/2)
        B = (1 - erf((z[1] - μ) / (sqrt(2) * σY)))
        return  λ(x)*A/B
    end

    # symmetric division with Gamma noise (same mean=0, variance=σX²)
    # function gamma_noise(σ)
    #     # Gamma(1, σ) has mean=σ, variance=σ²
    #     # Shift by -σ to get mean=0, variance=σ²
    #     return rand(Gamma(1, σ)) - σ *0.9
    # end
    h(z, x) = ([z[1] - log(2) + rand(Normal(0,0.05)),0],rand(Normal(0,σX),3))

    # constant x in cell-cycle
    L(x, t, dt) = 0.0

    model = GrowthModel(λ, β, h, L)

    return  λ, β, h, L
end





