
function M0()
    λ(x) = 1.0 + x[1]
    α = 0.5
    σY = 0.01
    σX = 0.01
    D = 0.01

    function β(x, z, t)
        y0 = z[1]- z[2]
        μ = log(2) - α * y0
        A = exp(-((z[1] - μ)^2) / (2 * σY^2))/(2 * π * σY^2)^(1/2)
        B = (1 - erf((z[1] - μ) / (sqrt(2) * σY)))
        return  λ(x)*A/B
    end

    # symmetric division 
    h(z, x) = ([z[1] - log(2) .+ rand(Normal(0,0.1)),0],x .+ rand(Normal(0,σX),3))

    # OU process dynamics of x
    L(x, t, dt) =   -x*dt .+ sqrt(dt*D) .*rand(Normal(0,1),3)


    return λ, β, h, L
end

function M1()
    λ(x) = 1.0 + x[1]
    α = 0.5
    σY = 0.01
    D = 0.01
    σX = 0.01

    function β(x, z, t)
        y0 = z[1]- z[2]
        μ = log(2) - α * y0
        A = exp(-((z[1] - μ)^2) / (2 * σY^2))/(2 * π * σY^2)^(1/2)
        B = (1 - erf((z[1] - μ) / (sqrt(2) * σY)))
        return λ(x)*A/B
    end

    # symmetric division 
    h(z, x) = ([z[1] - log(2)  .+ rand(Normal(0,0.1)),0],x .+ rand(Normal(0,σX),3))

    # OU process dynamics of x
    L(x, t, dt) = -0.5 .*x*dt .+ sqrt(dt*D) .*rand(Normal(0,1),3)

    return λ, β, h, L
end


function M2()
    λ(x) = 1.0 + x[1]
    α = 0.5
    σY = 0.01
    σX = 0.01

    function β(x, z, t)
        y0 = z[1]- z[2]
        μ = log(2) - α * y0
        A = exp(-((z[1] - μ)^2) / (2 * σY^2))/(2 * π * σY^2)^(1/2)
        B = (1 - erf((z[1] - μ) / (sqrt(2) * σY)))
        return  λ(x)*A/B
    end

    # symmetric division 
    h(z, x) = ([z[1] - log(2)  .+ rand(Normal(0,0.1)),0],x .+ rand(Normal(0,σX),3))

    # OU process dynamics of x
    L(x, t, dt) = 0.0

    model = GrowthModel(λ, β, h, L)

    return  λ, β, h, L
end





