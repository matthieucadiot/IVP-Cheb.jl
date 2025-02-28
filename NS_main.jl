using RadiiPolynomial, IntervalArithmetic, LinearAlgebra 

include("path_of_the_folder/list_functions_NS.jl")

N = 3
K = 10

N0 = 10
K0 = 15 ; S = Chebyshev(N0)⊗SinFourier(K0,1)⊗SinFourier(K0,1)

ν = interval(0.2)
h = interval(0.03)

β = Sequence(Chebyshev(N0)⊗SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1)),interval.(zeros(dimension(S))))
b = Sequence(SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1)),interval.(zeros(dimension(SinFourier(K0,1)⊗SinFourier(K0,1)))))
b[(2,1)] = interval(0.25)
b[(1,2)] = -interval(0.25)
b[(3,1)] = -interval(0.1)
b[(1,3)] = interval(0.1)
b[(3,2)] = -interval(0.15)
b[(2,3)] = interval(0.15)
β[(0,1:K0,1:K0)] = b[(:,:)]

km = 10
precis = 1e-14

U0 = Sequence(S,mid.(coefficients(β)))
U0 = interval.(newton_NS(mid.(U0),N0,K0,km,precis,mid.(β)))
U0 = Sequence(Chebyshev(N0)⊗SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1)), coefficients(U0))

r = proof_integration_NS_sinus_Z1_new(U0,h,ν,β,N0,K0,N,K,interval(0))

b = U0((1,nothing,nothing))
for k1 = 1:K0
    for k2=1:K0
        β[(0,k1,k2)] = b[(0,k1,k2)]
    end
end
U0 = Sequence(S,mid.(coefficients(β)))
U0 = interval.(newton_NS(mid.(U0),N0,K0,km,precis,mid.(β)))
U0 = Sequence(Chebyshev(N0)⊗SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1)), coefficients(U0))


for j = 1:73
    global r,U0,b
    r = proof_integration_NS_sinus_Z1_new(U0,h,ν,β,N0,K0,N,K,r)
    display(r)
    b = U0((1,nothing,nothing))
    for k1 = 1:K0
        for k2=1:K0
            β[(0,k1,k2)] = b[(0,k1,k2)]
        end
    end
    display(j)
        if sup(norm(b,1))<inf(interval(0.116))
            display("global existence is reached !")
            display(j)
            display(norm(b,1))
            display(r)
        end
    U0 = Sequence(S,mid.(coefficients(β)))
    U0 = interval.(newton_NS(mid.(U0),N0,K0,km,precis,mid.(β)))
    U0 = Sequence(Chebyshev(N0)⊗SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1)), coefficients(U0))
end


r = proof_integration_NS_sinus_Z1_new(U0,interval(0.72),ν,β,N0,K0,N,K,r)
if sup(norm(b,1))<inf(interval(0.115))
    display("global existence is reached !")
    display(norm(b,1))
    display(r)
end


######################## PLOTTING PART #################################################

using MATLAB

function PlotCoeffs3D(U0,t)
    #U0 is a sequence in 2D
    #a,b,c,d are the endpoints of the interval [a,b] × [c,d]

    dx = 2π/1001
    dy = 2π/1001
    y1 = 0:dx:2π
    y2 = 0:dy:2π
    m=length(y1)
    n=length(y2)
    U = zeros(m,n)
    for b₁ = 1:m
        for b₂ = 1:n
            U[b₁,b₂] = U0(t,y1[b₁],y2[b₂])
        end
    end
    U = real.(U)
    mat"
    h = surf($y2,$y1,$U)
    set(h,'LineStyle','none')"
end

U0 = Sequence(S,mid.(coefficients(U0)))
 PlotCoeffs3D(U0,1)


