using RadiiPolynomial, IntervalArithmetic, LinearAlgebra 

include("path_of_the_folder/list_functions_SH.jl")
N = 150 
K = 15

N0 = 150
K0 = 24 ; S = Chebyshev(N0)⊗CosFourier(K0,1)

β1 = interval(1)  
β2 = -interval(2)
β0 = interval(7.1)
h = interval(3)
νinf = interval(8.3)

β = Sequence(Chebyshev(N0)⊗CosFourier(K0,interval(1)),interval.(zeros(dimension(S))))
b = Sequence(CosFourier(K0,interval(1)),interval.(zeros(dimension(CosFourier(K0,1)))))
b[1] = interval(0.01)
β[(0,0:K0)] = interval.(b[(:)])


# code to compute an approximate solution
# km = 80
# precis = 1e-14

# U0 = Sequence(S,mid.(coefficients(β)))
# U0 = interval.(newton_NS(mid.(U0),N0,K0,km,precis,mid.(β)))
# U0 = Sequence(Chebyshev(N0)⊗CosFourier(K0,interval(1)), coefficients(U0))

U0 = load("U0_SH.jld2","U0") # load the approximate solution
r = proof_integration_SH_sinus_Z1_new(U0,h,β,N0,K0,N,K)



######################## PLOTTING PART #################################################
using MATLAB

function PlotCoeffs2D(U0)
    #U0 is a sequence in 2D
    #a,b,c,d are the endpoints of the interval [a,b] × [c,d]

    dx = 2/1001
    dy = 2π/1001
    y1 = -1:dx:1
    y2 = 0:dy:2π
    m=length(y1)
    n=length(y2)
    U = zeros(m,n)
    for b₁ = 1:m
        for b₂ = 1:n
            U[b₁,b₂] = U0(y1[b₁],y2[b₂])
        end
    end
    
    y1 = mid(h)/2*(y1.+1)
    mat"
    h = surf($y2,$y1,$U)
    set(h,'LineStyle','none')"
end

U0 = Sequence(Chebyshev(N0)⊗CosFourier(K0,1), mid.(coefficients(U0)))
PlotCoeffs2D(mid.(U0))

