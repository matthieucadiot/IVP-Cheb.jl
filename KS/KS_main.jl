

using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, TickTock
include("my_path_to_folder/list_functions_KS.jl")

IntervalArithmetic.matmul_mode() = IntervalArithmetic.MatMulMode{:fast}()

# values for the size of the sequence
N0 = 80
K0 = 31
# values for the size of the operators
N = 180
K = 20

# parameters of the equation and of the IVP
ν = -interval(2)
h = interval(2.244333563761175)
α = interval(0.127)
νinf = interval(0.25)/α # νinf is chosen so that μk is always positive (cf. equation (8))

S0 = Chebyshev(N0)⊗SinFourier(K0,interval(1)) # definition of the spaces
S = Chebyshev(N)⊗SinFourier(K,interval(1))

β = Sequence(S0,vec(interval.(zeros((N0+1)*(K0)))))
# the initial condition is a point on the periodic orbit in [40]
b = Sequence(SinFourier(K0,1),interval.(zeros(K0)))
b0 = [
2.012364292932942e-01;
1.289993211670137e+00;
2.011499954939851e-01;
-3.778743081553121e-01;
-4.229439951263216e-02;
4.316415107347703e-02;
6.937968915698671e-03;
-4.156922054455959e-03;
-7.942420514021223e-04;
3.316645095010771e-04;
7.936884299000886e-05;
-2.391471591898029e-05;
-7.079686365878430e-06;
1.568946982933196e-06;
5.845069061129325e-07;
-9.427015521193508e-08;
-4.537920144230167e-08;
5.108267040182597e-09;
3.353536537928777e-09;
-2.381786166952790e-10;
-2.377246547769026e-10;
8.056817619592433e-12;
1.625293162000136e-11;
1.164121054373202e-14;
-1.075623410628802e-12;
-3.685233168519958e-14;
6.906797853087279e-14;
4.698824785510242e-15;
-4.308234657272735e-15;
-4.452347717777985e-16;
2.610688409363968e-16]
b[1:K0] = b0[1:K0]
for k = 1:K0
    β[(0,k)] = b[k] # construction of the initial date in Fourier-Cheb
end

# the commented code below allows to compute an approximate solution thanks to Newton's method
# ℒ0 = LinearOperator(S0,S0,interval.(zeros(dimension(S0),dimension(S0))))
# Λ0 = LinearOperator(S0,S0,interval.(zeros(dimension(S0),dimension(S0))))
# for k = 1:K0
#     μk = h*interval(0.5)*(-interval(k)^2+ α*interval(k)^4 + νinf)
#     Λ0[(:,k),(:,k)]  = conversion(N0+1,N0+1)
#     ℒ0[(:,k),(:,k)] = complete_op(N0+1,N0+1,μk)
#     β[(0,k)] = b[k]
# end
# km = 10
# precis = 2e-14
# U0 = Sequence(Chebyshev(N0)⊗SinFourier(K0,1),mid.(coefficients(β)))
# U0 = int_cheb(U0,Chebyshev(N0)⊗CosFourier(K0,1),mid(ν),mid(νinf),mid.(ℒ0),mid.(Λ0),mid(h),mid.(β))
# U0 = Sequence(S0,interval.(coefficients(U0)))
# ℒ0 = Nothing ; Λ0 = Nothing

# construction of operators of interest
ℒ = LinearOperator(S,S,interval.(zeros(dimension(S),dimension(S))))
Λ = LinearOperator(S,S,interval.(zeros(dimension(S),dimension(S))))
Ω = Sequence(S,interval.(2*ones(dimension(S)))) # this sequences allows to define the weights in the norm
Ωinv = Sequence(S,interval.(interval(0.5)*interval.(ones(dimension(S)))))
for k = 1:K
    μk = h*interval(0.5)*(-interval(k)^2+ α*interval(k)^4 + νinf)
    Λ[(:,k),(:,k)]  = conversion(N+1,N+1)
    ℒ[(:,k),(:,k)] = complete_op(N+1,N+1,μk)
    Ω[(0,k)] = interval(1)
    Ωinv[(0,k)] = interval(1)
end

# computer assisted proof of existence. The value of r is the ball in which the solution is proven.
U0 = load("U0_KS.jld2","U0") # load the approximate solution
r = proof_integration_KS_sinus_Z1_new(U0,h,ν,β,N0,K0,N,K)

# plot of the solution
using MATLAB
function PlotCoeffs2D(U0)
    dx = 2/1001
    dy = 2π/1001
    y1 = -1:dx:1
    y2 = -π:dy:π
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
U0 = Sequence(Chebyshev(N0)⊗SinFourier(K0,1), mid.(coefficients(U0)))
 PlotCoeffs2D(U0)