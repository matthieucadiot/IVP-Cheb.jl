
using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, JLD2

IntervalArithmetic.matmul_mode() = IntervalArithmetic.MatMulMode{:fast}()


function tridiag_int(N,μ::Float64)
   
    b = 1.0*μ*ones(N-1)
    a = 2.0*collect(1:N)   
    return Tridiagonal(b,a,-b)
end

function tridiag_int(N,M,μ::Interval{Float64})
   if N==M
        b = μ*interval.(ones(N-1))
        a = interval(2)*interval.(collect(1:N))   
        return Tridiagonal(b,a,-b)
   elseif M>N
        T = interval(zeros(N,M))
        b = μ*interval.(ones(N-1))
        a = interval(2)*interval.(collect(1:N))   
        T[1:N,1:N] = Tridiagonal(b,a,-b)
        T[N,N+1] = -μ
        return T
   else
        display("need to choose the number of columns bigger")
   end
end


function complete_op(N,μ::Float64)
    T = zeros(N,N)
    T[1,1] = 1.0
    T[1,2:N] = 2.0*((-1).^(1:N-1))'
    T[2,1] = μ
    T[2:end,2:end] = Matrix(tridiag_int(N-1,μ))
    return T
end

function complete_op(N,M,μ::Interval{Float64})
    T = interval.(zeros(N,M))
    T[1,1] = interval(1.0)
    T[1,2:M] = interval(2.0)*interval.(((-1).^(1:M-1))')
    T[2,1] = μ
    T[2:end,2:end] = Matrix(tridiag_int(N-1,M-1,μ))
    return T
end

function conversion(N,M)
    T = interval.(zeros(N,M))
    for i=1:N
        for j=1:N
            if j==i+1
                T[i,j] = interval(1)
            elseif i==j+1
                T[i,j] =- interval(1) 
            end
        end
    end
    T[1,2]= interval(0)
    if N==M
        return T
    else
        T[N,N+1] = interval(1)
        return T
    end
end



function lambda_ext(U,K,N)
    S = Chebyshev(N+1)⊗CosFourier(K,1)
    F = Sequence(S,interval.(zeros((N+2)*(K+1))))
    for k ∈ 0:K
            for n ∈ 1:N-1
                F[(n,k)] = -U[(n-1,k)]+U[(n+1,k)]
            end

            F[(N,k)] = -U[(N-1,k)]
            F[(N+1,k)] = -U[(N,k)]
        
    end
    return F
end


function F(U,ℒ,Λ,h,β)
    return ℒ*U - mid(h)/2*Λ*project(U^3,space(U)) - β
end

function DF(U,ℒ,Λ)
    return ℒ - 3*mid(h)/2*Λ*project(Multiplication(U^2),space(U),space(U),Float64)
end


function newton_NS(U,N0,K0,km,precis,β)
S = Chebyshev(N0)⊗CosFourier(K0,1)
ℒ = LinearOperator(S,S,(zeros(dimension(S),dimension(S))))
Λ = LinearOperator(S,S,(zeros(dimension(S),dimension(S))))
for k = 0:K0
            μk = h*0.5*(mid(β1)*k^4 + mid(β2)*k^2 - mid(β0))
            Λ[(:,k),(:,k)]  = mid.(conversion(N0+1,N0+1))
            ℒ[(:,k),(:,k)] = mid.(complete_op(N0+1,N0+1,interval(μk)))
end
    F1 = F(U,ℒ,Λ,mid(h),mid.(β))
    nf = norm(coefficients(F1),Inf)
    p=1
    display(nf)
    while (p<km)&&(nf>precis)
        DF1 = DF(U,ℒ,Λ)
        U = U- DF1\F1
        F1 = F(U,ℒ,Λ,mid(h),mid.(β))
        nf = norm(coefficients(F1),Inf)
        display(nf)
        p = p+1
    end
    return U
end

function int_cheb(U,ℒ,Λ,h,β)
    U, _ = newton(U; maxiter = 20, tol = 1.0e-14) do U
        return F(U,ℒ,Λ,h,β), DF(U,ℒ,Λ)
     end
    return U 
end

function function_g(x,μ)
    return μ*(x*log(sqrt(x^2 + interval(1)) + x) -sqrt(x^2 + interval(1)) + interval(1))
end

function norm_first_column(N,μ)
    T = tridiag_int(N,N,μ)
    e0 = interval.(zeros(N)) ; e0[1] = interval(1)
    x = interval.(mid.(T)\mid.(e0))
    err = first_column_residual(N,μ)
    return interval(inf(norm(x,1) - norm(T*x-e0,1)-err),sup(norm(x,1) + norm(T*x-e0,1)+err))
end

function first_column_residual(N,μ)
    N = interval(N)
    if (sup(N/μ) <= 0.5)
        return exp(-function_g(N/μ,μ))/sqrt(interval(4)+μ^2) + μ/sqrt(interval(4)+μ^2)*(exp(-interval(0.49)*N^2/μ)*interval(0.5)*minimum([interval(1) sqrt(interval(π)/(interval(0.49)*μ))]) + interval(0.6)*exp(-interval(0.115)*μ) +interval(2)*exp(-0.55*μ)/μ )
    elseif (0.5<inf(N/μ))&&(sup(N/μ)<=1.1)
        return exp(-function_g(N/μ,μ))/sqrt(interval(4)+μ^2) + μ/sqrt(interval(4)+μ^2)*(exp(-interval(0.46)*N^2/μ)*minimum([interval(0.6) sqrt(interval(π)/(interval(4)*interval(0.46)*μ))])  +interval(2)*exp(-interval(0.55)*μ)/μ )
    else
        return exp(-function_g(N/μ,μ))/sqrt(interval(4)+μ^2)  + interval(2)*exp(-interval(0.55)*interval(N))/sqrt(interval(4)+μ^2)
    end
end

function fct_μk(k)
    return  h*interval(0.5)*(β1*interval(k)^4 + β2*interval(k)^2 - β0  + νinf)
end

function component1_Z∞1_Z∞2(V)
    Z∞1 = interval(0) ; Z∞2 = interval(0)
    Vc = Sequence(Chebyshev(N)⊗CosFourier(3K,interval(1)),interval.(zeros((N+1)*(3*K+1))))
    Vc[(:,0:K)] = V[(:,:)]
    d = interval.(2*ones(N+1)) ; d[1] = interval(1)
    for k1 = K+1:2*K
            S = interval(0)
            for k2 = K+1:K+k1
                        μk = fct_μk(k2) 
                        S += function_C0(μk)/(μk+interval(2))*norm(d.*Vc[(:,abs(k2-k1))] ,1)
            end
            Z∞1 = maximum([Z∞1 S])
    end
    for k1 = 1:K
            S = interval(0)
            for k2 =  K+1:K+k1
                    μk = fct_μk(k2) 
                    S += function_C0(μk)/(μk+interval(2))*norm(d.*Vc[(:,k2-k1)] ,1)
            end
            Z∞2 = maximum([Z∞2 S])
    end
    return Z∞1, Z∞2
end

function component2_Z∞1_Z∞2_ZN_ZA0(V,D_C1,D_A,N00)
    Z∞1 = interval(0) ; Z∞2 = interval(0); ZN = interval(0) ;  ZA0 = interval(0)
    Vc = Sequence(Chebyshev(N)⊗CosFourier(2K,interval(1)),interval.(zeros((N+1)*(2*K+1))))
    Vc[(:,0:K)] = V[(:,:)]
    d = interval.(2*ones(N+1)) ; d[1] = interval(1)
    
    # Z∞1 can be improved by using C1 instead of C0. Aut if μk increases quick enough, that might be enough
    for k1 =  K+1:2*K
        S = interval(0)
        for k2 = 0:K
            μk = fct_μk(k2) 
            S += function_C0(μk)/(μk+interval(2))*norm(d.*Vc[(:,k1-k2)],1)
        end
        Z∞1 = maximum([Z∞1 S])
    end

    Vc = Sequence(Chebyshev(N00+2N)⊗CosFourier(2K,interval(1)),interval.(zeros(((N00+2N)+1)*(2*K+1))))
    Vc[(0:N,0:K)] = V[(:,:)]

    for k1 = 0:K
        for n1 =  N+1:N00+N
        S13 = interval(0)
        S15 = interval(0)
            for k2 = 0:K
                for n2 = 0:minimum([N+n1 N00])
                    S13 += D_C1[(n2,k2)]*(abs(Vc[(abs(n1-n2),abs(k1-k2))]) + abs(Vc[(abs(n1-n2),k1+k2)]))
                    S15 += D_A[(n2,k2)]*(abs(Vc[(abs(n1-n2),abs(k1-k2))])  + abs(Vc[(abs(n1-n2),k1+k2)]))
                end
            end
            Z∞2 = maximum([Z∞2 S13])
            ZA0 = maximum([ZA0 S15])
            
        end

        for n1 = 0:N
            S14 = interval(0)
                for k2 = 0:K
                    for n2 = N+1:2N 
                        S14 += D_C1[(n2,k2)]*(abs(Vc[(n2-n1,abs(k1-k2))])+abs(Vc[(n2-n1,k1+k2)]))
                    end
                end
                ZN = maximum([ZN S14])
            end
    end

    Vc = Sequence(Chebyshev(N00+2N)⊗CosFourier(2K,interval(1)),interval.(zeros(((N00+2N)+1)*(2K+1))))
    Vc[(0:N,0:K)] = V[(:,:)]

    for k1 = K+1:2*K
        for n1 = 0:N00
        S15 = interval(0)
            for k2 = 0:K
                S15 += D_A[(0,k2)]*abs(Vc[(n1,k1-k2)])
                for n2 = 1:N
                    S15 += D_A[(n2,k2)]*(abs(Vc[(abs(n1-n2),k1-k2)])+abs(Vc[(n1+n2,k1-k2)]))
                end
                for n2 = N+1:N+n1
                    S15 += D_A[(n2,k2)]*abs(Vc[(abs(n1-n2),k1-k2)])
                end
            end
            ZA0 = maximum([ZA0 S15])
        end
    end

    return Z∞1, Z∞2, ZN, ZA0
end

function Tinv_V(V,μ,N)

    T = Matrix(tridiag_int(N,N,μ))
    Tinv = interval.(inv(mid.(T)))
    nT = opnorm(I - Tinv*T,1)

    α = interval(1)/interval(inf(interval(N)+sqrt(interval(N)^2 + μ^2)),sup(interval(N+1) + sqrt(interval((N+1)^2) + μ^2)))
    U = Tinv*V
    col_1 = Tinv[:,1]

    U = U - μ^2*U[N]*Tinv[:,N]*interval(inf(α)/(1+sup(μ)^2*inf(α)*sup(Tinv[N,N])), sup(α)/(1+inf(μ)^2*sup(α)*inf(Tinv[N,N]))) 
    return interval.(inf.(U .- abs.(U)*nT/(interval(1)-nT)), sup.(U .+ abs.(U)*nT/(interval(1)-nT))), abs(Tinv[N,N])-nT*abs(Tinv[N,N]), interval.(inf.(col_1- nT*abs.(col_1)),sup.(col_1+nT*abs.(col_1)))
end

function Linv_V(V,μ,N)
    # Computation of Tinv*V
    U, err_end, col_1 =  Tinv_V(V[2:end],μ,N)

    N00 = Int(floor(mid.(μ)/5)) + 100
    norm_col_1 =   norm_first_column(N00+1,μ)

    α = interval(1)/interval(inf(interval(N)+sqrt(interval(N)^2 + μ^2)),sup(interval(N+1) + sqrt(interval((N+1)^2) + μ^2)))
    err_tail = μ*abs(U[end])/(interval(1)+μ^2*α*err_end)*norm_col_1*sqrt(interval(4)+μ^2)/(interval(N) + sqrt(interval(N)^2+μ^2))
    
    return [(-interval(2)*norm_col_1*V[1] - dot(interval(2)*interval.(((-1).^(1:N))),U))/(interval(1)+interval(2)*μ*norm_col_1) ;interval(1)/(interval(1)+interval(2)*μ*norm_col_1)*(-V[1] + dot(interval(2)*μ*interval.(((-1).^(1:N))),U) )*col_1+U],err_tail 
end

# using the appendix, the inverse of T can be computed by computed the inverse of T(N) and by adding some error coming from the last column of the inverse. The function below takes care of this error
function err_inverse(Tinv,N,μ)
    S = interval(zeros(N,N))
    α = interval(1)/interval(inf(interval(N)+sqrt(interval(N)^2 + μ^2)),sup(interval(N+1) + sqrt(interval((N+1)^2) + μ^2)))
    for i = 1:N
        for n = 1:N 
            S[n,i] = -μ^2*α*Tinv[N,i]*Tinv[n,N]/(interval(1) + μ^2*α*Tinv[N,N]) 
        end
    end
return S
end

function Linv_Lamb(μ,N)

    # V as become Lambda_0*v
    T = Matrix(tridiag_int(N,N,μ))
    Tinv = interval.(inv(mid.(T)))
    nT = opnorm(I - Tinv*T,1)
    
    # Product Tinv \Lambda
    Tlamb = (Tinv*conversion(N+2,N+2)[2:N+1,2:N+2])[1:N,1:N] 

    # error from the approximate inverse
    Tlamb = interval.(inf.(Tlamb-abs.(Tlamb)*nT/(interval(1)-nT)), sup.(Tlamb+abs.(Tlamb)*nT/(interval(1)-nT)))

    # formula above equation (96)
    Tlamb = Tlamb + err_inverse(Tlamb,N,μ)

    # computation of the norm of the first column
    N00 = Int(floor(mid.(μ)/5)) + 100
    norm_col_1 = norm_first_column(N00+1,μ)
    theta = (interval.(2*(-1).^(1:N)))'

    # propagation of error from the term θ0^T Tinv Λ 
    T_err = [exp(function_g(interval(i-2)/μ,μ) - function_g(interval(N-2)/μ,μ))/(interval(N-1) + sqrt(interval(N-1)^2 + μ^2))+exp(function_g(interval(i-2)/μ,μ) - function_g(interval(N-1)/μ,μ))/(interval(N) + sqrt(interval(N)^2 + μ^2)) for i∈1:N]' 

    # error propagation in each column coming from the first column of Tinv
    err_col = function_C0(μ)/(interval(2)+μ)*μ/(interval(1)+interval(2)*norm_col_1*μ)*first_column_residual(N+2,μ)*ones(N+1) ; err_col[1] =  interval(1)/(interval(1)+interval(2)*norm_col_1*μ)*first_column_residual(N+2,μ)
    err_col[2:end] = err_col[2:end] + T_err'

    T_err = interval.(-sup.(T_err),sup.(T_err))

     return [-interval(2)*norm_col_1/(interval(1)+interval(2)*norm_col_1*μ) -theta*Tlamb/(interval(1)+interval(2)*norm_col_1*μ)+T_err;
             -Tinv[:,1]/(interval(1)+interval(2)*norm_col_1*μ)  (I + μ/(interval(1)+interval(2)*norm_col_1*μ)*Tinv[:,1]*theta)*Tlamb+μ/(interval(1)+interval(2)*norm_col_1*μ)*Tinv[:,1]*T_err], err_col
end


function function_C0(μ)
    
    if sup(μ)<=1
        C0_vec = load("C0_1.jld2","C0_vec")
        i = Int(floor(mid(μ)/0.01))+1 
        return C0_vec[i]*(μ+interval(2))/(interval(2) + interval(i*0.01) - abs(μ-interval(i*0.01))*C0_vec[i])
    elseif (inf(μ)>1)&&(sup(μ)<=20)
        C0_vec = load("C0_2.jld2","C0_vec")
        i = Int(floor((mid(μ)-1)/0.1))+1 ; μ0 = interval(1) + interval(i*0.1)
        return C0_vec[i]*(μ+interval(2))/(interval(2) + interval(μ0) - abs(μ-interval(μ0))*C0_vec[i])
    elseif (inf(μ)>20)&&(sup(μ)<=100)
        C0_vec = load("C0_3.jld2","C0_vec")
        i = Int(floor((mid(μ)-20)/0.5))+1 ; μ0 = interval(20) + interval(i*0.5)
        return C0_vec[i]*(μ+interval(2))/(interval(2) + interval(μ0) - abs(μ-interval(μ0))*C0_vec[i])
    elseif (inf(μ)>100)&&(sup(μ)<=200)
        C0_vec = load("C0_4.jld2","C0_vec")
        i = Int(floor((mid(μ)-100)))+1 ; μ0 = interval(100) + interval(i)
        return C0_vec[i]*(μ+interval(2))/(interval(2) + interval(μ0) - abs(μ-interval(μ0))*C0_vec[i])
    elseif (inf(μ)>200)&&(sup(μ)<305)
        C0_vec = load("C0_5.jld2","C0_vec")
        i = Int(floor((mid(μ)-200)/5))+1 ; μ0 = interval(200) + interval(i*5)
        return C0_vec[i]*(μ+interval(2))/(interval(2) + interval(μ0) - abs(μ-interval(μ0))*C0_vec[i])
    elseif (inf(μ)>=305)&&(sup(μ)<=400)
        μ1 = 305 ; dμ = 1
        C0_vec = load("C0_6.jld2","C0_vec")
        i = Int(floor((mid(μ)-μ1)/dμ))+1 ; μ0 = interval(μ1) + interval(i*dμ)
        return C0_vec[i]*(μ+interval(2))/(interval(2) + interval(μ0) - abs(μ-interval(μ0))*C0_vec[i])
    elseif (inf(μ)>400)&&(sup(μ)<410)
        μ0 = interval(400) 
        C0_vec = load("C0_6.jld2","C0_vec")
        return C0_vec[end]*(μ+interval(2))/(interval(2) + interval(μ0) - abs(μ-interval(μ0))*C0_vec[end])   
    elseif (inf(μ)>=410)&&(sup(μ)<700)
        μ1 = 410 ; dμ = 10
        C0_vec = load("C0_7.jld2","C0_vec")
        i = Int(floor((mid(μ)-μ1)/dμ))+1 ; μ0 = interval(μ1) + interval(i*dμ)
        return C0_vec[i]*(μ+interval(2))/(interval(2) + interval(μ0) - abs(μ-interval(μ0))*C0_vec[i])
    elseif (inf(μ)>=700)&&(sup(μ)<=710)
        μ0 = interval(700) 
        C0_vec = load("C0_7.jld2","C0_vec")
        return C0_vec[end]*(μ+interval(2))/(interval(2) + interval(μ0) - abs(μ-interval(μ0))*C0_vec[end])  
    elseif (inf(μ)>710)&&(sup(μ)<=1000)
        μ1 = 710 ; dμ = 15
        C0_vec = load("C0_8.jld2","C0_vec")
        i = Int(floor((mid(μ)-μ1)/dμ))+1 ; μ0 = interval(μ1) + interval(i*dμ)
        return C0_vec[i]*(μ+interval(2))/(interval(2) + interval(μ0) - abs(μ-interval(μ0))*C0_vec[i])
    elseif (inf(μ)>1000)&&(sup(μ)<=2000)
        μ1 = 1000 ; dμ = 50
        C0_vec = load("C0_9.jld2","C0_vec")
        i = Int(floor((mid(μ)-μ1)/dμ))+1 ; μ0 = interval(μ1) + interval(i*dμ)
        return C0_vec[i]*(μ+interval(2))/(interval(2) + interval(μ0) - abs(μ-interval(μ0))*C0_vec[i])
    elseif (inf(μ)>2000)&&(sup(μ)<=3000)
        μ1 = 2000 ; dμ = 100
        C0_vec = load("C0_10.jld2","C0_vec")
        i = Int(floor((mid(μ)-μ1)/dμ))+1 ; μ0 = interval(μ1) + interval(i*dμ)
        return C0_vec[i]*(μ+interval(2))/(interval(2) + interval(μ0) - abs(μ-interval(μ0))*C0_vec[i])
    elseif (inf(μ)>3000)&&(sup(μ)<=10000)
        return interval(3.5)
    elseif (inf(μ)>10000)&&(sup(μ)<=inf((interval(10002)+interval(10000)*interval(3.5))/interval(3.5)))
        μ0 = interval(3000)
        return minimum([interval(15.71) interval(3.5)*(μ+interval(2))/(interval(2) + μ0 - abs(μ-μ0)*interval(3.5))])
    else
        return interval(15.71)
    end
end

function function_C1(μ)
    if sup(μ)<=1
        C1_vec = load("C1_1.jld2","C1_vec")
        C0_vec = load("C0_1.jld2","C0_vec")
        i = Int(floor(mid(μ)/0.01))+1 ; μ0 = interval(i*0.01)
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[i]/(interval(1) - abs(μ-μ0)*C0_vec[i]/(μ0+interval(2)))
    elseif (inf(μ)>1)&&(sup(μ)<=20)
        C1_vec = load("C1_2.jld2","C1_vec")
        C0_vec = load("C0_2.jld2","C0_vec")
        i = Int(floor((mid(μ)-1)/0.1))+1 ; μ0 = interval(1) + interval(i*0.1)
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[i]/(interval(1) - abs(μ-μ0)*C0_vec[i]/(μ0+interval(2)))
    elseif (inf(μ)>20)&&(sup(μ)<=100)
        C1_vec = load("C1_3.jld2","C1_vec")
        C0_vec = load("C0_3.jld2","C0_vec")
        i = Int(floor((mid(μ)-20)/0.5))+1 ; μ0 = interval(20) + interval(i*0.5)
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[i]/(interval(1) - abs(μ-μ0)*C0_vec[i]/(μ0+interval(2)))
    elseif (inf(μ)>100)&&(sup(μ)<=200)
        C1_vec = load("C1_4.jld2","C1_vec")
        C0_vec = load("C0_4.jld2","C0_vec")
        i = Int(floor((mid(μ)-100)))+1 ; μ0 = interval(100) + interval(i)
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[i]/(interval(1) - abs(μ-μ0)*C0_vec[i]/(μ0+interval(2)))
    elseif (inf(μ)>200)&&(sup(μ)<305)
        C1_vec = load("C1_5.jld2","C1_vec")
        C0_vec = load("C0_5.jld2","C0_vec")
        i = Int(floor((mid(μ)-200)/5))+1 ; μ0 = interval(200) + interval(i*5) ; 
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[i]/(interval(1) - abs(μ-μ0)*C0_vec[i]/(μ0+interval(2)))
    elseif (inf(μ)>=305)&&(sup(μ)<=400)
        μ1 = 305 ; dμ = 1
        C1_vec = load("C1_6.jld2","C1_vec")
        C0_vec = load("C0_6.jld2","C0_vec")
        i = Int(floor((mid(μ)-μ1)/dμ))+1 ; μ0 = interval(μ1) + interval(i*dμ) ; 
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[i]/(interval(1) - abs(μ-μ0)*C0_vec[i]/(μ0+interval(2)))
    elseif (inf(μ)>400)&&(sup(μ)<410)
        C1_vec = load("C1_6.jld2","C1_vec")
        C0_vec = load("C0_6.jld2","C0_vec")
         μ0 = interval(400)
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[end]/(interval(1) - abs(μ-μ0)*C0_vec[end]/(μ0+interval(2)))
    elseif (inf(μ)>=410)&&(sup(μ)<700)
        μ1 = 410 ; dμ = 10
        C1_vec = load("C1_7.jld2","C1_vec")
        C0_vec = load("C0_7.jld2","C0_vec")
        i = Int(floor((mid(μ)-μ1)/dμ))+1 ; μ0 = interval(μ1) + interval(i*dμ) ; 
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[i]/(interval(1) - abs(μ-μ0)*C0_vec[i]/(μ0+interval(2)))
    elseif (inf(μ)>=700)&&(sup(μ)<=710)
        C1_vec = load("C1_7.jld2","C1_vec")
        C0_vec = load("C0_7.jld2","C0_vec")
         μ0 = interval(700)
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[end]/(interval(1) - abs(μ-μ0)*C0_vec[end]/(μ0+interval(2)))
    elseif (inf(μ)>710)&&(sup(μ)<=1000)
        μ1 = 710 ; dμ = 15
        C1_vec = load("C1_8.jld2","C1_vec")
        C0_vec = load("C0_8.jld2","C0_vec")
        i = Int(floor((mid(μ)-μ1)/dμ))+1 ; μ0 = interval(μ1) + interval(i*dμ) ; 
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[i]/(interval(1) - abs(μ-μ0)*C0_vec[i]/(μ0+interval(2)))
    elseif (inf(μ)>1000)&&(sup(μ)<=2000)
        μ1 = 1000 ; dμ = 50
        C1_vec = load("C1_9.jld2","C1_vec")
        C0_vec = load("C0_9.jld2","C0_vec")
        i = Int(floor((mid(μ)-μ1)/dμ))+1 ; μ0 = interval(μ1) + interval(i*dμ) ; 
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[i]/(interval(1) - abs(μ-μ0)*C0_vec[i]/(μ0+interval(2)))
    elseif (inf(μ)>2000)&&(sup(μ)<=3000)
        μ1 = 2000 ; dμ = 100
        C1_vec = load("C1_10.jld2","C1_vec")
        C0_vec = load("C0_10.jld2","C0_vec")
        i = Int(floor((mid(μ)-μ1)/dμ))+1 ; μ0 = interval(μ1) + interval(i*dμ) ; 
        return (interval(1)+μ)/(interval(1)+μ0)*C1_vec[i]/(interval(1) - abs(μ-μ0)*C0_vec[i]/(μ0+interval(2)))
    elseif (inf(μ)>3000)&&(sup(μ)<=10000)
        return interval(3.5)
    elseif (inf(μ)>10000)&&(sup(μ)<=inf((interval(10002)+interval(10000)*interval(3.5))/interval(3.5)))
        μ0 = interval(3000)
        return minimum([interval(15.71) (interval(1)+μ)/(interval(1)+μ0)*interval(3.5)/(interval(1) - abs(μ-μ0)*interval(3.5)/(μ0+interval(2)))])
    else
        return interval(15.71)
    end
end





function proof_integration_SH_sinus_Z1_new(U0,h,β,N0,K0,N,K)

S = Chebyshev(N)⊗CosFourier(K,interval(1))
ℒ = LinearOperator(S,S,interval.(zeros(dimension(S),dimension(S))))
Λ = LinearOperator(S,S,interval.(zeros(dimension(S),dimension(S))))
Ω = interval(2)*Sequence(S,interval.(ones(dimension(S))))

for k = 0:K
    μk = fct_μk(k)
     Λ[(:,k),(:,k)]  = conversion(N+1,N+1)
    ℒ[(:,k),(:,k)] = complete_op(N+1,N+1,μk)
    Ω[(0,k)] = interval(1)
end
    A = interval.(inv(mid.(DF(project(mid.(U0),S),mid.(ℒ),mid.(Λ))))*mid.(ℒ))
    norm_A = opnorm(A,1)
    display("A created")

    # Computation of the bound Z1
    N00 = maximum([3*N+2 300])
    D_C1 = Sequence(Chebyshev(N00+N)⊗CosFourier(K,interval(1)),interval.(zeros((N00+N+1)*(K+1))))
    D_A = Sequence(Chebyshev(N00+N)⊗CosFourier(K,interval(1)),interval.(zeros((N00+N+1)*(K+1))))

    LU0 = Sequence(Chebyshev(N0)⊗CosFourier(K0,interval(1)),interval.(zeros((N0+1)*(K0+1))))
    for k = 0:K0
        μk = h*interval(0.5)*(β1*interval(k)^4 + β2*interval(k)^2 - β0)
        LU0[(:,k)]= complete_op(N0+1,N0+1,μk)*U0[(:,k)]
    end

    FU0 = LU0 - h*interval(0.5)*lambda_ext(U0^3,3K0,3N0) - β
    for k = 1:K0
        μk = h*interval(0.5)*(β1*interval(k)^4 + β2*interval(k)^2 - β0)
            FU0[(N0+1,k)] += μk*U0[(N0,k)]
    end

    A_LΛ = LinearOperator(Chebyshev(3*N)⊗CosFourier(K,interval(1)),S,interval.(zeros(dimension(S),dimension(Chebyshev(3N)⊗CosFourier(K,interval(1))))))

    L_FU0 = Sequence(Chebyshev(N)⊗CosFourier(K,interval(1)),interval.(zeros((N+1)*(K+1))))
    Y1 = interval(0)

    norm_N00_N1 = interval(0) 
    for k = 0:K

        μk = fct_μk(k)
        L_Λ, err_col = Linv_Lamb(μk,N00+N)
        # D_C10 = sum(abs,L_Λ[N+2:end,1:2*N+1],dims=1)+err_col[1:2*N+1]'    ; D_C10[1] = interval(2)*D_C10[1] 
        D_C10 = sum(abs,L_Λ[1:end,N+1:N00+N+1],dims=1)+err_col[N+1:N00+N+1]'    ; D_C10[1] = interval(2)*D_C10[1] 
        D_C1[(N:N00+N,k)] = D_C10
        
        g = interval.(zeros(maximum([N+1 3*N0+2])))
        g[1:3*N0+2],err_tail = Linv_V(FU0[(:,k)],μk,3*N0+1)
        L_FU0[(:,k)] = g[1:N+1]

        # display(norm(g[1:N+1],1))
        # in this bound we also incorporate the defect of the product L
        Y1 += interval(4)*(err_tail + norm(g[N+2:3*N0+2],1))

        # estimation of the sum of the columns from 1 to N after N00   
        norm_N00_N1 = maximum([norm_N00_N1 norm(L_Λ[1:N+2,N00],1)+err_col[N00]])

        A_LΛ[(:,:),(:,k)] = A[(:,:),(:,k)]*L_Λ[1:N+1,1:3N+1]
        D_A[(0:3N,k)] = sum(abs, A_LΛ[(:,:),(:,k)],dims=1) ; D_A[(0,k)] = interval(2)*D_A[(0,k)]
        DA0 = sum(abs,A[(:,:),(:,k)]*L_Λ[1:N+1,3N+2:N00+N+1],dims=1) ; 
        D_A[(3N+1:N00+N,k)] = DA0
        
    end

    for k = K+1:K0
            μk = fct_μk(k)       
            g = interval.(zeros(maximum([N+1 3*N0+2])))
            g[1:3*N0+2],err_tail = Linv_V(FU0[(:,k)],μk,3*N0+1)
            # in this bound we also incorporate the defect of the product L
            Y1 += interval(4)*(err_tail + norm(g[1:3*N0+2],1))
    end

    display("computations done")

    DQ = interval(3)*project(Multiplication(U0^2),S,Chebyshev(3*N)⊗CosFourier(K,interval(1)))
    Z0 = opnorm(I - A + h*interval(0.5)*A_LΛ*DQ,1)

    # display("value of Z0")
    # display(Z0)
    
    Z∞11, Z∞21 = component1_Z∞1_Z∞2(project(interval(3)*U0^2 - νinf,S))

    W0K = interval(0)
    for k = K+1:1000 
        μk = fct_μk(k)
        W0K = maximum([W0K function_C0(μk)/(interval(2)+μk)])
    end
    Z∞11 = maximum([Z∞11 norm(project(interval(3)*U0^2 - νinf,S),1)*W0K])

    # display("first Z1 bounds")
    # display(Z∞11)
    # display(Z∞21)
   
    Z∞12, Z∞22, ZN1, ZA01 = component2_Z∞1_Z∞2_ZN_ZA0(project(interval(3)*U0^2 - νinf,S),D_C1,D_A,N00)

    W0N00 = interval(0)
    for k=0:K 
        μk = fct_μk(k)
        W0N00 = maximum([W0N00 function_C1(μk)/(interval(N00+2) + sqrt(interval(N00+1)^2 + μk^2))])
    end
    Z∞22 = maximum([Z∞22 norm(project(interval(3)*U0^2 - νinf,S),1)*W0N00])
    ZA01 = maximum([ZA01 norm_A*norm(project(interval(3)*U0^2 - νinf,S),1)*norm_N00_N1])

    # display("second Z1 bounds")
    # display(Z∞12)
    # display(Z∞22)
    # display(ZN1)
    # display(ZA01)
    
    Z∞1 = h*interval(0.5)*(Z∞11 + Z∞12) 
    Z∞2 = h*interval(0.5)*(Z∞21 + Z∞22)
    Z∞ = maximum([Z∞1 Z∞2])
    ZN = h*interval(0.5)*(Z∞21 + ZN1)
    ZA0 = h*interval(0.5)*ZA01 

    Z1_1 = Z∞^2 + (interval(1) + Z∞)*ZN*ZA0 + ZA0*Z∞^2
    Z1_2 = Z0*(interval(1) + ZN + ZN*Z∞) + ZN*ZA0*(interval(1) + Z∞)

    # we also add the defect coming from the different size of truncation for the operators and the sequences. It is essentially the norm of A times the difference in ell 1 norm.
    Z1 = maximum([Z1_1 Z1_2]) + h*interval(0.5)*maximum([norm_A*(interval(1)+ZN*(interval(1)+Z∞)) (interval(1)+Z∞)*(interval(1)+ZA0)])*norm(interval(3)*U0^2-project(interval(3)*U0^2,S),1)

    display("value of Z1")
    display(Z1)
    # display(Z∞1)
    # display(Z∞2)
    # display(ZN)
    # display(ZA0)

####### Computation of the bound Y #####################


Y0 = norm(A*L_FU0,1) 

QU =  h*interval(0.5)*U0^3
for k=K0+1:3K0 
        μk = fct_μk(k)
        Y1 += interval(4)*function_C0(μk)/(interval(2)+μk)*norm(QU[(:,k)],1)
end

Y = Y0*(interval(1) + ZN*(interval(1) + Z∞)) + Y1*(interval(1) +Z∞)*(ZA0 + interval(1)) 

display("Y bound")
display(Y)

###############  computation of Z2   ##################################

r_star = interval(0.000001)
Z2 = h*(interval(6)*r_star + interval(3)*norm(U0^2,1))*maximum([norm_A*(interval(1)+ZN*(interval(1)+Z∞)) (interval(1)+Z∞)*(interval(1)+ZA0)])

display("Z2 bound")
display(Z2)

    if inf((interval(1)-Z1)^2) > sup(interval(2)*Z2*Y)
        r = (interval(1)-interval(sup(Z1)) - sqrt((interval(1)-interval(sup(Z1)))^2-interval(2)*Z2*Y))/Z2
        if sup(Z1 + Z2*abs(r)) < 1
            if sup(r)<0.000001
                display("proof successful for r =")
                display(r)
                return r
            else 
                display("r star too small")
                return Nan
            end
        else
            display("second condition not verified")
            return naN
        end
    else
        display("first condition not verified")
        return naN
    end  
    
    
end



