
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
    S = Chebyshev(N+1)⊗SinFourier(K,interval(1))⊗SinFourier(K,interval(1))
    F = Sequence(S,interval.(zeros((N+2)*(K^2))))
    for k1 ∈ 1:K
        for k2 ∈ 1:K
            for n ∈ 1:N-1
                F[(n,k1,k2)] = -U[(n-1,k1,k2)]+U[(n+1,k1,k2)]
            end

            F[(N,k1,k2)] = -U[(N-1,k1,k2)]
            F[(N+1,k1,k2)] = -U[(N,k1,k2)]
        end
    end
    return F
end



function lambda_ext(U,K,N,id)
    S = Chebyshev(N+1)⊗SinFourier(K,1)⊗SinFourier(K,1)
    F = Sequence(S,(zeros((N+2)*(K^2))))
    for k1 ∈ 1:K
        for k2 ∈ 1:K
            for n ∈ 1:N-1
                F[(n,k1,k2)] = -U[(n-1,k1,k2)]+U[(n+1,k1,k2)]
            end

            F[(N,k1,k2)] = -U[(N-1,k1,k2)]
            F[(N+1,k1,k2)] = -U[(N,k1,k2)]
        end
    end
    return F
end

function F(U0,N0,K0,β)
    LU0 = Sequence(Chebyshev(N0)⊗SinFourier(K0,1)⊗SinFourier(K0,(1)),(zeros((N0+1)*(K0)^2)))
    D1U0 = Sequence(Chebyshev(N0)⊗CosFourier(K0,(1))⊗SinFourier(K0,(1)),(zeros((N0+1)*(K0)*(K0+1))))
    D2U0 = Sequence(Chebyshev(N0)⊗SinFourier(K0,1)⊗CosFourier(K0,(1)),(zeros((N0+1)*(K0)*(K0+1))))
    for k1 = 1:K0
        for k2 = 1:K0
            LU0[(:,k1,k2)]= complete_op(N0+1,mid.(fct_μk(k1,k2)))*U0[(:,k1,k2)]
            D1U0[(:,k1,k2)] = k1/(k1^2 + k2^2)*U0[(:,k1,k2)]
            D2U0[(:,k1,k2)] = -k2/((k1)^2 + (k2)^2)*U0[(:,k1,k2)]
        end
    end

    return  project(LU0 + mid(h)*0.5*lambda_ext(Derivative(0,1,0)*(D2U0*U0) + Derivative(0,0,1)*(D1U0*U0),2K0,2N0,0) - mid.(β),Chebyshev(N0)⊗SinFourier(K0,1)⊗SinFourier(K0,(1)))
end


function DF(ν::Float64,U,S,S1,S2,ℒ,Λ,D1,D2,Der_1,Der_2)

    DF0 = (Der_1).*(mid.(conversion_op(project(Multiplication(U),S2, S1,Float64),S,S,[3 2])).*D2'+  mid.(conversion_op(project(Multiplication(mid.(conversion_sin_cos(D2.*U,S2,3))),S, S1,Float64),S,S,[0 2]))) 

    DF0 = DF0 + (Der_2).*(mid.(conversion_op(project(Multiplication(U),S1, S2,Float64),S,S,[2 3])).*D1' +  mid.(conversion_op(project(Multiplication(mid.(conversion_sin_cos(D1.*U,S1,2))),S, S2,Float64),S,S, [0 3])))

    return mid.(ℒ - mid(h)/2*Λ*DF0)
end


function newton_NS(U,N0,K0,km,precis,β)

    
S = Chebyshev(N0)⊗SinFourier(K0,1)⊗SinFourier(K0,1)
S1 = Chebyshev(N0)⊗CosFourier(K0,1)⊗SinFourier(K0,1)
S2 = Chebyshev(N0)⊗SinFourier(K0,1)⊗CosFourier(K0,1)

ℒ = LinearOperator(S,S,(zeros(dimension(S),dimension(S))))
Λ = LinearOperator(S,S,(zeros(dimension(S),dimension(S))))
D1 = Sequence(S,(zeros(dimension(S))))
D2  = Sequence(S,(zeros(dimension(S))))
Der_1  = Sequence(S,(zeros(dimension(S))))
Der_2  = Sequence(S,(zeros(dimension(S))))


for k1 = 1:K0
    for k2=1:K0
            μk = mid(h)*(0.5)*mid(ν)*((k1)^2+(k2)^2)
            Λ[(:,k1,k2),(:,k1,k2)]  = mid.(conversion(N0+1,N0+1))
            ℒ[(:,k1,k2),(:,k1,k2)] = mid.(complete_op(N0+1,N0+1,interval(μk)))

            D1[(:,k1,k2)] = k1/((k1)^2+(k2)^2)*((ones(N0+1)))
            D2[(:,k1,k2)] = -k2/((k1)^2+(k2)^2)*((ones(N0+1)))
            Der_1[(:,k1,k2)] = k1*((ones(N0+1)))
            Der_2[(:,k1,k2)] = k2*((ones(N0+1)))
    end
end
D1 = coefficients(D1) ; D2 = coefficients(D2) ; Der_1 = coefficients(Der_1) ; Der_2 = coefficients(Der_2)

    F1 = F(U,N0,K0,β)
    nf = norm(coefficients(F1),Inf)
    p=1
    display(nf)
    while (p<km)&&(nf>precis)
        DF1 = DF(mid.(ν),U,S,S1,S2,ℒ,Λ,D1,(D2),(Der_1),(Der_2))
        U = U- DF1\F1
        F1 = F(U,N0,K0,β)
        nf = norm(coefficients(F1),Inf)
        display(nf)
        p = p+1
    end
    return U
end


    

function conversion_sin_cos(U,S2,id)
    m = order(S2)

    if id==2
        V = Sequence(S2,interval.(zeros((m[1]+1)*(m[2]+1)*m[3])))
        for m2 = 1:m[2]
            V[(:,m2,:)] = U[(:,m2,:)]
        end
        return V
    else
        V = Sequence(S2,interval.(zeros((m[1]+1)*(m[3]+1)*m[2])))
        for m3 = 1:m[3]
            V[(:,:,m3)] = U[(:,:,m3)]
        end
        return V
    end
end


function conversion_cos_sin(U,S2,id)
    m = order(S2)

    if id==2
        V = Sequence(S2,interval.(zeros((m[1]+1)*(m[2])*m[3])))
        for m2 = 1:m[2]
            V[(:,m2,:)] = U[(:,m2,:)]
        end
        return V
    else
        V = Sequence(S2,interval.(zeros((m[1]+1)*(m[3])*m[2])))
        for m3 = 1:m[3]
            V[(:,:,m3)] = U[(:,:,m3)]
        end
        return V
    end
end

function conversion_op(M,S1,S2,id)

    oK1 = order(S1)[2] ; oK2 = order(S2)[2]
    if (id[1]==2)&&(id[2] ==3)
        V = LinearOperator(S1,S2,interval.(zeros(dimension(S2),dimension(S1))))
        for k2=1:oK1 
            for k3 = 1:oK2 
                V[(:,:,k3),(:,k2,:)] = M[(:,:,k3),(:,k2,:)]
            end
        end
        return V
    elseif (id[1]==3)&&(id[2] ==2)
        V = LinearOperator(S1,S2,interval.(zeros(dimension(S2),dimension(S1))))
        for k2=1:oK1 
            for k3 = 1:oK2 
                V[(:,k3,:),(:,:,k2)] = M[(:,k3,:),(:,:,k2)]
            end
        end
        return V
    elseif (id[1]==0)&&(id[2] ==2)
        V = LinearOperator(S1,S2,interval.(zeros(dimension(S2),dimension(S1))))
            for k3 = 1:oK2 
                V[(:,k3,:),(:,:,:)] = M[(:,k3,:),(:,:,:)]
            end
        return V
    elseif (id[1]==0)&&(id[2] ==3)
        V = LinearOperator(S1,S2,interval.(zeros(dimension(S2),dimension(S1))))
            for k3 = 1:oK2 
                V[(:,:,k3),(:,:,:)] = M[(:,:,k3),(:,:,:)]
            end
        return V
    else
        return Nan
    end
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



function fct_μk(k1,k2)
    return  h*interval(0.5)*ν*(interval(k1)^2+interval(k2)^2)
end

function non_L(n1,n2,k1,k2)
    return ExactReal(k1*n2-k2*n1)*( interval(1)/ExactReal(k1^2 + k2^2) - interval(1)/(ExactReal((n1- k1)^2 + (n2 - k2)^2)) )
end



function component_Zinf1(V)
    Z∞1 = interval(0) 
   Vc = Sequence(Chebyshev(N)⊗CosFourier(3K,interval(1))⊗CosFourier(3K,interval(1)),interval.(zeros((N+1)*(3*K+1)^2)))
   Vc[(:,1:K,1:K)] = V[(:,:,:)]
   d = interval.(2*ones(N+1)) ; d[1] = interval(1)
   
   for n1 = K+1:2*K
       for n2 =  1:2*K
           S = interval(0)
           for k1 =  n1-K:K+n1
               for k2 = maximum([1 n2-K]):K+n2
                   if (n1!=k1)||(n2!=k2)
                       μk = fct_μk(k1,k2) 
                       S += function_C0(μk)/(μk+interval(2))*norm(non_L(n1,n2,k1,k2)*d.*Vc[(:,abs(n1-k1),abs(n2-k2))],1)
                   end
               end
           end
           Z∞1 = maximum([Z∞1 S])
       end
   end
   return Z∞1
end



function component_Zinf2(V,D_C1,N00)
    Z∞2 = interval(0)
    
    d = interval.(2*ones(N+1)) ; d[1] = interval(1)
    Vc = Sequence(Chebyshev(N00+2N)⊗CosFourier(3K,interval(1))⊗CosFourier(3K,interval(1)),interval.(zeros(((N00+2N)+1)*(3*K+1)^2)))
    Vc[(0:N,1:K,1:K)] = V[(:,:,:)]

    for k1 = 1:K
        for k2 = 1:K
            S = interval(0)
            for j1 =  K+1:K+k1
                for j2 = 1:K+k2
                    μk = fct_μk(j1,j2) 
                    S += function_C1(μk)/(μk+interval(2))*norm(d.*( non_L(k1,k2,j1,j2)*Vc[(0:N,j1-k1,abs(k2-j2))] + non_L(k1,k2,j1,-j2)*Vc[(0:N,j1-k1,j2+k2)]) ,1)
                end
            end
            for j1 = 1:K
                for j2 =  K+1:K+k2
                    μk = fct_μk(j1,j2) 
                    S += function_C1(μk)/(μk+interval(2))*norm(d.*( non_L(k1,k2,j1,j2)*Vc[(0:N,abs(j1-k1),j2-k2)] + non_L(k1,k2,-j1,j2)*Vc[(0:N,j1+k1,abs(k2-j2))]) ,1)
                end
            end

            for n1 =  N+1:N00
                S13 = interval(0)
                
                for j1 = 1:k1-1
                    for j2 = 1:K
                        for n2 = N+1:N+n1
                            # we add up all the combinations (j1-k1)(j2-k2), (j1+k1)(j2-k2), (j1-k1)(j2+k2), (j1+k1)(j2+k2)
                            S13 += D_C1[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                        end
                    end
                end
                for j1 = k1+1:K
                    for j2 = 1:K
                        for n2 = N+1:N+n1
                            # we add up all the combinations (j1-k1)(j2-k2), (j1+k1)(j2-k2), (j1-k1)(j2+k2), (j1+k1)(j2+k2)
                            S13 += D_C1[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                        end
                    end
                end
                j1=k1
                for j2 = 1:k2-1
                    for n2 = N+1:N+n1
                        # we add up all the combinations (j1-k1)(j2-k2), (j1+k1)(j2-k2), (j1-k1)(j2+k2), (j1+k1)(j2+k2)
                        S13 += D_C1[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                    end
                end
                for j2 = k2+1:K
                    for n2 = N+1:N+n1
                        # we add up all the combinations (j1-k1)(j2-k2), (j1+k1)(j2-k2), (j1-k1)(j2+k2), (j1+k1)(j2+k2)
                        S13 += D_C1[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                    end
                end
                Z∞2 = maximum([Z∞2 S13+S])
            end 
        end
    end
    return Z∞2
end



function component_ZA0(V,D_A,N00)
    ZA0 = interval(0)
    Vc = Sequence(Chebyshev(N00+2N)⊗CosFourier(2K,interval(1))⊗CosFourier(2K,interval(1)),interval.(zeros(((N00+2N)+1)*(2*K+1)^2)))
    Vc[(0:N,1:K,1:K)] = V[(:,:,:)]

    for k1 = 1:K
        for k2 = 1:K
            for n1 =  N+1:N00
                S13 = interval(0)
                S15 = interval(0)
                for j1 = 1:k1-1
                    for j2 = 1:K
                        for n2 = 0:N+n1
                            # we add up all the combinations (j1-k1)(j2-k2), (j1+k1)(j2-k2), (j1-k1)(j2+k2), (j1+k1)(j2+k2)
                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                        end
                    end
                end
                for j1 = k1+1:K
                    for j2 = 1:K
                        for n2 = 0:N+n1
                            # we add up all the combinations (j1-k1)(j2-k2), (j1+k1)(j2-k2), (j1-k1)(j2+k2), (j1+k1)(j2+k2)
                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                        end
                    end
                end
                j1=k1
                for j2 = 1:k2-1
                    for n2 = 0:N+n1
                        S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                    end
                end
                for j2 = k2+1:K
                    for n2 = 0:N+n1
                        S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                    end
                end
                ZA0 = maximum([ZA0 S15])
            end 
        end
    end

    Vc = Sequence(Chebyshev(N00+2N)⊗CosFourier(3K,interval(1))⊗CosFourier(3K,interval(1)),interval.(zeros(((N00+2N)+1)*(3*K+1)^2)))
    Vc[(0:N,1:K,1:K)] = V[(:,:,:)]

    for k1 = K+1:2*K
        for k2 = 1:2K
            for n1 = 0:N00
                S15 = interval(0)
                for j1 = 1:K
                    for j2 = 1:K
                        S15 += D_A[(0,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(n1,abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(n1,abs(j1-k1),j2+k2)])
                        for n2 = 1:N
                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)] + non_L(j1,j2,k1,k2)*Vc[(n1+n2,abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(n1+n2,abs(j1-k1),j2+k2)])
                        end
                        for n2 = N+1:N+n1
                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)])
                        end
                    end 
                end
                ZA0 = maximum([ZA0 S15])
            end
        end
    end


    for k1 = 1:K
        for k2 = K+1:2K
            for n1 = 0:N00
                S15 = interval(0)
                for j1 = 1:K
                    for j2 = 1:K
                        S15 += D_A[(0,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(n1,abs(j1-k1),k2-j2)]  + non_L(j1,j2,-k1,k2)*Vc[(n1,j1+k1,k2-j2)])
                        for n2 = 1:N
                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),k2-j2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,k2-j2)] + non_L(j1,j2,k1,k2)*Vc[(n1+n2,abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,-k1,k2)*Vc[(n1+n2,j1+k1,k2-j2)])
                        end
                        for n2 = N+1:N+n1
                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,k2-j2)])
                        end
                    end 
                end
                ZA0 = maximum([ZA0 S15])
            end
        end
    end

    return ZA0 
end




function component_ZN(V,D_C1,N00)

    ZN = interval(0) 
    d = interval.(2*ones(N+1)) ; d[1] = interval(1)
    Vc = Sequence(Chebyshev(N00+2N)⊗CosFourier(3K,interval(1))⊗CosFourier(3K,interval(1)),interval.(zeros(((N00+2N)+1)*(3*K+1)^2)))
    Vc[(0:N,1:K,1:K)] = V[(:,:,:)]

    for k1 = 1:K 
        for k2 = 1:K 
            S = interval(0) 
            for j1 =  K+1:K+k1
                for j2 = 1:K+k2
                    μk = fct_μk(j1,j2) 
                    S += function_C1(μk)/(μk+interval(2))*norm(d.*( non_L(k1,k2,j1,j2)*Vc[(0:N,j1-k1,abs(k2-j2))] + non_L(k1,k2,j1,-j2)*Vc[(0:N,j1-k1,j2+k2)]) ,1)
                end
            end
            for j1 = 1:K
                for j2 =  K+1:K+k2
                    μk = fct_μk(j1,j2) 
                    S += function_C1(μk)/(μk+interval(2))*norm(d.*( non_L(k1,k2,j1,j2)*Vc[(0:N,abs(j1-k1),j2-k2)] + non_L(k1,k2,-j1,j2)*Vc[(0:N,j1+k1,abs(k2-j2))]) ,1)
                end
            end
            for n1 = 0:N
                S14 = interval(0)
                    for j1 = 1:K 
                        for j2 = 1:K
                            if (j1!=k1)||(j2!=k2)
                                for n2 = N+1:N+n1 
                                    S14 += D_C1[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(n2-n1,abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(n2-n1,abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(n2-n1,j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(n2-n1,j1+k1,j2+k2)] )
                                end
                            end
                        end
                    end 
                ZN = maximum([ZN S14+S])
            end
        end
        
    end

    return ZN 
end 




function component1_Z∞1(V)
    Z∞1 = interval(0) 
     Z∞2 = interval(0)
    Vc = Sequence(Chebyshev(N)⊗CosFourier(3K,interval(1))⊗CosFourier(3K,interval(1)),interval.(zeros((N+1)*(3*K+1)^2)))
    Vc[(:,1:K,1:K)] = V[(:,:,:)]
    d = interval.(2*ones(N+1)) ; d[1] = interval(1)
    
    for n1 = K+1:2*K
        for n2 =  1:2*K
            S = interval(0)
            for k1 =  n1-K:K+n1
                for k2 = maximum([1 n2-K]):K+n2
                    if (n1!=k1)||(n2!=k2)
                        μk = fct_μk(k1,k2) 
                        S += function_C0(μk)/(μk+interval(2))*norm(non_L(n1,n2,k1,k2)*d.*Vc[(:,abs(n1-k1),abs(n2-k2))],1)
                    end
                end
            end
            Z∞1 = maximum([Z∞1 S])
        end
    end

    # for n1 = 1:K
    #     for n2 =  K+1:4K
    #         S = interval(0)
    #         for k1 =  1:n1-1
    #             for k2 = K+1:K+n2
    #                 μk = fct_μk(k1,k2) 
    #                 S += function_C0(μk)/(μk+interval(2))*(norm(non_L(n1,n2,k1,k2)*d.*Vc[(:,abs(n1-k1),abs(n2-k2))] ,1) +  norm(non_L(n1,n2,-k1,k2)*d.*Vc[(:,n1+k1,abs(n2-k2))],1))
    #             end
    #         end
    #         for k1 =  n1+1:2K
    #             for k2 = K+1:2K
    #                 μk = fct_μk(k1,k2) 
    #                 S += function_C0(μk)/(μk+interval(2))*(norm(non_L(n1,n2,k1,k2)*d.*Vc[(:,abs(n1-k1),abs(n2-k2))] ,1) +  norm(non_L(n1,n2,-k1,k2)*d.*Vc[(:,n1+k1,abs(n2-k2))],1))
    #             end
    #         end
    #         k1=n1
    #         for k2 = K+1:n2-1
    #             μk = fct_μk(k1,k2) 
    #             S += function_C0(μk)/(μk+interval(2))*(norm(non_L(n1,n2,k1,k2)*d.*Vc[(:,abs(n1-k1),abs(n2-k2))] ,1) +  norm(non_L(n1,n2,-k1,k2)*d.*Vc[(:,n1+k1,abs(n2-k2))],1))
    #         end
    #         for k2 = n2+1:2K
    #             μk = fct_μk(k1,k2) 
    #             S += function_C0(μk)/(μk+interval(2))*(norm(non_L(n1,n2,k1,k2)*d.*Vc[(:,abs(n1-k1),abs(n2-k2))] ,1) +  norm(non_L(n1,n2,-k1,k2)*d.*Vc[(:,n1+k1,abs(n2-k2))],1))
    #         end
    #         for k1 = K+1:2K
    #             for k2 =  1:K
    #                     μk = fct_μk(k1,k2) 
    #                     S += function_C0(μk)/(μk+interval(2))*norm(non_L(n1,n2,k1,k2)*d.*Vc[(:,abs(n1-k1),abs(n2-k2))] ,1)
    #             end
    #         end
    #         Z∞1 = maximum([Z∞1 S])
    #     end
    # end
    
    for n1 = 1:K
        for n2 = 1:K
            S = interval(0)
            for k1 =  K+1:K+n1
                for k2 = 1:K+n2
                    μk = fct_μk(k1,k2) 
                    S += function_C0(μk)/(μk+interval(2))*norm(d.*( non_L(n1,n2,k1,k2)*Vc[(:,k1-n1,abs(n2-k2))] + non_L(n1,n2,k1,-k2)*Vc[(:,k1-n1,n2+k2)]) ,1)
                end
            end
            for k1 = 1:K
                for k2 =  K+1:K+n2
                    μk = fct_μk(k1,k2) 
                    S += function_C0(μk)/(μk+interval(2))*norm(d.*( non_L(n1,n2,k1,k2)*Vc[(:,abs(k1-n1),k2-n2)] + non_L(n1,n2,-k1,k2)*Vc[(:,k1+n1,abs(n2-k2))]) ,1)
                end
            end
            Z∞2 = maximum([Z∞2 S])
        end
    end

    return Z∞1, Z∞2
end





function component2_Z∞1_Z∞2_ZN_ZA0(V,D_C1,D_A,N00)
     Z∞2 = interval(0); ZN = interval(0) ;  ZA0 = interval(0)
    Vc = Sequence(Chebyshev(N)⊗CosFourier(3K,interval(1))⊗CosFourier(3K,interval(1)),interval.(zeros((N+1)*(3*K+1)^2)))
    Vc[(:,1:K,1:K)] = V[(:,:,:)]
    d = interval.(2*ones(N+1)) ; d[1] = interval(1)
    
    # Z∞1 can be improved by using C1 instead of C0. Aut if μk increases quick enough, that might be enough
    # we also use the symmetries of the sinus to reduce the number of loops
    # for n1 =  K+1:2*K
    #     for n2 =  1:2K
    #         S = interval(0)
    #         for k1 = 1:K
    #             for k2 = 1:K
    #                 μk = fct_μk(k1,k2) 
    #                 S += function_C0(μk)/(μk+interval(2))*norm(d.*( non_L(n1,n2,k1,k2)*Vc[(:,n1-k1,abs(n2-k2))] + non_L(n1,n2,k1,-k2)*Vc[(:,n1-k1,n2+k2)]) ,1)
    #             end
    #         end
    #         Z∞1 = maximum([Z∞1 S])
    #     end
    # end

    # for n1 =  1:K
    #     for n2 =  K+1:2K
    #         S = interval(0)
    #         for k1 = 1:K
    #             for k2 = 1:K
    #                 μk = fct_μk(k1,k2) 
    #                 S += function_C0(μk)/(μk+interval(2))*(norm(d.*(k1*(((n2-k2)/((n1-k1)^2+(n2-k2)^2) + n2/(n1^2+n2^2))*Vc[(:,abs(n1-k1),n2-k2)] )),1) + norm(d.*(k1*(((n2-k2)/((n1+k1)^2+(n2-k2)^2) + n2/(n1^2+n2^2))*Vc[(:,n1+k1,n2-k2)] )),1) )
    #             end
    #         end
    #         Z∞1 = maximum([Z∞1 S])
    #     end
    # end

    Vc = Sequence(Chebyshev(N00+2N)⊗CosFourier(2K,interval(1))⊗CosFourier(2K,interval(1)),interval.(zeros(((N00+2N)+1)*(2*K+1)^2)))
    Vc[(0:N,1:K,1:K)] = V[(:,:,:)]

    for k1 = 1:K
        for k2 = 1:K
            for n1 =  N+1:N00
                S13 = interval(0)
                S15 = interval(0)
                for j1 = 1:k1-1
                    for j2 = 1:K
                        for n2 = 0:N+n1
                            # we add up all the combinations (j1-k1)(j2-k2), (j1+k1)(j2-k2), (j1-k1)(j2+k2), (j1+k1)(j2+k2)
                            S13 += D_C1[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )

                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                        end
                    end
                end
                for j1 = k1+1:K
                    for j2 = 1:K
                        for n2 = 0:N+n1
                            # we add up all the combinations (j1-k1)(j2-k2), (j1+k1)(j2-k2), (j1-k1)(j2+k2), (j1+k1)(j2+k2)
                            S13 += D_C1[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )

                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                        end
                    end
                end
                j1=k1
                for j2 = 1:k2-1
                    for n2 = 0:N+n1
                        # we add up all the combinations (j1-k1)(j2-k2), (j1+k1)(j2-k2), (j1-k1)(j2+k2), (j1+k1)(j2+k2)
                        S13 += D_C1[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )

                        S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                    end
                end
                for j2 = k2+1:K
                    for n2 = 0:N+n1
                        # we add up all the combinations (j1-k1)(j2-k2), (j1+k1)(j2-k2), (j1-k1)(j2+k2), (j1+k1)(j2+k2)
                        S13 += D_C1[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )

                        S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(abs(n1-n2),j1+k1,j2+k2)] )
                    end
                end
                Z∞2 = maximum([Z∞2 S13])
                ZA0 = maximum([ZA0 S15])
            end 
        
            for n1 = 0:N
                S14 = interval(0)
                for j1 = 1:K 
                    for j2 = 1:K
                        if (j1!=k1)||(j2!=k2)
                        for n2 = N+1:N+n1 
                            S14 += D_C1[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(n2-n1,abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(n2-n1,abs(j1-k1),j2+k2)]  + non_L(j1,j2,-k1,k2)*Vc[(n2-n1,j1+k1,abs(j2-k2))]  + non_L(j1,j2,-k1,-k2)*Vc[(n2-n1,j1+k1,j2+k2)] )
                        end
                        end
                    end
                end
                ZN = maximum([ZN S14])
            end
        end
    end

    Vc = Sequence(Chebyshev(N00+2N)⊗CosFourier(3K,interval(1))⊗CosFourier(3K,interval(1)),interval.(zeros(((N00+2N)+1)*(3*K+1)^2)))
    Vc[(0:N,1:K,1:K)] = V[(:,:,:)]

    for k1 = K+1:2*K
        for k2 = 1:2K
            for n1 = 0:N00
                S15 = interval(0)
                for j1 = 1:K
                    for j2 = 1:K
                        S15 += D_A[(0,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(n1,abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(n1,abs(j1-k1),j2+k2)])
                        for n2 = 1:N
                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)] + non_L(j1,j2,k1,k2)*Vc[(n1+n2,abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(n1+n2,abs(j1-k1),j2+k2)])
                        end
                        for n2 = N+1:N+n1
                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,k1,-k2)*Vc[(abs(n1-n2),abs(j1-k1),j2+k2)])
                        end
                    end 
                end
                ZA0 = maximum([ZA0 S15])
            end
        end
    end


    for k1 = 1:K
        for k2 = K+1:2K
            for n1 = 0:N00
                S15 = interval(0)
                for j1 = 1:K
                    for j2 = 1:K
                        S15 += D_A[(0,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(n1,abs(j1-k1),k2-j2)]  + non_L(j1,j2,-k1,k2)*Vc[(n1,j1+k1,k2-j2)])
                        for n2 = 1:N
                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),k2-j2)]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,k2-j2)] + non_L(j1,j2,k1,k2)*Vc[(n1+n2,abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,-k1,k2)*Vc[(n1+n2,j1+k1,k2-j2)])
                        end
                        for n2 = N+1:N+n1
                            S15 += D_A[(n2,j1,j2)]*abs( non_L(j1,j2,k1,k2)*Vc[(abs(n1-n2),abs(j1-k1),abs(j2-k2))]  + non_L(j1,j2,-k1,k2)*Vc[(abs(n1-n2),j1+k1,k2-j2)])
                        end
                    end 
                end
                ZA0 = maximum([ZA0 S15])
            end
        end
    end

    return Z∞2, ZN, ZA0
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

    # first column of Linv
    L_col1 = interval.(zeros(N+1))
    L_col1[1] = interval(1)/(interval(1)+interval(2)*norm_col_1*μ)
    L_col1[2:end] = -μ*Tinv[:,1]/(interval(1)+interval(2)*norm_col_1*μ)
    # error from the approximate inverse (that could be better but ok in practice)
    L_col1 = interval.(inf.(L_col1-abs.(L_col1)*nT/(interval(1)-nT)), sup.(L_col1+abs.(L_col1)*nT/(interval(1)-nT)))

    # propagation of error from the term θ0^T Tinv Λ 
    T_err = [exp(function_g(interval(i-2)/μ,μ) - function_g(interval(N-2)/μ,μ))/(interval(N-1) + sqrt(interval(N-1)^2 + μ^2))+exp(function_g(interval(i-2)/μ,μ) - function_g(interval(N-1)/μ,μ))/(interval(N) + sqrt(interval(N)^2 + μ^2)) for i∈1:N]' 

    # error propagation in each column coming from the first column of Tinv
    err_col = function_C0(μ)/(interval(2)+μ)*μ/(interval(1)+interval(2)*norm_col_1*μ)*first_column_residual(N+2,μ)*ones(N+1) ; err_col[1] =  interval(1)/(interval(1)+interval(2)*norm_col_1*μ)*first_column_residual(N+2,μ)
    err_col[2:end] = err_col[2:end] + T_err'

    T_err = interval.(-sup.(T_err),sup.(T_err))

     return [-interval(2)*norm_col_1/(interval(1)+interval(2)*norm_col_1*μ) -theta*Tlamb/(interval(1)+interval(2)*norm_col_1*μ)+T_err;
             -Tinv[:,1]/(interval(1)+interval(2)*norm_col_1*μ)  (I + μ/(interval(1)+interval(2)*norm_col_1*μ)*Tinv[:,1]*theta)*Tlamb+μ/(interval(1)+interval(2)*norm_col_1*μ)*Tinv[:,1]*T_err], err_col, L_col1, nT
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




function norm_A_L_col_k(k1,k2,U0,D1U0,D2U0,Lcol0,N00)
    d = interval(2*ones(N00)) ; d[1] = interval(1)
    D1Lcol = Sequence(Chebyshev(N00)⊗CosFourier(K0,interval(1))⊗SinFourier(K0,interval(1)),interval.(zeros((N00+1)*(K0)*(K0+1))))
    D2Lcol= Sequence(Chebyshev(N00)⊗SinFourier(K0,interval(1))⊗CosFourier(K0,interval(1)),interval.(zeros((N00+1)*(K0)*(K0+1))))

    D2Lcol[(:,k1,k2)] = -interval(k2)/(interval(k1)^2 + interval(k2)^2)*Lcol0[(:,k1,k2)]
    D1Lcol[(:,k1,k2)] = interval(k1)/(interval(k1)^2 + interval(k2)^2)*Lcol0[(:,k1,k2)]
    
    col = Lcol0 - h*interval(0.5)*lambda_ext(Derivative(0,1,0)*(D2U0*Lcol0 + U0*D2Lcol) + Derivative(0,0,1)*(D1U0*Lcol0 + U0*D1Lcol),2K0,2N0)
    return norm(d.*col,1)
end 



function norm_A_L_col0(U0,D1U0,D2U0,N00,Z10)
    norm_col = interval(0)
    for k1 = 1:2K0 
        for k2 = 1:2K0 
            μ = fct_μk(k1,k2)
            T = Matrix(tridiag_int(N00,N00,μ))
            Tinv = interval.(inv(mid.(T)))
            nT = opnorm(I - Tinv*T,1)

            N = Int(floor(mid.(μ)/5)) + 100
            norm_col_1 = norm_first_column(N+1,μ)
            

            # first column of Linv
            L_col = interval.(zeros(N00+1))
            L_col[1] = interval(1)/(interval(1)+interval(2)*norm_col_1*μ)
            L_col[2:end] = -μ*Tinv[:,1]/(interval(1)+interval(2)*norm_col_1*μ)

            err_col =  interval(1)/(interval(1)+interval(2)*norm_col_1*μ)*first_column_residual(N00+2,μ)

            L_col0 = Sequence(Chebyshev(N00)⊗SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1)),interval.(zeros((N00+1)*(K0)^2)))
            L_col0[(:,k1,k2)] = L_col
            norm_col = maximum([norm_A_L_col_k(k1,k2,U0,D1U0,D2U0,L_col0,N00)/(interval(1)-nT)+Z10*(err_col) norm_col])
            display(k1)
            display(k2)
        end 
    end 

    return norm_col

end 






function proof_integration_NS_sinus_Z1_new(U0,h,ν,β,N0,K0,N,K,rold)

    
S = Chebyshev(N)⊗SinFourier(K,interval(1))⊗SinFourier(K,interval(1))
S1 = Chebyshev(N)⊗CosFourier(K,interval(1))⊗SinFourier(K,interval(1))
S2 = Chebyshev(N)⊗SinFourier(K,interval(1))⊗CosFourier(K,interval(1))

ℒ = LinearOperator(S,S,interval.(zeros(dimension(S),dimension(S))))
Λ = LinearOperator(S,S,interval.(zeros(dimension(S),dimension(S))))
D1 = Sequence(S,interval.(zeros(dimension(S))))
D2  = Sequence(S,interval.(zeros(dimension(S))))
Der_1 = Sequence(S,interval.(zeros(dimension(S))))
Der_2  = Sequence(S,interval.(zeros(dimension(S))))
Ω = interval(2)*Sequence(S,interval.(ones(dimension(S))))

for k1 = 1:K
    for k2=1:K
            μk = h*interval(0.5)*ν*(interval(k1)^2+interval(k2)^2)
            Λ[(:,k1,k2),(:,k1,k2)]  = conversion(N+1,N+1)
            ℒ[(:,k1,k2),(:,k1,k2)] = complete_op(N+1,N+1,μk)
            D1[(:,k1,k2)] = interval(k1)/(interval(k1)^2+interval(k2)^2)*interval.((ones(N+1)))
            D2[(:,k1,k2)] = -interval(k2)/(interval(k1)^2+interval(k2)^2)*interval.((ones(N+1)))
            Der_1[(:,k1,k2)] = interval(k1)*interval.((ones(N+1)))
            Der_2[(:,k1,k2)] = interval(k2)*interval.((ones(N+1)))
            Ω[(0,k1,k2)] = interval(1)
    end
end
D1 = coefficients(D1) ; D2 = coefficients(D2) ; Der_1 = coefficients(Der_1) ; Der_2 = coefficients(Der_2); Ω = coefficients(Ω)

Der_10  = Sequence(Chebyshev(N0)⊗SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1)),interval.(zeros(dimension(Chebyshev(N0)⊗SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1))))))
Der_20  = Sequence(Chebyshev(N0)⊗SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1)),interval.(zeros(dimension(Chebyshev(N0)⊗SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1))))))
for k1 = 1:K0
    for k2=1:K0
            Der_10[(:,k1,k2)] = interval(k1)*interval.((ones(N0+1)))
            Der_20[(:,k1,k2)] = interval(k2)*interval.((ones(N0+1)))
    end
end
Der_10 = coefficients(Der_10) ; Der_20 = coefficients(Der_20)
Δ2 = sqrt.(Der_10.^2 + Der_20.^2) ; norm_U0_Δ2 = norm(Δ2.*U0,1) ; norm_U0_inv_Δ2 = norm_U0_Δ2 = norm(U0./Δ2,1)
Der_10 = Nothing ; Der_20 = Nothing


Sext = Chebyshev(2N)⊗SinFourier(K,interval(1))⊗SinFourier(K,interval(1))
Der_1  = Sequence(Sext,interval.(zeros(dimension(Sext))))
Der_2  = Sequence(Sext,interval.(zeros(dimension(Sext))))

for k1 = 1:K
    for k2=1:K
            Der_1[(:,k1,k2)] = interval(k1)*interval.((ones(2N+1)))
            Der_2[(:,k1,k2)] = interval(k2)*interval.((ones(2N+1)))
    end
end
Der_1 = coefficients(Der_1) ; Der_2 = coefficients(Der_2);

    # Computation of the bound Z1
    N00 = maximum([2*N+2 50])

    LU0 = Sequence(Chebyshev(N0)⊗SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1)),interval.(zeros((N0+1)*(K0)^2)))
    D1U0 = Sequence(Chebyshev(N0)⊗CosFourier(K0,interval(1))⊗SinFourier(K0,interval(1)),interval.(zeros((N0+1)*(K0)*(K0+1))))
    D2U0 = Sequence(Chebyshev(N0)⊗SinFourier(K0,interval(1))⊗CosFourier(K0,interval(1)),interval.(zeros((N0+1)*(K0)*(K0+1))))
    for k1 = 1:K0
        for k2 = 1:K0
            LU0[(:,k1,k2)]= complete_op(N0+1,N0+1,fct_μk(k1,k2))*U0[(:,k1,k2)]
            D1U0[(:,k1,k2)] = interval(k1)/(interval(k1)^2 + interval(k2)^2)*U0[(:,k1,k2)]
            D2U0[(:,k1,k2)] = -interval(k2)/(interval(k1)^2 + interval(k2)^2)*U0[(:,k1,k2)]
        end
    end

    FU0 = LU0 + h*interval(0.5)*lambda_ext(Derivative(0,1,0)*(D2U0*U0) + Derivative(0,0,1)*(D1U0*U0),2K0,2N0) - β

    
    for k1 = 1:K0
        for k2 = 1:K0
            μk = fct_μk(k1,k2)
            FU0[(N0+1,k1,k2)] += μk*U0[(N0,k1,k2)]
        end
    end

    L_FU0 = Sequence(Chebyshev(N)⊗SinFourier(K0,interval(1))⊗SinFourier(K0,interval(1)),interval.(zeros((N+1)*(K0)^2)))
    Y1 = interval(0)
    for k1 = 1:K0
        for k2 = 1:K0
        μk = fct_μk(k1,k2)
        
        g = interval.(zeros(maximum([N+1 2*N0+2])))
        g[1:2*N0+2],err_tail = Linv_V(FU0[(:,k1,k2)],μk,2*N0+1)
        L_FU0[(:,k1,k2)] = g[1:N+1]

        # in this bound we also incorporate the defect of the product L
        Y1 += interval(8)*(err_tail + norm(g[N+2:2*N0+2],1))
        end
    end

    display("computations done")

    
    # Computation of Z1
    # the maximum of h|k|/(2+νh|k|^2) at reached at |k| = sqrt(2/(νh)) \leq 20 in our case
    W0k = interval(0)
    for k1=1:20
        for k2=1:20
            μk = fct_μk(k1,k2)
            W0k = maximum([W0k function_C0(μk)*sqrt(interval(k1)^2+interval(k2)^2)/((interval(2)+μk))])
        end
    end

    # display("components Z∞")
    Z10 = h*interval(0.5)*(norm_U0_inv_Δ2*W0k+norm_U0_Δ2*W0k*interval(0.25))
    Z1 = interval(sup(Z10^2))
    display("value of Z1")
    display(Z1)
    # display(Z∞1)
    # display(Z∞2)
    # display(ZN)
    # display(ZA0)

    # norm_A_L_col = norm_A_L_col0(U0,D1U0,D2U0,N00,Z10)
    # display("norm A col L 0")
    # display(norm_A_L_col)

####### Computation of the bound Y #####################


Y0 = (interval(1)+Z10)*norm(L_FU0,1) 

QU =  h*interval(0.5)*(Derivative(0,1,0)*(D2U0*U0) + Derivative(0,0,1)*(D1U0*U0))
for k1=K0+1:2K0 
    for k2 = 1:K0
        μk = fct_μk(k1,k2)
        Y1 += interval(8)*function_C0(μk)/(interval(2)+μk)*norm(QU[(:,k1,k2)],1)
    end
end

for k1=K0+1:2K0 
    for k2 = K0+1:2K0
        μk = fct_μk(k1,k2)
        Y1 += interval(8)*function_C0(μk)/(interval(2)+μk)*norm(QU[(:,k1,k2)],1)
    end
end

# the term with rold is the propagation of error from the previous step as described in appendix.
# first we compute the max of \|A(ℒkinv)_col(1)\|

max_A_L_col1 = interval(1)+Z10
Y = interval(sup(Y0 + Y1 + rold*max_A_L_col1))

display("Y bound")
display(Y)
# display(Y0)
# display(Y1)


#################################################


# ##############################################

# Z21 = WN
# Z22 = h*interval(0.5)*nk2
# Z23 = h*interval(0.5)*Z∞2
# Z24 = h*interval(0.5)*ZN

# Z2 = interval(sup(Z21*(interval(1) + nk1) + Z22 + Z23 + Z24))


Z2 = interval(sup(h*interval(0.5)*(norm_U0_inv_Δ2*W0k+interval(0.25)*norm_U0_Δ2*W0k)*(interval(1)+Z10)))
display("Z2 bound")
display(Z2)


    if inf((interval(1)-Z1)^2) > sup(interval(2)*Z2*Y)
        r = (interval(1)-interval(sup(Z1)) - sqrt((interval(1)-interval(sup(Z1)))^2-interval(2)*Z2*Y))/Z2
        if sup(Z1 + Z2*abs(r)) < 1
            display("proof successful for r =")
            display(r)
            return interval(sup(r))
        else
            display("second condition not verified")
            return naN
        end
    else
        display("first condition not verified")
        return naN
    end  
    
    
end



