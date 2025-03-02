using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, JLD2

IntervalArithmetic.matmul_mode() = IntervalArithmetic.MatMulMode{:fast}()


function complete_op(N,M,μ::Interval{Float64})
     T = interval.(zeros(N,M))
     T[1,1] = interval(1)
     T[1,2:M] = interval(2)*interval.(((-1).^(1:M-1))')
     T[2,1] = μ
     T[2:end,2:end] = Matrix(tridiag_int(N-1,M-1,μ))
     return T
 end


 function complete_op(N,M,μ::Float64)
    T = zeros(N,M)
    T[1,1] = 1
    T[1,2:M] = 2*((-1).^(1:M-1))'
    T[2,1] = μ
    T[2:end,2:end] = Matrix(tridiag_int(N-1,M-1,μ))
    return T
end

function complete_op(N,M,μ::Int)
    T = zeros(N,M)
    T[1,1] = 1
    T[1,2:M] = 2*((-1).^(1:M-1))'
    T[2,1] = μ
    T[2:end,2:end] = Matrix(tridiag_int(N-1,M-1,μ))
    return T
end

function tridiag_int(N,M,μ::Interval{Float64})
    if N==M
         b = μ*interval.(ones(N-1))
         a = interval(2)*interval.(collect(1:N))   
         return Tridiagonal(b,a,-b)
    elseif M>N
         T = interval.(zeros(N,M))
         b = μ*interval.(ones(N-1))
         a = interval(2)*interval.(collect(1:N))   
         T[1:N,1:N] = Tridiagonal(b,a,-b)
         T[N,N+1] = -μ
         return T
    else
         display("need to choose the number of columns bigger")
    end
 end


 function tridiag_int(N,M,μ::Float64)
    if N==M
         b = μ*ones(N-1)
         a = 2*interval.(collect(1:N))   
         return Tridiagonal(b,a,-b)
    elseif M>N
         T = zeros(N,M)
         b = μ*ones(N-1)
         a = 2*collect(1:N)  
         T[1:N,1:N] = Tridiagonal(b,a,-b)
         T[N,N+1] = -μ
         return T
    else
         display("need to choose the number of columns bigger")
    end
 end


 function tridiag_int(N,M,μ::Int)
    if N==M
         b = μ*ones(N-1)
         a = 2.0*collect(1:N)   
         return Tridiagonal(b,a,-b)
    elseif M>N
         T = zeros(N,M)
         b = μ*ones(N-1)
         a = 2*collect(1:N)  
         T[1:N,1:N] = Tridiagonal(b,a,-b)
         T[N,N+1] = -μ
         return T
    else
         display("need to choose the number of columns bigger")
    end
 end

 function conversion(N,M)
     T = interval.(zeros(N,M))
     for i=1:N
         for j=1:N
             if j==i+1
                 T[i,j] = interval(1)
             elseif i==j+1
                 T[i,j] =-interval(1) 
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
 



###################### Computation of the desired constants #################################################


 μ_start = 0 ;  μ_end = 1 ; dμ = 0.01 ; N = 50
 J = μ_start:dμ:μ_end  ; C0_vec = interval.(zeros(size(J))) ; C1_vec = interval.(zeros(size(J)))

# μ_start = 1 ;  μ_end = 20 ; dμ = 0.1 ; N = 100
# J = μ_start:dμ:μ_end  ; C0_vec = interval.(zeros(size(J))) ; C1_vec = interval.(zeros(size(J)))

# μ_start = 20 ;  μ_end = 100 ; dμ = 0.5 ; N = 300
# J = μ_start:dμ:μ_end  ; C0_vec = interval.(zeros(size(J))) ; C1_vec = interval.(zeros(size(J)))

# μ_start = 100 ;  μ_end = 200 ; dμ = 1 ; N = 600
# J = μ_start:dμ:μ_end  ; C0_vec = interval.(zeros(size(J))) ; C1_vec = interval.(zeros(size(J)))

# μ_start = 200 ;  μ_end = 300 ; dμ = 5 ; N = 900
# J = μ_start:dμ:μ_end  ; C0_vec = interval.(zeros(size(J))) ; C1_vec = interval.(zeros(size(J)))

# μ_start = 305 ;  μ_end = 400 ; dμ = 1 ; N = 1200
# J = μ_start:dμ:μ_end  ; C0_vec = interval.(zeros(size(J))) ; C1_vec = interval.(zeros(size(J)))

#  μ_start = 410 ;  μ_end = 700 ; dμ = 10 ; N = 2100
#  J = μ_start:dμ:μ_end  ; C0_vec = interval.(zeros(size(J))) ; C1_vec = interval.(zeros(size(J)))

#  μ_start = 710 ;  μ_end = 1000 ; dμ = 15 ; N = 3000
#  J = μ_start:dμ:μ_end  ; C0_vec = interval.(zeros(size(J))) ; C1_vec = interval.(zeros(size(J)))

#  μ_start = 1000 ;  μ_end = 2000 ; dμ = 50 ; N = 6000
#  J = μ_start:dμ:μ_end  ; C0_vec = interval.(zeros(size(J))) ; C1_vec = interval.(zeros(size(J)))

#  μ_start = 2000 ;  μ_end = 3000 ; dμ = 100 ; N = 9000
#  J = μ_start:dμ:μ_end  ; C0_vec = interval.(zeros(size(J))) ; C1_vec = interval.(zeros(size(J)))

# I = 1:dx:x_lim ; C0 = interval.(zeros(size(I))) ; p = 1
# D = Diagonal(interval(2)*interval.(0:N))
# θ = -interval(2.0)*interval.(((-1).^(1:N+1))')
# θ[1,1] = interval(1)

# 𝒟 = Matrix(D) ; 𝒟[1,:] = θ[1,:]

W0 = interval.(2*ones(N+1)); W0[1] = interval(1) ; W0inv = interval(0.5)*interval.(ones(N+1)); W0inv[1] = interval(1) 
D0 = interval(0:N) 


p=1
for μ ∈ J
    global p
    μ = interval(μ)
    D = (D0 .+ interval(1) + sqrt( μ^2 .+ D.^2)).*W0inv
    L = complete_op(N+1,N+1,μ)
    Linv = interval.(inv(mid.(L)))


    nL = opnorm(I - W0.*LinearOperator(Linv)*LinearOperator(L).*W0inv',1)
     v = W0.*Linv[:,1] ; nv = norm(v,1)
     ρ1 = interval(1)/(interval(N+1))*(nv + interval(1))
     ρ2 = norm(W0.*Linv[:,end] + interval(1)/interval(N+2)*v,1) + interval(1)/interval(N+2)
     ρ3 = interval(2)*nv/(interval((N+2)*(N+3))) + interval(1)/interval(N+1) + interval(1)/interval(N+3)

     ρ = interval(0.5)*μ*maximum([ρ1 ρ2 ρ3])

    nvD = norm(D.*Linv[:,1],1)

    δ1 = nvD/(interval(N+1)*(μ+interval(2)*interval(N+1))) + interval(1)/(interval(2)*interval(N+1))
    δ2 = norm(D.*Linv[:,end],1)/(interval(2)*interval(N+2) + μ) + interval(1)/(interval(2)*interval(N+3))
    δ3 = interval(1)/(interval(N+1))

    δ = μ*maximum([δ1 δ2 δ3])

    if (sup(δ)>1)||(sup(ρ)>1)
        display("The value of N is not big enough")
        brk = breaknow
    end

    LV =  LinearOperator(Linv)*LinearOperator(conversion(N+1,N+1))
    C00 = (μ+interval(2))*opnorm(W0.*LV.*W0inv',1)
    C11 = opnorm(W0.*LV.*D',1)

    err1 =  (μ+ interval(2))*( nv/interval(N+1) + interval(1)/(interval(2)*interval(N+1)) + norm(W0.*LV[:,N+1],1) )
    # We use that D = (L - e0θ^T + \mu \Lambda)
    err2 = maximum([err1+interval(1)/(interval(N+1)) μ/(interval(2)*interval(N+1))+(μ+interval(2))/(interval(2*(N+2)))+μ/(interval(2*(N+3)))])

    
    C0 = interval(1)/(interval(1)-nL)*maximum([C00 interval(1)/(interval(1)-ρ)*err1 interval(1)/(interval(1)-ρ)*ρ3 interval(1)/(interval(1)-ρ)*ρ1])
    C1 = interval(1)/(interval(1)-nL)*maximum([C11 interval(1)/(interval(1)-δ)*err2 interval(2)])

    C0_vec[p] = C0 ; C1_vec[p] = C1

    display("value of C0")
    display(mid(C0))

    display("value of C1")
    display(mid(C11))
    display(mid(C1))

    display("value of p")
    display(p)
    p = p+1
end
    
# save the obtained constants
# jldsave("C0_10.jld2" ; C0_vec)
# jldsave("C1_10.jld2" ; C1_vec)


# μ = 0:0.01:6000

# C0 = zeros(length(μ))
# C1 = zeros(length(μ))

# p = 1
# for k ∈ μ
#     global μ,p,C0,C1
#     # C0[p] = mid.(function_C0(k)) 
#     C1[p] = mid.(function_C1(k)) 
#     p = p+1
# end

# using MATLAB

#     # mat"
#     # plot($μ,$C0)"

#  mat"
#      plot($μ,$C1)"