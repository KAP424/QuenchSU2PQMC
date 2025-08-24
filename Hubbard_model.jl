
struct _Hubbard_Para
    Lattice::String
    t::Float64
    U::Vector{Float64}
    Rate::Float64
    site::Vector{Int64}
    Ns::Int64
    Nt::Int64
    K::Array{Float64,2}
    BatchSize::Int64
    WrapTime::Int64
    Δt::Float64
    α::Vector{Float64}
    γ::Vector{Float64}
    η::Vector{Float64}
    Pt::Array{Float64,2}
    HalfeK::Array{Float64,2}
    eK::Array{Float64,2}
    HalfeKinv::Array{Float64,2}
    eKinv::Array{Float64,2}
end

function Hubbard_Para(t,U0,Uend,Rate,Lattice::String,site,Δt,BatchSize)
    Nt::Int64=2*Int(round((Uend-U0)/Δt*Rate))
    U=U0.+vcat( collect(1:cld(Nt,2)),collect(cld(Nt,2):-1:1) ).*Rate*Δt
    
    α=zeros(Nt)
    for i in 1:Nt
        α[i]=sqrt(Δt*U[i]/2)
    end

    WrapTime::Int64=div(BatchSize,2)
    
    γ::Vector{Float64}=[1+sqrt(6)/3,1+sqrt(6)/3,1-sqrt(6)/3,1-sqrt(6)/3]
    η::Vector{Float64}=[sqrt(2*(3-sqrt(6))),-sqrt(2*(3-sqrt(6))),sqrt(2*(3+sqrt(6))),-sqrt(2*(3+sqrt(6)))]
    
    K::Array{Float64,2}=t.*K_Matrix(Lattice,site)
    Ns::Int64=size(K)[1]

    μ=0.0
    if Lattice=="HoneyComb"
        K+=μ*diagm(repeat([-1, 1], div(Ns, 2)))
    elseif Lattice=="SQUARE"
        for i in 1:Ns
            x,y=i_xy(Lattice,site,i)
            K[i,i]+=μ*(-1)^(x+y)
        end
    end
    # K[K .!= 0] .+=( rand(size(K)...) * 0.1)[K.!= 0]
    # K=(K+K')./2

    E,V=eigen(K)
    
    HalfeK::Array{Float64,2}=V*diagm(exp.(-Δt.*E./2))*V'
    eK::Array{Float64,2}=V*diagm(exp.(-Δt.*E))*V'
    HalfeKinv::Array{Float64,2}=V*diagm(exp.(Δt.*E./2))*V'
    eKinv::Array{Float64,2}=V*diagm(exp.(Δt.*E))*V'

    # 从无序费米液体基态-->反铁磁序直积态
    # U:0->Θ*Rate
    Pt=V[:,1:div(Ns,2)]
    

    return _Hubbard_Para(Lattice,t,U,Rate,site,Ns,Nt,K,BatchSize,WrapTime,Δt,α,γ,η,Pt,HalfeK,eK,HalfeKinv,eKinv)

end

