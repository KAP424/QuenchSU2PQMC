include("E:/桌面/JuliaDQMC/JuliaDQMC/tool.jl")
push!(LOAD_PATH, "E:/桌面/JuliaDQMC/JuliaDQMC/quench/")
using LinearAlgebra
using Plots
using DelimitedFiles 
using LaTeXStrings
using Statistics
using KAPDQMC  
using Random

# t=1.0;   Lattice="HoneyComb"    
# U=3.8;     Δt=0.1;     Θ=0.1;
# BatchSize=100;  L=6;
# site=[L,L];

# model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"H0")

# # Half
# indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
# # println(indexA)

# # HalfHalf
# indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))

# rng=MersenneTwister(Threads.threadid())
# s1=Initial_s(model,rng)
# s2=Initial_s(model,rng)
# lambda=0.0;            N=2;
# ss=ctrl_SCEEicr("E:/桌面/",model,indexA,indexB,1,lambda,N,[s1[:,:],s2[:,:]],false)
# println(typeof(ss))
# ctrl_SCEEicr("E:/桌面/",model,indexA,indexB,1,lambda,N,ss,true)


t=1.0;      V_quench=1.0;    Lattice="HoneyComb"    
Δt=0.1;     Θ=1;
BatchSize=10;  L=3;
site=[L,L];
model=Quench_Hubbard(t,V_quench,Lattice,site,Δt,Θ,BatchSize)

# Half
indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
# println(indexA)

# HalfHalf
indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))

rng=MersenneTwister(Threads.threadid())
s1=Initial_s(model,rng)
s2=Initial_s(model,rng)

lambda=0.0
N=1
ss=ctrl_EEicr("E:/桌面/JuliaDQMC/JuliaDQMC/quench/",model,indexA,1,lambda,N,[s1[:,:],s2[:,:]],true)



println(model.U)
println(model.α)

Nt=40
NO=div(Nt,2)
for lt in 1:Nt
    print(Int(NO+0.5-abs(lt-NO-0.5)))
end