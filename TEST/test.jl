push!(LOAD_PATH,"E:/桌面/JuliaDQMC/code/QuencSU2PQMC/")
using DelimitedFiles

using KAPDQMC
using LinearAlgebra
using Random
rng=MersenneTwister(time_ns())

t=1;   Lattice="HoneyComb"    
U0=0;   Uend=5;     Δt=0.1;    Rate=1;
BatchSize=10;
  
path="E:/桌面/JuliaDQMC/code/QuenchSU2PQMC/TEST/"

L=3
site=[L,L]

model=Hubbard_Para(t,U0,Uend,Rate,Lattice,site,Δt,BatchSize)

s=Initial_s(model,rng)

Gt=Gτ(model,s,1)

Gt,G0,Gt0,G0t=G4(model,s,1,div(model.Nt,2))

println(model.Nt)
# Half
indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
# println(indexA)

# HalfHalf
indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))

λ=0.0
Nλ=2
# ss=ctrl_SCEEicr(path,model,indexA,indexB,1,λ,Nλ,[s[:,:],s[:,:]],false)
ss=ctrl_SCEEicr(path,model,indexA,indexB,1,λ,Nλ,[s[:,:],s[:,:]],true)



# s=phy_update(path,model,s,10,false)
# s=phy_update(path,model,s,500,true)


