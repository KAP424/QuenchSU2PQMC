module KAPDQMC
    using Base.Filesystem
    using LinearAlgebra
    using DelimitedFiles
    using Random
    using Statistics

    include("Geometry.jl")
    export K_Matrix,xy_i,i_xy,area_index

    include("Hubbard_model.jl")
    export Hubbard_Para,_Hubbard_Para

    include("GreenMatrix.jl")
    export Gτ,G4,Initial_s,GroverMatrix
    # ,G12FF

    # include("phy_measure.jl")
    # export EK,NN,Magnetism,CzzofSpin

    # include("phy_update.jl")
    # export phy_update

    # include("EE_update.jl")
    # export ctrl_EEicr
    # ,EE_dir,EEICR

    include("SCEE.jl")
    export ctrl_SCEEicr
    
end
