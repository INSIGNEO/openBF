#=
Copyright 2022 INSIGNEO Institute for in silico Medicine

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=#


__precompile__(true)
module openBF

    export Vessel, Heart, Blood, runSimulation
    using YAML
    using ArgParse
    using StaticArrays
    using DelimitedFiles
    using LinearAlgebra
    using Printf


    """
    Heart type
    """
    struct Heart
        inlet_type :: String
        cardiac_T :: Float64
        input_data :: Array{Float64,2}
        inlet_number :: Int
    end


    """
    Blood type
    """
    struct Blood
        mu  :: Float64
        rho :: Float64
        rho_inv :: Float64
    end


    """
    Vessel type
    """
    mutable struct Vessel
        label :: String

        #Topology
        ID :: Int
        sn :: Int
        tn :: Int

        #Inlet
        inlet :: Bool
        heart :: Heart

        #Numerical constants
        M       :: Int
        dx      :: Float64
        invDx   :: Float64
        halfDx  :: Float64

        #Physical constants
        beta            :: Vector{Float64}
        gamma           :: Vector{Float64}
        s_15_gamma      :: Vector{Float64}
        gamma_ghost     :: Vector{Float64}
        A0              :: Vector{Float64}
        s_A0            :: Vector{Float64}
        inv_A0          :: Vector{Float64}
        s_inv_A0        :: Vector{Float64}
        Pext            :: Float64
        viscT           :: Float64
        wallE           :: Vector{Float64}
        wallVa          :: Vector{Float64}
        wallVb          :: Vector{Float64}

        #Iterative solution
        A :: Vector{Float64}
        Q :: Vector{Float64}
        u :: Vector{Float64}
        c :: Vector{Float64}
        P :: Vector{Float64}

        A_t :: Array{Float64,2}
        Q_t :: Array{Float64,2}
        u_t :: Array{Float64,2}
        c_t :: Array{Float64,2}
        P_t :: Array{Float64,2}

        A_l :: Array{Float64,2}
        Q_l :: Array{Float64,2}
        u_l :: Array{Float64,2}
        c_l :: Array{Float64,2}
        P_l :: Array{Float64,2}

        #Riemann invariants
        W1M0 :: Float64
        W2M0 :: Float64

        #Ghost cells
        U00A :: Float64
        U00Q :: Float64
        U01A :: Float64
        U01Q :: Float64

        UM1A :: Float64
        UM1Q :: Float64
        UM2A :: Float64
        UM2Q :: Float64

        #Result file names
        last_P_name :: String
        last_Q_name :: String
        last_A_name :: String
        last_c_name :: String
        last_u_name :: String

        out_P_name :: String
        out_Q_name :: String
        out_A_name :: String
        out_c_name :: String
        out_u_name :: String

        #Saving nodes
        node2 :: Int
        node3 :: Int
        node4 :: Int

        #Peripheral boundary condition
        Rt :: Float64
        R1 :: Float64
        R2 :: Float64
        Cc :: Float64
        Pc :: Float64

        #Slope
        slope :: Vector{Float64}

        #MUSCLArrays
        fluxA :: Vector{Float64}
        fluxQ :: Vector{Float64}
        uStar :: Array{Float64,2}

        vA :: Vector{Float64}
        vQ :: Vector{Float64}

        dU :: Array{Float64,2}

        slopesA :: Vector{Float64}
        slopesQ :: Vector{Float64}

        Al :: Vector{Float64}
        Ar :: Vector{Float64}

        Ql :: Vector{Float64}
        Qr :: Vector{Float64}

        Fl :: Vector{Float64}
        Fr :: Vector{Float64}

        #Outlet type
        outlet :: String

        dTaudx :: Vector{Float64}
        dA0dx :: Vector{Float64}
    end

    include("initialise.jl")

    include("boundary_conditions.jl")
    include("solver.jl")

    include("junctions.jl")
    include("conjunctions.jl")
    include("bifurcations.jl")
    include("anastomosis.jl")

    include("IOutils.jl")
    include("check_convergence.jl")

    include("program.jl")
end
