#=
Copyright 2020 INSIGNEO Institute for in silico Medicine

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


(VERSION >= v"0.4") && __precompile__(true)
module openBF

    export Vessel, Heart, Blood, runSimulation
    using YAML
    using ArgParse
    using StaticArrays
    using DelimitedFiles
    using LinearAlgebra
    using Printf
    # IMPORT PACKAGE FOR HYPERGEOMETRIC FUNCTIONS
    using HypergeometricFunctions
    import HypergeometricFunctions: drummond2F1


    """
    Heart type
    """
    mutable struct Heart
        inlet_type :: String
        cardiac_T :: Float64
        input_data :: Array{Float64,2}
        inlet_number :: Int
    end


    """
    Blood type
    """
    mutable struct Blood
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
        beta            :: Array{Float64,1}
        gamma           :: Array{Float64,1}
        s_15_gamma      :: Array{Float64,1}
        gamma_ghost     :: Array{Float64,1}
        A0              :: Array{Float64,1}
        s_A0            :: Array{Float64,1}
        inv_A0          :: Array{Float64,1}
        s_inv_A0        :: Array{Float64,1}
        Pext            :: Float64
        viscT           :: Float64
        wallE           :: Array{Float64,1}
        wallVa          :: Array{Float64,1}
        wallVb          :: Array{Float64,1}

        #Iterative solution
        A :: Array{Float64,1}
        Q :: Array{Float64,1}
        u :: Array{Float64,1}
        c :: Array{Float64,1}
        P :: Array{Float64,1}

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
        slope :: Array{Float64,1}

        #MUSCLArrays
        flux :: Array{Float64,2}
        uStar :: Array{Float64,2}

        vA :: Array{Float64,1}
        vQ :: Array{Float64,1}

        dU :: Array{Float64,2}

        slopesA :: Array{Float64,1}
        slopesQ :: Array{Float64,1}

        Al :: Array{Float64,1}
        Ar :: Array{Float64,1}

        Ql :: Array{Float64,1}
        Qr :: Array{Float64,1}

        Fl :: Array{Float64,2}
        Fr :: Array{Float64,2}

        # Catheter
        Ac :: Float64               # Area of the catheter
        corrRI :: Float64           # Correction factor to Riemann Invariants

        #Outlet type
        outlet :: String
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
