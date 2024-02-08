#=
Copyright 2015-2024 INSIGNEO Institute for in silico Medicine

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

function compute_pressure(n::Network)
    compute_pressure.(values(n.vessels), Ref(n.blood.rho), Ref("temp"))
    compute_pressure.(values(n.vessels), Ref(n.blood.rho), Ref("last"))
end
function compute_pressure(v::Vessel, rho::Float64, ext::String)
    A = readdlm(v.label * "_A.$(ext)")

    P = zeros(Float64, size(A))
    P[:,1] = A[:,1]

    for (i, j) in enumerate([1, v.node2, v.node3, v.node4, v.M])
        @. P[:,i+1] = pressure(A[:,i+1], v.A0[j], v.beta[j], v.Pext)

        if v.Cv[1] != 0.0
            P[:,i+1] .+= v.Cv[j] .* rho ./ A[:,i+1] .* diff([A[:,i+1];A[end,i+1]])./diff([A[:,1];A[end,1]+(A[end,1]-A[end-1,1])])
        end
    end

    open(v.label * "_P.$(ext)", "w") do io
        writedlm(io, P)
    end
end

load_yaml_config(yaml_config::String) = YAML.load_file(yaml_config)

function get_conv_error(n::Network)
    error, edge = findmax(get_conv_error, n.vessels)
    error, n.vessels[edge].label
end

function get_conv_error(v::Vessel)
    current = readdlm(v.label * "_P.temp")
    prev = readdlm(v.label * "_P.last")
    sqrt(mean(current[:, 4] - prev[:, 4]) .^ 2) / 133.332
end

getprog(cycle::Int64, verbose::Bool) =
    verbose ?
    ProgressUnknown(desc = "Solving cycle #$cycle:", spinner = true, showspeed = true) :
    nothing
function get_inlet_file(config)
    if ~haskey(config, "inlet_file")
        return config["project_name"]*"_inlet.dat"
    end
    config["inlet_file"]
end
function run_simulation(
    yaml_config::String;
    verbose::Bool = true,
    out_files::Bool = false,
    conv_ceil::Bool = false,
)
    verbose && println("Loading config...")
    config = load_yaml_config(yaml_config)

    tosave = config["write_results"]
    isempty(tosave) && (tosave = ["P", "A"])
    "P" ∉ tosave && push!(tosave, "P")
    "A" ∉ tosave && push!(tosave, "A")

    verbose && println("project name: $(config["project_name"])")

    project_name = config["project_name"]
    blood = Blood(config["blood"])
    heart = Heart(get_inlet_file(config))

    # TODO: results dir from config
    ########################################
    results_dir = project_name * "_results"
    isdir(results_dir) && rm(results_dir, recursive = true)
    ~isdir(results_dir) && mkdir(results_dir)
    cp(yaml_config, joinpath(results_dir, yaml_config), force = true)
    cd(results_dir)
    #########################################

    network = Network(
        config["network"],
        blood,
        heart,
        config["solver"]["Ccfl"],
        verbose = verbose,
    )
    total_time = config["solver"]["cycles"] * heart.cardiac_period
    jump = config["solver"]["jump"]
    checkpoints = range(0, stop = heart.cardiac_period, length = jump)

    verbose && println("\nStart simulation")

    current_time = 0.0
    passed_cycles = 0
    prog = getprog(passed_cycles, verbose)
    counter = 1
    conv_error = floatmax()
    @time while true

        # step
        dt = calculate_Δt(network)
        solve!(network, dt, current_time)
        update_ghost_cells!(network)
        verbose && next!(prog)

        if current_time >= checkpoints[counter]
            flush_to_temp(current_time, network, tosave)
            counter += 1
        end

        # at the end of the cardiac cycle
        if (current_time - heart.cardiac_period * passed_cycles) >= heart.cardiac_period &&
           (current_time - heart.cardiac_period * passed_cycles + dt) > heart.cardiac_period

            # check convergence
            if passed_cycles > 0
                compute_pressure(network)
                conv_error, error_loc = get_conv_error(network)
            end

            move_temp_to_last(network, tosave)
            out_files && append_last_to_out(network, tosave)

            if verbose
                if passed_cycles > 0
                    finish!(
                        prog,
                        showvalues = [("RMSE (mmHg)", conv_error), ("@", error_loc)],
                    )
                else
                    finish!(prog)
                end
                println()
            end

            checkpoints = checkpoints .+ heart.cardiac_period
            passed_cycles += 1
            prog = getprog(passed_cycles, verbose)
            counter = 1
        end

        # at the end of the simulation
        if current_time >= total_time ||
           passed_cycles == config["solver"]["cycles"] ||
           conv_error < config["solver"]["convergence_tolerance"]
            verbose && finish!(prog)
            break
        end
        current_time += dt
    end
    cd("..")
end
