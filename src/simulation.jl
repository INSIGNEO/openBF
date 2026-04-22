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


# -------------------- deprecated ???? -----------------------------------------
function compute_pressure(n::Network)
    compute_pressure.(n.vessels_vec, Ref(n.blood.rho), Ref("temp"))
    compute_pressure.(n.vessels_vec, Ref(n.blood.rho), Ref("last"))
end
function compute_pressure(v::Vessel, rho::Float64, ext::String)
    A = readdlm(v.label * "_A.$(ext)")

    P = zeros(Float64, size(A))
    P[:,1] = A[:,1]

    for (i, j) in enumerate([1, v.node2, v.node3, v.node4, v.M])
        @. P[:,i+1] = pressure(A[:,i+1], v.A0[j], v.beta[j], v.Pext)

        if v.viscoelastic
            P[:,i+1] .+= v.Cv[j] .* rho ./ A[:,i+1] .* diff([A[:,i+1];A[end,i+1]])./diff([A[:,1];A[end,1]+(A[end,1]-A[end-1,1])])
        end
    end

    open(v.label * "_P.$(ext)", "w") do io
        writedlm(io, P)
    end
end
# ------------------------------------------------------------------------------

load_yaml_config(yaml_config::String) = YAML.load_file(yaml_config)

function get_conv_error(n::Network)
    error, idx = findmax(get_conv_error, n.vessels_vec)
    error, n.vessels_vec[idx].label
end

function get_conv_error(v::Vessel)
    current = @view v.waveforms["P"][:, 4]
    prev    = @view v.waveforms_prev["P"][:, 4]
    sqrt(mean(abs2.(current .- prev))) / 133.332
end

function _conv_error_field(n::Network, field::String)
    err = 0.0
    for v in n.vessels_vec
        haskey(v.waveforms, field) || continue
        current = @view v.waveforms[field][:, 4]
        prev    = @view v.waveforms_prev[field][:, 4]
        e = sqrt(mean(abs2.(current .- prev)))
        e > err && (err = e)
    end
    err
end

getprog(cycle::Int64, verbose::Bool) =
    verbose ?
    ProgressUnknown(desc = "Solving cycle #$cycle:", spinner = true, showspeed = true) :
    nothing
function get_inlet_file(config)::String
    if ~haskey(config, "inlet_file")
        return config["project_name"]::String * "_inlet.dat"
    end
    config["inlet_file"]
end

function preamble(yaml_config, verbose, savedir)

    verbose && println("Loading config...")
    config = load_yaml_config(yaml_config)

    project_name = config["project_name"]
    verbose && println("project name: $project_name")

    if savedir == ""
        results_dir = get(config, "output_directory", project_name * "_results")
    else
        results_dir = savedir
    end
    isdir(results_dir) && rm(results_dir, recursive = true)
    ~isdir(results_dir) && mkpath(results_dir)
    yaml_config_name = last(splitpath(yaml_config))
    cp(yaml_config, joinpath(results_dir, yaml_config_name), force = true)

    inlet_file = get_inlet_file(config)
    try
        cp(inlet_file, joinpath(results_dir, inlet_file), force = true)
    catch
        cp(joinpath([splitpath(yaml_config)[1:end-1]; inlet_file]), joinpath(results_dir, inlet_file), force=true)
    end
    cd(results_dir)

    config
end

function run_simulation(
    yaml_config::String;
    verbose::Bool = true,
    out_files::Bool = false,
    save_stats::Bool = false,
    savedir::String = "",
    no_tui::Bool = false,
    _observer = nothing,
)
    initial_dir = pwd()
    config::Dict{String, Any} = preamble(yaml_config, verbose, savedir)

    blood = Blood(config["blood"])
    heart = Heart(get_inlet_file(config)::String)

    tempsave::Vector{String} = deepcopy(get(config, "write_results", []))
    "P" ∉ tempsave && push!(tempsave, "P")
    network = Network(
        config["network"],
        blood,
        heart,
        config["solver"]["Ccfl"]::Float64,
        config["solver"]["jump"]::Int64,
        tempsave,
        verbose = verbose,
    )
    total_time = float(config["solver"]["cycles"]) * heart.cardiac_period
    jump::Int64 = config["solver"]["jump"]
    checkpoints = range(0, stop = heart.cardiac_period, length = jump)

    verbose && println("\nStart simulation")

    observer = if _observer !== nothing
        _observer
    elseif tui_should_run(!no_tui)
        TUIObserver(network.vessels_vec, network.is_outlet;
                    sim_name = config["project_name"],
                    t_final  = total_time,
                    n_cycles = config["solver"]["cycles"])
    else
        nothing
    end
    observer !== nothing && Base.exit_on_sigint(false)
    start_render!(observer)

    # suppress ProgressMeter when TUI is active (Phase 8 routes output to log pane)
    loop_verbose = verbose && (observer === nothing)

    current_time = 0.0
    passed_cycles = 0
    prog = getprog(passed_cycles, loop_verbose)
    counter = 1
    counter_prog = 0
    conv_error::Float64 = floatmax()
    converged = false

    try
        stats = @timed while true

            # step
            dt = calculateΔt(network)
            solve!(network, dt, current_time)
            update_ghost_cells!(network)
            record_step!(observer, network.vessels_vec, current_time, dt)
            loop_verbose && (counter_prog += 1) % 100 == 0 && next!(prog)

            if current_time >= checkpoints[counter]
                save_waveforms.(counter, current_time, network.vessels_vec)
                counter += 1
            end

            # at the end of the cardiac cycle
            if (current_time - heart.cardiac_period * passed_cycles) >= heart.cardiac_period &&
               (current_time - heart.cardiac_period * passed_cycles + dt) > heart.cardiac_period

                # check convergence
                conv_error_A = 0.0
                conv_error_Q = 0.0
                if passed_cycles > 0
                    # compute_pressure(network)
                    conv_error, error_loc = get_conv_error(network)
                    conv_error_A = _conv_error_field(network, "A")
                    conv_error_Q = _conv_error_field(network, "Q")
                end
                record_cycle!(observer, passed_cycles, conv_error, conv_error_A, conv_error_Q)

                flush_waveforms.(network.vessels_vec)
                swap_waveforms!.(network.vessels_vec)
                out_files && append_last_to_out.(network.vessels_vec)

                if loop_verbose
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
                prog = getprog(passed_cycles, loop_verbose)
                counter = 1
            end

            # at the end of the simulation
            if current_time >= total_time ||
               passed_cycles == config["solver"]["cycles"] ||
               conv_error < config["solver"]["convergence_tolerance"]
                loop_verbose && finish!(prog)
                converged = true
                break
            end
            current_time += dt
        end
        stop_render!(observer)
        observer !== nothing && Base.exit_on_sigint(true)
        tokeep::Vector{String} = get(config, "write_results", [])
        cleanup.(network.vessels_vec, Ref(tokeep))

        verbose && printstats(stats, converged, passed_cycles)
        save_stats && savestats(stats, converged, passed_cycles, config["project_name"])
    catch e
        stop_render!(observer)
        observer !== nothing && Base.exit_on_sigint(true)
        println("\nAn error occurred, terminating simulation.\n")
        println(e)
    end

    cd(initial_dir)
end

function cleanup(v::Vessel, tokeep)
    ~v.tosave && rm.(glob("$(v.label)_*"))
    for l in ("P", "A", "Q", "u")
        l ∉ tokeep && rm.(glob("$(v.label)_$l.*"), force=true)
    end
end

function printstats(s, converged, passed_cycles)
    if converged
        println("Simulation converged after $passed_cycles cardiac cycles")
    else
        println("Simulation not converged, stopping after $passed_cycles cardiac cycles")
    end
    @printf "Elapsed time: %5.2fs\n" s.time
    @printf "Memory: %5.2fGB\n" s.bytes*1e-9
    gctime = @sprintf "GC time: %5.2f" s.gctime
    println("$(gctime)%")
end

function savestats(s, converged, passed_cycles, simulation_name)
    open("$simulation_name.conv", "w") do io
        println(io, "$converged")
        println(io, "$passed_cycles")
        println(io, "$(s.time)")
        println(io, "$(s.bytes)")
        println(io, "$(s.gctime)")
    end
end
