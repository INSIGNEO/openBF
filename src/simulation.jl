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

function run_simulation(yaml_config::String;
    verbose::Bool = true,
    out_files::Bool = false,
    conv_ceil::Bool = false,
)
    verbose && println("Loading config...")
    config = load_yaml_config(yaml_config)

    project_name = config["project_name"]
    blood = Blood(config["blood"])
    heart = Heart(config["project_name"])

    # TODO: results dir from config
    ########################################
    results_dir = project_name * "_results"
    isdir(results_dir) && rm(results_dir, recursive=true)
    ~isdir(results_dir) && mkdir(results_dir)
    cp(yaml_config, joinpath(results_dir, yaml_config), force=true)
    cd(results_dir)
    #########################################

    network = Network(config["network"], blood, heart, config["solver"]["Ccfl"], verbose=verbose)
    total_time = config["solver"]["cycles"] * heart.cardiac_period
    jump = config["solver"]["jump"]
    checkpoints = range(0, stop = heart.cardiac_period, length = jump)

    verbose && println("\nStart simulation")
    # TODO: log stats

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
            flush_to_temp(current_time, network)
            counter += 1
        end

        # at the end of the cardiac cycle
        if (current_time - heart.cardiac_period * passed_cycles) >= heart.cardiac_period &&
           (current_time - heart.cardiac_period * passed_cycles + dt) > heart.cardiac_period

            # check convergence
            if passed_cycles > 0
                conv_error, error_loc = get_conv_error(network)
            end

            move_temp_to_last(network)
            out_files && append_last_to_out(network)

            if verbose
                if passed_cycles > 0
                    finish!(prog, showvalues = [("RMSE (mmHg)", conv_error), ("@", error_loc)])
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
            finish!(prog)
            break
        end
        current_time += dt
    end
    cd("..")
end
