load_yaml_config(yaml_config::String) = YAML.load_file(yaml_config)

function get_conv_error(n::Network)
    error, edge = findmax(get_conv_error, n.vessels)
    error, n.vessels[edge].label
end

function get_conv_error(v::Vessel)
    current = readdlm(v.label * "_P.temp")
    prev = readdlm(v.label * "_P.last")
    sqrt(sum((current[:, 4] - prev[:, 4]) .^ 2)/v.M) / 133.332
end

function step!(n::Network, dt::Float64, current_time::Float64)
    dt = calculate_Δt(n)
    solve!(n, dt, current_time)
    update_ghost_cells!(n)
end

getprog(cycle::Int64, verbose::Bool) =
    verbose ?
    ProgressUnknown(desc = "Solving cycle #$cycle:", spinner = true, showspeed = true) :
    nothing

function run_simulation(
    yaml_config::String;
    verbose::Bool = false,
    out_files::Bool = false,
    conv_ceil::Bool = false,
)
    verbose && println("Loading config...")
    config = load_yaml_config(yaml_config)

    project_name = config["project_name"]
    blood = Blood(config["blood"])

    # TODO: multiple inlets
    heart = Heart(config["project_name"])

    # TODO: results dir from config
    ###############################################
    results_dir = project_name * "_results"
    ~isdir(results_dir) && mkdir(results_dir)
    # TODO: handle absolute paths
    cp(yaml_config, joinpath(results_dir, yaml_config), force=true)
    cd(results_dir)
    ###############################################

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
    dt = calculate_Δt(network)

    passed_cycles = 0
    prog = getprog(passed_cycles, verbose)

    counter = 1
    conv_error = floatmax()
    @time while true
        step!(network, dt, current_time)

        if current_time >= checkpoints[counter]
            flush_to_temp(current_time, network, config)
            counter += 1
        end

        verbose && next!(prog)

        if (current_time - heart.cardiac_period * passed_cycles) >= heart.cardiac_period &&
           (current_time - heart.cardiac_period * passed_cycles + dt) > heart.cardiac_period

            # check convergence
            if passed_cycles > 0
                conv_error, error_loc = get_conv_error(network)
            end

            move_temp_to_last(network, config)
            out_files && append_last_to_out(network, config)

            passed_cycles > 0 &&
                verbose &&
                finish!(prog, showvalues = [("RMSE (mmHg)", conv_error), ("@", error_loc)])
            verbose && println()

            checkpoints = checkpoints .+ heart.cardiac_period
            passed_cycles += 1
            counter = 1
            prog = getprog(passed_cycles, verbose)
        end
        current_time += dt
        if current_time >= total_time ||
           passed_cycles == config["solver"]["cycles"] ||
           conv_error < config["solver"]["convergence_tolerance"]
            break
        end
    end
    cd("..")
end
