const _SPINNER = ('|', '/', '-', '\\')

function _make_header(obs::TUIObserver, term_cols::Int)
    s       = obs.scalars
    spinner = _SPINNER[(obs.frame_count % 4) + 1]
    elapsed = time() - obs.start_time
    n       = length(obs.waveforms)
    idx     = obs.active_site[]

    content = @sprintf(
        " %c  %-18s  t=%7.4f/%-7.4fs  cycle %2d/%-2d  wall %5.1fs  %6.0f stp/s  [%d/%d] ",
        spinner, obs.sim_name, s.t, obs.t_final,
        s.passed_cycles, obs.n_cycles,
        elapsed, obs.steps_per_sec, idx, n)

    Panel(content; width = term_cols, padding = (0,0,0,0), style = "blue dim")
end
