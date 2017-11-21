#=
Copyright 2017 INSIGNEO Institute for in silico Medicine

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

# Convergence is checked starting from the simulated third cardiac cycle.
# Waveforms are first loaded and subdivided depending on the number of cycles
# already passed. This task is accomplished by function
# [`extractWaveform`](check_convergence.html#extractWaveform). Waveforms are
# then resampled in order to contain a smaller number of points. This would
# speed up the following error calculation. These two steps are taken by
# [`resampleAndCalculateError`](check_convergence.html#resampleAndCalculateError)
# function. Loading, re-sampling, and error calculation operations are done
# for all the quantities calculated by `openBF`. All these functions are
# wrapped into
# [`checkAllQuantities`](check_convergence.html#checkAllQuantities) which is
# called inside [main.jl](main.html#check_convergence) loop.

# *function* __`checkAllQuantities`__ $\rightarrow$ `max_error::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# --------------- ------------------------------------------------------------
# `v`             `::Vessel` vessel of which the error is being checked.
#
# `passed_cycles` `::Int` number of cardiac cycles already passed.
#
# `n_pts`         `::Int` number of points in the resampled waveform.
# ----------------------------------------------------------------------------
# <a name="checkAllQuantities"></a>
function checkAllQuantities(v             :: Vessel,
                            passed_cycles :: Int64,
                            n_pts         :: Int64)
  # --------------------------------------------------------------------------
  # Functioning
  # --------------------------------------------------------------------------
  # `er` is an array which will contain the percentage error for each quantity
  # in the current vessel. `qs` collection contains all the suffixes needed to
  # build `.out` filenames for each quantity.
  # --------------------------------------------------------------------------
  er = zeros(Float64, 5)

  qs = ["_P.out", "_Q.out", "_A.out", "_c.out", "_u.out"]
  # For each quantity `q` listed in `qs`, the `.out` filename is created by
  # joining the `label` inside `v` and the suffix `q`. Then the desired
  # waveforms are loaded by
  # [`extractWaveform`](check_convergence.html#extractWaveform), re-sampled,
  # and the error is calculated by
  # [`resampleAndCalculateError`](check_convergence.html#resampleAndCalculateError).
  # Eventually, the maximum error scored among all the quantities is returned.
  #
  # --------------------------------------------------------------------------
  # Returns:
  # ------------- ------------------------------------------------------------
  # `max_error`   `::Float` maximum error within all the quantities in the
  #               current vessel.
  # --------------------------------------------------------------------------
  idx = 1
  for q in qs
    filename = join([v.label, q])

    W = extractWaveform(filename, passed_cycles)

    er[idx] = resampleAndCalculateError(v, W, passed_cycles, n_pts)

    idx += 1
  end

  return maximum(er)
end

# *function* __`checkAllQuantities`__ $\rightarrow$ `max_error::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# --------------- ------------------------------------------------------------
# `vs`             `::Array{Vessel, 1}` collection of vessels.
#
# `passed_cycles` `::Int` number of cardiac cycles already passed.
#
# `n_pts`         `::Int` number of points in the resampled waveform.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# To compute the error for all the vessels, `checkAllQuantities` is
# re-defined to handle also a vector of `Vessel`s. `checkAllQuantities` is
# called recursively.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ------------- --------------------------------------------------------------
# `max_error`   `::Float` maximum error among all the vessels.
# ----------------------------------------------------------------------------
# <a name="checkAllQuantities2"></a>
function checkAllQuantities(vs            :: Array{Vessel, 1},
                            passed_cycles :: Int64,
                            n_pts         :: Int64)

  ers = zeros(Float64, length(vs))

  idx = 1
  for v in vs

    ers[idx] = checkAllQuantities(v, passed_cycles, n_pts)
    idx += 1

  end

  return maximum(ers)
end

# *function* __`extractWaveform`__ $\rightarrow$ `Fx::Array{Array{Float, 2}}`
#
# ----------------------------------------------------------------------------
# Parameters:
# ------------- --------------------------------------------------------------
# `filename`    `::String` current time in seconds.
#
# `cycles`      `::Int` $\Delta t$, current time step.
# ----------------------------------------------------------------------------
# <a name="extractWaveform"></a>
function extractWaveform(filename :: String, cycles :: Int64)
  # --------------------------------------------------------------------------
  # Functioning
  # --------------------------------------------------------------------------
  # `filename` is a string referring to the `.out` file containing the
  # waveform to be imported. This is read by `readdlm` function and stored
  # into
  # `F` array. `F` is a timesteps $\times$ 6 array; the first column contains
  # the time variable, the remaining columns are one for each monitored node
  # along the current vessel.
  # --------------------------------------------------------------------------
  F = readdlm(filename)
  # These three lines find the cardiac period `t` by knowing how many
  # cycles have already passed.
  sf1 = int(length(F[:,1])*(cycles-1)/cycles)
  T1 = F[sf1:end,1] .- F[sf1,1]
  t = T1[end]
  # `Fx` is an empty array which will contain one array for each cardiac
  # cycle.
  Fx = Array{Float64, 2}[]
  # This loop skims through `F`, separates waveforms, and add them to `Fx`.
  # `si` is the waveform source index, the index from which the waveform
  # starts. `ti` is the waveform terminal index, the index at which the
  # waveforms ends because a new cardiac cycle is due to start.
  si = 1
  ti = 1
  idx = 1
  for i in 1:length(F[:,1])-1
    if F[i,1] < t*idx && F[i+1, 1] >= t*idx
      ti = i
      push!(Fx, F[si:ti,:])
      si = ti
      idx += 1
    end
  end
  # --------------------------------------------------------------------------
  # Returns:
  # ------------- ------------------------------------------------------------
  # `Fx`           `::Array{Array{Float, 2}}`
  # --------------------------------------------------------------------------
  return Fx
end

# *function* __`resampleAndCalculateError`__ $\rightarrow$
# `max_error::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ----------------- ----------------------------------------------------------
# `W`               `::Array{Array{Float, 2}}` collection of all waveforms
#                   during the passed cardiac cycles for the current
#                   quantity.
#
# `passed_cycles`   `::Int` number of cardiac cycles already simulated.
#
# `n_pts`           `::Int` number of points to be contained in the
#                   re-sampled waveform.
# ----------------------------------------------------------------------------
# <a name="resampleAndCalculateError"></a>
function resampleAndCalculateError(v, W :: Array{Array{Float64, 2}, 1},
                                   passed_cycles :: Int64,
                                   n_pts         :: Int64)
  # --------------------------------------------------------------------------
  # Functioning
  # --------------------------------------------------------------------------
  # `F` will contain the re-sampled waveform, `er` is an empty variable
  # for the maximum error.
  # --------------------------------------------------------------------------
  F = zeros(Float64, n_pts, 2)
  er = 0
  # The error is computed for all the monitor nodes. Their values are stored
  # in `W`'s columns from 2 to 6. Each waveform is re-sampled by using
  # [`resample`](http://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.resample.html#scipy.signal.resample)
  # function from `scipy.signal` library. In order to compute the error, only
  # the last and the previous waveforms are analysed.
  for j in 2:6

    F[:,1], t = signal.resample(W[end-1][:,j], n_pts, W[end-1][:,1])
    F[:,2], t = signal.resample(W[end  ][:,j], n_pts, W[end  ][:,1])

    er = abs(F[:,1] .- F[:,2])./F[:,1]*100
  end
  # --------------------------------------------------------------------------
  # Returns:
  # ----------- --------------------------------------------------------------
  # `max_error` `::Float` maximum error on the current quantity along the
  #             the current vessel.
  # --------------------------------------------------------------------------
  return maximum(er)
end

function checkConvergence(vessels :: Array{Vessel, 1},
  convergence_tollerance :: Float64)

  for v in vessels
    err = checkConvergence(v, convergence_tollerance)
    if err > convergence_tollerance
      return err
    end
  end
  return err
end

function checkConvergence(v :: Vessel,
  convergence_tollerance :: Float64)

  qs = ["_P", "_Q"]

  for q in qs
    filename_last = join([v.label, q, ".last"])
    filename_temp = join([v.label, q, ".temp"])

    err = checkQuantity(filename_last, filename_temp, convergence_tollerance)
    if err > convergence_tollerance
      return err
    end
  end
  return err
end

function checkQuantity(last, temp, convergence_tollerance)

  s1 = readdlm(last)
  s2 = readdlm(temp)
  n_pts = 1000

  r1 = n_pts/length(s1[:,2])
  r2 = n_pts/length(s2[:,2])
  for j in 2:5
    sr1 = resample(s1[:,j], r1)
    sr2 = resample(s2[:,j], r2)
    for i in 20:length(sr1)-20
      err = abs(100*(sr2[i]-sr1[i])/sr1[i])
      if err > convergence_tollerance
        # println("\n",err)
        return err
      end
    end
  end
  return 9.99
end

function checkConvergence(edge_list, vessels :: Array{Vessel, 1}, passed_cycles :: Int64)
    qs = "_A"
    err = zeros(size(edge_list)[1]).+100
    if passed_cycles >= 2

        for i in 1:size(edge_list)[1]
            v = vessels[i]
            lbl = v.label

            filename_last = join([lbl, qs, ".last"])
            filename_temp = join([lbl, qs, ".temp"])
            w_last = readdlm(filename_last)
            w_temp = readdlm(filename_temp)

            if length(w_last[:,1]) == length(w_temp[:,1])
                err[i] = maximum(abs.((w_last[2:end,:].-w_temp[2:end,:])./w_last[2:end,:])*100)
            end

        end
    end
    return maximum(err)
end
