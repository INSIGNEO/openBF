#=
Copyright 2018 INSIGNEO Institute for in silico Medicine

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

function joinVessels(b :: Blood, vessels...)
  if length(vessels) == 2
    solveConjunction(b, vessels[1], vessels[2])

  elseif length(vessels) == 3
    solveBifurcation(vessels[1], vessels[2], vessels[3])
  end
end

function newtonRaphson(J :: Array{Float64,2}, v1 :: Vessel, v2 :: Vessel, v3 :: Vessel,
                       U :: Array, k :: Array, funW, funF)
    W = funW(U, k)
    F = funF(v1, v2, v3, U, k, W)

    nr_toll_U = 1e-5
    nr_toll_F = 1e-5

    while true
      dU = J\(-F)
      U_new = U + dU

      if any(isnan(dot(F,F)))
        # println(F)
        # @printf "error at bifurcation with vessels %s, %s, and %s \n" v1.label v2.label v3.label
        error("Newton-Raphson doens't converge!")
      end

      u_ok = 0
      f_ok = 0
      for i in 1:length(dU)
        if abs(dU[i]) <= nr_toll_U || abs(F[i]) <= nr_toll_F
          u_ok += 1
          f_ok += 1
        end
      end

      if u_ok == length(dU) || f_ok == length(dU)
        return U_new
      else
        U = U_new
        W = calculateWstarBifurcation(U, k)
        F = calculateFofUBifurcation(v1, v2, v3, U, k, W)
      end
    end
end
