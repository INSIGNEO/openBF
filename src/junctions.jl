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

# In `openBF` two types of junctions are provided: conjunctions
# <div style="text-align:center">
# <img src="images/conj.pdf.png" width="225"></div>
# and bifurcations.
# <div style="text-align:center">
# <img src="images/bif.pdf.png" width="225"></div>
#
# `joinVessels` function calls the right solver depending on the number of
# vessels involved in the junction.

# *function* __`joinVessels`__
#
# ----------------------------------------------------------------------------
# Parameters:
# ------------------- --------------------------------------------------------
# `b`                 `::Blood` model matrix row.
#
# `vessels...`        Collection of`::Vessel` data structures involved in the
#                     junction system.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# A collection containing two vessels refers to a conjunction and is solved
# by [`solveConjunction`](conjunctions.html#solveConjunction) function.
# When there are three vessels the junction is a bifurcation and it is solved
# by [`solveBifurcation`](bifurcations.html#solveBifurcation).
# ----------------------------------------------------------------------------
# <a name="joinVessels"></a>
function joinVessels(b :: Blood, vessels...)

  if length(vessels) == 2
    solveConjunction(b, vessels[1], vessels[2])

  elseif length(vessels) == 3
    solveBifurcation(vessels[1], vessels[2], vessels[3])

  end

end
