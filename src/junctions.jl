#= Copyright (C) 2017 Alessandro Melis.

  This file is part of openBF.

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with openBF.  If not, see <http://www.gnu.org/licenses/>.
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
