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

# These functions handle I/O for result files.

# *function* __`openTempFiles`__
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------- -------------------------------------------------------------
# `v`            `::Vessel` vessel data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning:
# ----------------------------------------------------------------------------
# `.temp` files for the current vessel are `open`ed in the current directory
# with `w`riting permits.
# ----------------------------------------------------------------------------
# <a name="openTempFiles"></a>
function openTempFiles(v :: Vessel)

	v.temp_P      = open(v.temp_P_name,      "w")
	v.temp_Q      = open(v.temp_Q_name,      "w")
	v.temp_A      = open(v.temp_A_name,      "w")
	v.temp_c      = open(v.temp_c_name,      "w")
	v.temp_u      = open(v.temp_u_name,      "w")

end

# *function* __`openTempFiles`__
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------- -------------------------------------------------------------
# `vessels`      `::Array{Vessel, 1}` collection of vessel data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning:
# ----------------------------------------------------------------------------
# `.temp` files are `open`ed in the current directory with `w`riting permits
# for all the vessels in the array provided.
# ----------------------------------------------------------------------------
# <a name="openTempFiles2"></a>
function openTempFiles(vessels :: Array{Vessel, 1})

	for v in vessels
		openTempFiles(v)
	end

end

# *function* __`closeTempFiles`__
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------- -------------------------------------------------------------
# `v`            `::Vessel` vessel data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning:
# ----------------------------------------------------------------------------
# `.temp` files for the current vessel are `close`d in the current directory.
# ----------------------------------------------------------------------------
# <a name="closeTempFiles"></a>
function closeTempFiles(v :: Vessel)

	close(v.temp_P)
	close(v.temp_Q)
	close(v.temp_A)
	close(v.temp_c)
	close(v.temp_u)

end

# *function* __`closeTempFiles`__
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------- -------------------------------------------------------------
# `vessels`      `::Array{Vessel, 1}` collection of vessel data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning:
# ----------------------------------------------------------------------------
# `.temp` files are `close`d in the current directory
# for all the vessels in the array provided.
# ----------------------------------------------------------------------------
# <a name="closeTempFiles2"></a>
function closeTempFiles(vessels :: Array{Vessel, 1})

	for v in vessels
		closeTempFiles(v)
	end

end

function openCloseLastFiles(v :: Vessel)

	v.last_P = open(v.last_P_name, "w")
	v.last_Q = open(v.last_Q_name, "w")
	v.last_A = open(v.last_A_name, "w")
	v.last_c = open(v.last_c_name, "w")
	v.last_u = open(v.last_u_name, "w")

	close(v.last_P)
	close(v.last_Q)
	close(v.last_A)
	close(v.last_c)
	close(v.last_u)

end

function openCloseLastFiles(vessels :: Array{Vessel, 1})

	for v in vessels
		openCloseLastFiles(v)
	end

end

# *function* __`saveTempData`__
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------- -------------------------------------------------------------
# `t`            `::Float` current simulation time.
#
# `v`            `::Vessel` vessel data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning:
# ----------------------------------------------------------------------------
# `.temp` files are written with informations contained in the vessel data
# structure. For each quantity, 5 nodes are saved: the first, the middle, the
# last, and two intermediate nodes in between those three.
# ----------------------------------------------------------------------------
# <a name="saveTempData"></a>
# function saveTempData(t :: Float64, v :: Vessel)

# 	write(v.temp_P, t, " ", v.P[1], " ", v.P[v.node2], " ", v.P[v.node3], " ", v.P[v.node4], " ", v.P[end], "\n")
# 	write(v.temp_Q, t, " ", v.Q[1], " ", v.Q[v.node2], " ", v.Q[v.node3], " ", v.Q[v.node4], " ", v.Q[end], "\n")
# 	write(v.temp_A, t, " ", v.A[1], " ", v.A[v.node2], " ", v.A[v.node3], " ", v.A[v.node4], " ", v.A[end], "\n")
# 	write(v.temp_c, t, " ", v.c[1], " ", v.c[v.node2], " ", v.c[v.node3], " ", v.c[v.node4], " ", v.c[end], "\n")
# 	write(v.temp_u, t, " ", v.u[1], " ", v.u[v.node2], " ", v.u[v.node3], " ", v.u[v.node4], " ", v.u[end], "\n")

# 	# write(v.temp_Q, join((t, v.Q[1], v.Q[v.node2], v.Q[v.node3],
# 	# 												 v.Q[v.node4], v.Q[end]), " "), "\n")

# 	# write(v.temp_A, join((t, v.A[1], v.A[v.node2], v.A[v.node3],
# 	# 												 v.A[v.node4], v.A[end]), " "), "\n")

# 	# write(v.temp_c, join((t, v.c[1], v.c[v.node2], v.c[v.node3],
# 	# 												 v.c[v.node4], v.c[end]), " "), "\n")

# 	# write(v.temp_u, join((t, v.u[1], v.u[v.node2], v.u[v.node3],
# 	# 												 v.u[v.node4], v.u[end]), " "), "\n")

# end

function saveTempData(t :: Float64, v :: Vessel)

	println(v.temp_P, t, " ", v.P[1], " ", v.P[v.node2], " ", v.P[v.node3], " ", v.P[v.node4], " ", v.P[end])
	println(v.temp_A, t, " ", v.A[1], " ", v.A[v.node2], " ", v.A[v.node3], " ", v.A[v.node4], " ", v.A[end])
	println(v.temp_Q, t, " ", v.Q[1], " ", v.Q[v.node2], " ", v.Q[v.node3], " ", v.Q[v.node4], " ", v.Q[end])
	println(v.temp_c, t, " ", v.c[1], " ", v.c[v.node2], " ", v.c[v.node3], " ", v.c[v.node4], " ", v.c[end])
	println(v.temp_u, t, " ", v.u[1], " ", v.u[v.node2], " ", v.u[v.node3], " ", v.u[v.node4], " ", v.u[end])

end

# *function* __`saveTempData`__
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------- -------------------------------------------------------------
# `t`            `::Float` current simulation time.
#
# `vessels`      `::Array{Vessel, 1}` array containing vessel data structures.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning:
# ----------------------------------------------------------------------------
# `.temp` files are written with informations contained in the vessel data
# structure for all the vessels provided in the array `vessels`.
# ----------------------------------------------------------------------------
# <a name="saveTempData2"></a>
function saveTempData(t :: Float64, vessels :: Array{Vessel, 1})

	for v in vessels
		saveTempData(t, v)
	end

end

# *function* __`transferTempToOut`__
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------- -------------------------------------------------------------
# `v`            `::Vessel` vessel data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning:
# ----------------------------------------------------------------------------
# `.temp` files are concatenated to `.out` files by `appender.sh` script.
# ----------------------------------------------------------------------------
# <a name="transferTempToOut"></a>
function transferTempToOut(v :: Vessel)

	tempP = v.temp_P_name
	tempQ = v.temp_Q_name
	tempA = v.temp_A_name
	tempc = v.temp_c_name
	tempu = v.temp_u_name

	outP  = v.out_P_name
	outQ  = v.out_Q_name
	outA  = v.out_A_name
	outc  = v.out_c_name
	outu  = v.out_u_name

	run(`sh appender.sh $tempP $outP`)
	run(`sh appender.sh $tempQ $outQ`)
	run(`sh appender.sh $tempA $outA`)
	run(`sh appender.sh $tempc $outc`)
	run(`sh appender.sh $tempu $outu`)

end

function transferLastToOut(v :: Vessel)

	lastP = v.last_P_name
	lastQ = v.last_Q_name
	lastA = v.last_A_name
	lastc = v.last_c_name
	lastu = v.last_u_name

	outP  = v.out_P_name
	outQ  = v.out_Q_name
	outA  = v.out_A_name
	outc  = v.out_c_name
	outu  = v.out_u_name

	run(`sh appender.sh $lastP $outP`)
	run(`sh appender.sh $lastQ $outQ`)
	run(`sh appender.sh $lastA $outA`)
	run(`sh appender.sh $lastc $outc`)
	run(`sh appender.sh $lastu $outu`)

end

# *function* __`transferTempToOut`__
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------- -------------------------------------------------------------
# `vessels`      `::Array{Vessel, 1}` array containing vessel data structures.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning:
# ----------------------------------------------------------------------------
# `.temp` files are concatenated to `.out` files by `appender.sh` script for
# each vessel in `vessels`.
# ----------------------------------------------------------------------------
# <a name="transferTempToOut2"></a>
function transferTempToOut(vessels :: Array{Vessel, 1})

	for v in vessels
		transferTempToOut(v)
	end

end

function transferLastToOut(vessels :: Array{Vessel, 1})

	for v in vessels
		transferLastToOut(v)
	end

end

function transferTempToLast(v :: Vessel)

	tempP = v.temp_P_name
	tempQ = v.temp_Q_name
	tempA = v.temp_A_name
	tempc = v.temp_c_name
	tempu = v.temp_u_name

	lastP  = v.last_P_name
	lastQ  = v.last_Q_name
	lastA  = v.last_A_name
	lastc  = v.last_c_name
	lastu  = v.last_u_name

	run(`sh appender.sh $tempP $lastP`)
	run(`sh appender.sh $tempQ $lastQ`)
	run(`sh appender.sh $tempA $lastA`)
	run(`sh appender.sh $tempc $lastc`)
	run(`sh appender.sh $tempu $lastu`)

end

function transferTempToLast(vessels :: Array{Vessel, 1})

	for v in vessels
		transferTempToLast(v)
	end

end
