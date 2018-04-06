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


"""
	openTempFiles(vessels :: Array{Vessel, 1})
"""
function openTempFiles(vessels :: Array{Vessel, 1})
	for v in vessels
		openTempFiles(v)
	end
end


function openTempFiles(v :: Vessel)
	v.temp_P = open(v.temp_P_name, "w")
	v.temp_Q = open(v.temp_Q_name, "w")
	v.temp_A = open(v.temp_A_name, "w")
	v.temp_c = open(v.temp_c_name, "w")
	v.temp_u = open(v.temp_u_name, "w")
end


"""
	closeTempFiles(vessels :: Array{Vessel, 1})
"""
function closeTempFiles(vessels :: Array{Vessel, 1})
	for v in vessels
		closeTempFiles(v)
	end
end


function closeTempFiles(v :: Vessel)
	close(v.temp_P)
	close(v.temp_Q)
	close(v.temp_A)
	close(v.temp_c)
	close(v.temp_u)
end


"""
	openCloseLastFiles(vessels :: Array{Vessel, 1})

initialise empty `.last` files.
"""
function openCloseLastFiles(vessels :: Array{Vessel, 1})
	for v in vessels
		openCloseLastFiles(v)
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


"""
	saveTempData(t :: Float64, vessels :: Array{Vessel, 1})

Write the current time solution in `.temp` files.
"""
function saveTempData(t :: Float64, vessels :: Array{Vessel, 1})
	for v in vessels
		saveTempData(t, v)
	end
end


function saveTempData(t :: Float64, v :: Vessel)
	println(v.temp_P, t, " ", v.P[1], " ", v.P[v.node2], " ", v.P[v.node3], " ", v.P[v.node4], " ", v.P[end])
	println(v.temp_A, t, " ", v.A[1], " ", v.A[v.node2], " ", v.A[v.node3], " ", v.A[v.node4], " ", v.A[end])
	println(v.temp_Q, t, " ", v.Q[1], " ", v.Q[v.node2], " ", v.Q[v.node3], " ", v.Q[v.node4], " ", v.Q[end])
	println(v.temp_c, t, " ", v.c[1], " ", v.c[v.node2], " ", v.c[v.node3], " ", v.c[v.node4], " ", v.c[end])
	println(v.temp_u, t, " ", v.u[1], " ", v.u[v.node2], " ", v.u[v.node3], " ", v.u[v.node4], " ", v.u[end])
end


"""
	transferTempToOut(vessels :: Array{Vessel, 1})

Append `.temp` content to `.out` file.
"""
function transferTempToOut(vessels :: Array{Vessel, 1})
	for v in vessels
		transferTempToOut(v)
	end
end


function transferTempToOut(v :: Vessel)
	tempP = v.temp_P_name
	tempQ = v.temp_Q_name
	tempA = v.temp_A_name
	tempc = v.temp_c_name
	tempu = v.temp_u_name
	temps = [tempP, tempQ, tempA, tempc, tempu]

	outP  = v.out_P_name
	outQ  = v.out_Q_name
	outA  = v.out_A_name
	outc  = v.out_c_name
	outu  = v.out_u_name
	outs = [outP, outQ, outA, outc, outu]

	for (a, b) in zip(outs, temps)
		out_file = open(a, "a")
		temp_file = open(b, "r")
		write(out_file, temp_file)
		close(out_file)
		close(temp_file)
	end
end


"""
	transferLastToOut(vessels :: Array{Vessel, 1})

Move `.last` content to `.out` files.
"""
function transferLastToOut(vessels :: Array{Vessel, 1})
	for v in vessels
		transferLastToOut(v)
	end
end


function transferLastToOut(v :: Vessel)
	lastP = v.last_P_name
	lastQ = v.last_Q_name
	lastA = v.last_A_name
	lastc = v.last_c_name
	lastu = v.last_u_name
	lasts = [lastP, lastQ, lastA, lastc, lastu]

	outP  = v.out_P_name
	outQ  = v.out_Q_name
	outA  = v.out_A_name
	outc  = v.out_c_name
	outu  = v.out_u_name
	outs = [outP, outQ, outA, outc, outu]

	for (a, b) in zip(outs, lasts)
		out_file = open(a, "a")
		last_file = open(b, "r")
		write(out_file, last_file)
		close(out_file)
		close(last_file)
	end
end


"""
	transferTempToLast(vessels :: Array{Vessel, 1})

Move `.temp` content to `.last` files.
"""
function transferTempToLast(vessels :: Array{Vessel, 1})
	for v in vessels
		transferTempToLast(v)
	end
end


function transferTempToLast(v :: Vessel)
	tempP = v.temp_P_name
	tempQ = v.temp_Q_name
	tempA = v.temp_A_name
	tempc = v.temp_c_name
	tempu = v.temp_u_name
	temps = [tempP, tempQ, tempA, tempc, tempu]

	lastP  = v.last_P_name
	lastQ  = v.last_Q_name
	lastA  = v.last_A_name
	lastc  = v.last_c_name
	lastu  = v.last_u_name
	lasts = [lastP, lastQ, lastA, lastc, lastu]

	for (a, b) in zip(temps, lasts)
		temp_file = open(a, "r")
		last_file = open(b, "w")
		write(last_file, temp_file)
		close(temp_file)
		close(last_file)
	end
end


"""
	cleanTemps(vessels :: Array{Vessel, 1})

Delete `.temp` files from the results folder.
"""
function cleanTemps(vessels :: Array{Vessel, 1})
	for v in vessels
		cleanTemps(v)
	end
end


function cleanTemps(v :: Vessel)
	tempP = v.temp_P_name
	tempQ = v.temp_Q_name
	tempA = v.temp_A_name
	tempc = v.temp_c_name
	tempu = v.temp_u_name
	temps = [tempP, tempQ, tempA, tempc, tempu]

	for temp in temps
		rm(temp)
	end
end


"""
	cleanOuts(vessels :: Array{Vessel, 1})

Delete `.out` files from the results folder.
"""
function cleanOuts(vessels :: Array{Vessel, 1})
	for v in vessels
		cleanOuts(v)
	end
end


function cleanOuts(v :: Vessel)
	outP = v.out_P_name
	outQ = v.out_Q_name
	outA = v.out_A_name
	outc = v.out_c_name
	outu = v.out_u_name
	outs = [outP, outQ, outA, outc, outu]

	for out in outs
		rm(out)
	end
end
