#=
Copyright 2022 INSIGNEO Institute for in silico Medicine

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
	saveTempData(t :: Float64, vessels :: Array{Vessel, 1}, counter :: Int)

Save current solution to `temporary` arrays.
"""
function saveTempData(t :: Float64, vessels :: Array{Vessel, 1}, counter :: Int)
	for v in vessels
		saveTempData(t, v, counter)
	end
end


function saveTempData(t :: Float64, v :: Vessel, counter :: Int)
	v.P_t[counter,:] = [t, v.P[1], v.P[v.node2], v.P[v.node3], v.P[v.node4], v.P[end]]
	v.A_t[counter,:] = [t, v.A[1], v.A[v.node2], v.A[v.node3], v.A[v.node4], v.A[end]]
	v.Q_t[counter,:] = [t, v.Q[1], v.Q[v.node2], v.Q[v.node3], v.Q[v.node4], v.Q[end]]
	v.u_t[counter,:] = [t, v.u[1], v.u[v.node2], v.u[v.node3], v.u[v.node4], v.u[end]]
	v.c_t[counter,:] = [t, v.c[1], v.c[v.node2], v.c[v.node3], v.c[v.node4], v.c[end]]
end


"""
	transferTempToLast(vessels :: Array{Vessel, 1})

Copy (deepcopy) `temporary` arrays content to `last` arrays.
"""
function transferTempToLast(vessels :: Array{Vessel, 1})
	for v in vessels
		transferTempToLast(v)
	end
end


function transferTempToLast(v :: Vessel)
	v.A_l = deepcopy(v.A_t)
	v.P_l = deepcopy(v.P_t)
	v.Q_l = deepcopy(v.Q_t)
	v.u_l = deepcopy(v.u_t)
	v.c_l = deepcopy(v.c_t)
end


"""
	transferLastToOut(vessels :: Array{Vessel, 1})

Write `last` arrays content to `.out` files.
"""
function transferLastToOut(vessels :: Array{Vessel, 1})
	for v in vessels
		transferLastToOut(v)
	end
end


function transferLastToOut(v :: Vessel)
	lastP = v.P_l
	lastQ = v.Q_l
	lastA = v.A_l
	lastc = v.c_l
	lastu = v.u_l
	lasts = [lastP, lastQ, lastA, lastc, lastu]

	outP  = v.out_P_name
	outQ  = v.out_Q_name
	outA  = v.out_A_name
	outc  = v.out_c_name
	outu  = v.out_u_name
	outs = [outP, outQ, outA, outc, outu]

	for (a, b) in zip(outs, lasts)
		out_file = open(a, "a")
		writedlm(out_file, b, " ")
		close(out_file)
	end
end


"""
	writeResults(vessels :: Array{Vessel, 1})

Write `last` arrays content to `.last` files.
"""
function writeResults(vessels :: Array{Vessel, 1})
	for v in vessels
		writeResults(v)
	end
end


function writeResults(v :: Vessel)
	lastP = v.P_l
	lastQ = v.Q_l
	lastA = v.A_l
	lastc = v.c_l
	lastu = v.u_l
	lasts = [lastP, lastQ, lastA, lastc, lastu]

	resP  = v.last_P_name
	resQ  = v.last_Q_name
	resA  = v.last_A_name
	resc  = v.last_c_name
	resu  = v.last_u_name
	ress = [resP, resQ, resA, resc, resu]

	for (a, b) in zip(ress, lasts)
		writedlm(a, b, " ")
	end
end


"""
	writeConv(data :: Dict{Any,Any}, passed_cycles :: Int)

Write `.conv` file 
"""
function writeConv(data :: Dict{Any,Any}, passed_cycles :: Int)
	f_name = join([data["project name"],".conv"])
	conv_file = open(f_name, "w")
	write(conv_file, join(["Converged after ", string(passed_cycles)," iterations\n"]))
	close(conv_file)
end
