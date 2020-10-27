#=
Copyright 2020 INSIGNEO Institute for in silico Medicine

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
    solveAnastomosis(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel)

Solve the non-linear system at the anastomosis node between three vessels.
"""
function solveAnastomosis(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel)
    @fastmath @inbounds U0 = @SArray [v1.u[end],
         v2.u[end],
         v3.u[1],
         sqrt(sqrt(v1.A[end])),
         sqrt(sqrt(v2.A[end])),
         sqrt(sqrt(v3.A[1]))]

    @inbounds k1 = v1.s_15_gamma[end]
    @inbounds k2 = v2.s_15_gamma[end]
    @inbounds k3 = v3.s_15_gamma[1]
    k = @SArray [k1, k2, k3]

    J = calculateJacobianAnastomosis(v1, v2, v3, U0, k)
    U = newtonRaphson([v1, v2, v3], J, U0, k, calculateWstarAnastomosis, calculateFAnastomosis)

    updateAnastomosis(U, v1, v2, v3)
end


"""
    calculateJacobianAnastomosis(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel, U, k)

Return the Jacobian for anastomosis equations.
"""
function calculateJacobianAnastomosis(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel, U, k)
    @fastmath @inbounds g1 = sqrt(1-v1.Ac/U[4]^4)	# MODIFIED THIS LINE
    @fastmath @inbounds g2 = sqrt(1-v2.Ac/U[5]^4)	# MODIFIED THIS LINE
    @fastmath @inbounds g3 = sqrt(1-v3.Ac/U[6]^4)	# MODIFIED THIS LINE
    @fastmath @inbounds U43 = U[4]*U[4]*U[4]
    @fastmath @inbounds U53 = U[5]*U[5]*U[5]
    @fastmath @inbounds U63 = U[6]*U[6]*U[6]

    @fastmath @inbounds J14 =  4.0*k[1]*(v1.Ac/(2*U[4]*g1) + g1)	# MODIFIED THIS LINE
    @fastmath @inbounds J25 =  4.0*k[2]*(v2.Ac/(2*U[5]*g2) + g2) 	# MODIFIED THIS LINE
    @fastmath @inbounds J36 = -4.0*k[3]*(v3.Ac/(2*U[6]*g3) + g3)	# MODIFIED THIS LINE

    @fastmath @inbounds J41 =  U[4]*U43 - v1.Ac				# MODIFIED THIS LINE
    @fastmath @inbounds J42 =  U[5]*U53 - v2.Ac				# MODIFIED THIS LINE
    @fastmath @inbounds J43 = -U[6]*U63 + v3.Ac				# MODIFIED THIS LINE
    @fastmath @inbounds J44 =  4.0*U[1]*U43
    @fastmath @inbounds J45 =  4.0*U[2]*U53
    @fastmath @inbounds J46 = -4.0*U[3]*U63

    @fastmath @inbounds J54 =  2.0*v1.beta[end]*U[4]*v1.s_inv_A0[end]
    @fastmath @inbounds J56 = -2.0*v3.beta[1]*U[6]*v3.s_inv_A0[1]

    @fastmath @inbounds J65 =  2.0*v2.beta[end]*U[5]*v2.s_inv_A0[end]
    @fastmath @inbounds J66 = -2.0*v3.beta[1]*U[6]*v3.s_inv_A0[1]

    return @SArray [1.0 0.0 0.0 J14 0.0 0.0;
                    0.0 1.0 0.0 0.0 J25 0.0;
                    0.0 0.0 1.0 0.0 0.0 J36;
                    J41 J42 J43 J44 J45 J46;
                    0.0 0.0 0.0 J54 0.0 J56;
                    0.0 0.0 0.0 0.0 J65 J66]
end


"""
    calculateWstarAnastomosis(U, k)

Return the Riemann invariants at the anastomosis node.
"""
function calculateWstarAnastomosis(U, k, vessels :: Array{Vessel,1})
	# MODIFIED THIS FUNCTION
    v1 = vessels[1]						# MODIFIED THIS LINE
    v2 = vessels[2]						# MODIFIED THIS LINE
    v3 = vessels[3]						# MODIFIED THIS LINE
    @fastmath @inbounds g1 = sqrt(1-v1.Ac/U[4]^4)		# MODIFIED THIS LINE
    @fastmath @inbounds g2 = sqrt(1-v2.Ac/U[5]^4)		# MODIFIED THIS LINE
    @fastmath @inbounds g3 = sqrt(1-v3.Ac/U[6]^4)		# MODIFIED THIS LINE
    @fastmath @inbounds W1 = U[1] + 4*k[1]*U[4]*g1		# MODIFIED THIS LINE
    @fastmath @inbounds W2 = U[2] + 4*k[2]*U[5]*g2		# MODIFIED THIS LINE
    @fastmath @inbounds W3 = U[3] - 4*k[3]*U[6]*g3		# MODIFIED THIS LINE

    return @SArray [W1, W2, W3]
end


"""
    calculateFAnastomosis(vessels :: Array{Vessel,1}, U, k, W)
"""
function calculateFAnastomosis(vessels :: Array{Vessel,1}, U, k, W)
    v1 = vessels[1]
    v2 = vessels[2]
    v3 = vessels[3]

    @fastmath @inbounds g1 = sqrt(1-v1.Ac/U[4]^4)		# MODIFIED THIS LINE
    @fastmath @inbounds g2 = sqrt(1-v2.Ac/U[5]^4)		# MODIFIED THIS LINE
    @fastmath @inbounds g3 = sqrt(1-v3.Ac/U[6]^4)		# MODIFIED THIS LINE

    @fastmath @inbounds U42 = U[4]*U[4]
    @fastmath @inbounds U52 = U[5]*U[5]
    @fastmath @inbounds U62 = U[6]*U[6]

    @fastmath @inbounds f1 = U[1] + 4*k[1]*U[4]*g1 - W[1]	# MODIFIED THIS LINE

    @fastmath @inbounds f2 = U[2] + 4*k[2]*U[5]*g2 - W[2]	# MODIFIED THIS LINE

    @fastmath @inbounds f3 = U[3] - 4*k[3]*U[6]*g3 - W[3]	# MODIFIED THIS LINE

    @fastmath @inbounds f4 = U[1]*(U42*U42 - v1.Ac) + U[2]*(U52*U52 - v2.Ac) - U[3]*(U62*U62 - v3.Ac)	# MODIFIED THIS LINE

    f5 = v1.beta[end]*(U42*v1.s_inv_A0[end] - 1.0) -
        (v3.beta[1]*(U62*v3.s_inv_A0[1] - 1.0))

    f6 = v2.beta[1]*(U52*v2.s_inv_A0[end] - 1.0) -
        (v3.beta[1]*(U62*v3.s_inv_A0[1] - 1.0))

    return @SArray [f1, f2, f3, f4, f5, f6]
end


"""
    updateAnastomosis(U, v1 :: Vessel, v2 :: Vessel, v3 :: Vessel)

Update the values at the anastomosis node for the three vessels.
"""
function updateAnastomosis(U, v1 :: Vessel, v2 :: Vessel, v3 :: Vessel)
    @inbounds v1.u[end] = U[1]
    @inbounds v2.u[end] = U[2]
    @inbounds v3.u[1] = U[3]

    @fastmath @inbounds v1.A[end] = U[4]*U[4]*U[4]*U[4]
    @fastmath @inbounds v2.A[end] = U[5]*U[5]*U[5]*U[5]
    @fastmath @inbounds v3.A[1] = U[6]*U[6]*U[6]*U[6]

    @fastmath @inbounds v1.Q[end] = v1.u[end]*(v1.A[end] - v1.Ac)	# MODIFIED THIS LINE
    @fastmath @inbounds v2.Q[end] = v2.u[end]*(v2.A[end] - v2.Ac)	# MODIFIED THIS LINE
    @fastmath @inbounds v3.Q[1] = v3.u[1]*(v3.A[1] - v3.Ac)		# MODIFIED THIS LINE

    @inbounds v1.P[end] = pressure(v1.A[end], v1.A0[end], v1.beta[end], v1.Pext)
    @inbounds v2.P[end] = pressure(v2.A[end], v2.A0[end], v2.beta[end], v2.Pext)
    @inbounds v3.P[1] = pressure(v3.A[1], v3.A0[1], v3.beta[1], v3.Pext)

    @fastmath @inbounds v1.c[end] = waveSpeed(v1.A[end], v1.gamma[end], v1.Ac)	# MODIFIED THIS LINE
    @fastmath @inbounds v2.c[end] = waveSpeed(v2.A[end], v2.gamma[end], v2.Ac)	# MODIFIED THIS LINE
    @fastmath @inbounds v3.c[1] = waveSpeed(v3.A[1], v3.gamma[1], v3.Ac)	# MODIFIED THIS LINE
end
