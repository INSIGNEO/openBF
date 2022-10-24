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
    solveConjunction(b :: Blood, v1 :: Vessel, v2 :: Vessel)

Solve the non-linear system at the conjunction node between two vessels.
"""
function solveConjunction(b :: Blood, v1 :: Vessel, v2 :: Vessel)
    U0 = @SArray [v1.u[end],
         v2.u[1],
         sqrt(sqrt(v1.A[end])),
         sqrt(sqrt(v2.A[1]))]

    k1 = v1.s_15_gamma[end]
    k2 = v2.s_15_gamma[1]
    k3 = b.rho
    k  = @SArray [k1, k2, k3]

    J = calculateJacobianConjunction(v1, v2, U0, k)
    U = newtonRaphson([v1, v2], J, U0, k, calculateWstarConjunction, calculateFConjunction)

    updateConjunction(U, v1, v2)
end


"""
    calculateJacobianConjunction(v1 :: Vessel, v2 :: Vessel, U, k)

Return the Jacobian for conjunction equations.
"""
function calculateJacobianConjunction(v1 :: Vessel, v2 :: Vessel, U, k)
    U33 = U[3]*U[3]*U[3]
    U43 = U[4]*U[4]*U[4]

    J13 =  4.0*k[1]
    J24 = -4.0*k[2]

    J31 =  U33*U[3]
    J32 = -U43*U[4]
    J33 =  4.0*U[1]*U33
    J34 = -4.0*U[2]*U43

    J41 =  k[3]*U[1]
    J42 = -k[3]*U[2]
    J43 =  2.0*v1.beta[end]*U[3]*v1.s_inv_A0[end]
    J44 = -2.0*v2.beta[1]*U[4]*v2.s_inv_A0[1]

    return @SArray [1.0 0.0 J13 0.0;
                    0.0 1.0 0.0 J24;
                    J31 J32 J33 J34;
                    J41 J42 J43 J44]
end


"""
    calculateWstarConjunction(U, k)

Return the Riemann invariants at the conjunction node.
"""
function calculateWstarConjunction(U, k)
    W1 = U[1] + 4.0*k[1]*U[3]
    W2 = U[2] - 4.0*k[2]*U[4]

    return @SArray [W1, W2]
end


"""
    calculateFConjunction(vessels :: Array{Vessel,1}, U, k, W)
"""
function calculateFConjunction(vessels :: Array{Vessel,1}, U, k, W)
    v1 = vessels[1]
    v2 = vessels[2]

    U32 = U[3]*U[3]
    U42 = U[4]*U[4]

    f1 = U[1] + 4.0*k[1]*U[3] - W[1]

    f2 = U[2] - 4.0*k[2]*U[4] - W[2]

    f3 = U[1]*U32*U32 - U[2]*U42*U42

    f4 = 0.5*k[3]*U[1]*U[1] + v1.beta[end]*(U32*v1.s_inv_A0[end] - 1.0) -
        (0.5*k[3]*U[2]*U[2] + v2.beta[1]*(U42*v2.s_inv_A0[1] - 1.0))

    return @SArray [f1, f2, f3, f4]
end


"""
    updateConjunction(U, v1 :: Vessel, v2 :: Vessel)

Update the values at the conjunction node for the two vessels.
"""
function updateConjunction(U, v1 :: Vessel, v2 :: Vessel)
    v1.u[end] = U[1]
    v2.u[1] = U[2]

    v1.A[end] = U[3]*U[3]*U[3]*U[3]
    v1.Q[end] = v1.u[end]*v1.A[end]

    v2.A[1] = U[4]*U[4]*U[4]*U[4]
    v2.Q[1] = v2.u[1]*v2.A[1]

    v1.P[end] = pressure(v1.A[end], v1.A0[end], v1.beta[end], v1.Pext)
    v2.P[1] = pressure(v2.A[1], v2.A0[1], v2.beta[1], v2.Pext)

    v1.c[end] = waveSpeed(v1.A[end], v1.gamma[end])
    v2.c[1] = waveSpeed(v2.A[1], v2.gamma[1])
end
