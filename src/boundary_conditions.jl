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
    setInletBC(t :: Float64, dt :: Float64, v :: Vessel,  h  :: Heart)

Flow or pressure value is set from the inlet file stored in the `Heart` structure.
"""
function setInletBC(t :: Float64, dt :: Float64, v :: Vessel)
    h = v.heart

  	if h.inlet_type == "Q"
		v.Q[1] = inputFromData(t, h)
  	else
		v.P[1] = inputFromData(t, h)
  	end

	inletCompatibility(dt, v, h)
end


"""
    inputFromData(t :: Float64, h :: Heart)

Inlet flow waveform is stored in `h` data structure inside `input_data` matrix.
`input_data` has two columns and as many rows as time steps.
"""
function inputFromData(t :: Float64, h :: Heart)
	idt = h.input_data[:,1]
	idq = h.input_data[:,2]

	t_hat = div(t,h.cardiac_T)

	t -= t_hat*h.cardiac_T

	idx = 0
	@inbounds for i in 1:length(idt)
		if ((t >= idt[i]) && (t <= idt[i+1]))
			idx=i
			break
		end
	end

	@inbounds qu = idq[idx] + (t - idt[idx]) * (idq[idx+1] - idq[idx]) /
       (idt[idx+1] - idt[idx])

	return qu
end


"""
    inletCompatibility(dt :: Float64, v :: Vessel, h :: Heart)

Compute compatibility relations at inlet node.
"""
function inletCompatibility(dt :: Float64, v :: Vessel, h :: Heart)
	W11, W21 = riemannInvariants(1, v)
	W12, W22 = riemannInvariants(2, v)

	@fastmath @inbounds W11 += (W12 - W11)*(v.c[1] - v.u[1])*dt*v.invDx
	@fastmath @inbounds W21 = 2.0*v.Q[1]/v.A[1] - W11

	v.u[1], v.c[1] = inverseRiemannInvariants(W11, W21)

	if h.inlet_type == "Q"
		@fastmath @inbounds v.A[1] = v.Q[1]/v.u[1]
		v.P[1] = pressure(v.A[1], v.A0[1], v.beta[1], v.Pext)
	else
		v.A[1] = areaFromPressure(v.P[1], v.A0[1], v.beta[1], v.Pext)
		@fastmath @inbounds v.Q[1] = v.u[1]*v.A[1]
	end

end


"""
    riemannInvariants(i :: Int, v :: Vessel)

Calculate Riemann invariants at the node `i` from `u` and `c`.
"""
function riemannInvariants(i :: Int, v :: Vessel)
  @fastmath @inbounds W1 = v.u[i] - 4.0*v.c[i]
  @fastmath @inbounds W2 = v.u[i] + 4.0*v.c[i]

  return W1, W2
end


"""
    inverseRiemannInvariants(W1 :: Float64, W2 :: Float64)

Calculate `u` and `c` given `W1` and `W2`
"""
function inverseRiemannInvariants(W1 :: Float64, W2 :: Float64)
  @fastmath u = 0.5*(W1 + W2)
  @fastmath c = (W2 - W1)*0.125

  return u, c
end


"""
    areaFromPressure(P :: Float64, A0 :: Float64, beta :: Float64, Pext :: Float64)

Inverse constitutive equation. This is used only when a pressure inlet-time-function is
imposed (not recommended).
"""
function areaFromPressure(P :: Float64, A0 :: Float64, beta :: Float64, Pext :: Float64)
   return A0 * ((P-Pext)/beta + 1.0)*((P-Pext)/beta + 1.0)
end


"""
    setOutletBC(dt :: Float64, v :: Vessel)

Outlet boundary condition is applied to the last node of the involved vessel. Two
boundary conditions are available. Either a reflection coefficient is provided or a
three-element windkessel is coupled.
"""
function setOutletBC(dt :: Float64, v :: Vessel)
    if v.outlet == "reflection"
        v.P[end] = 2.0*v.P[end-1] - v.P[end-2]
		outletCompatibility(dt, v)
    elseif v.outlet == "wk3"
        threeElementWindkessel(dt, v)
	end
end


"""
    function outletCompatibility(dt :: Float64, v :: Vessel)

Outlet compatibility relations compute all the quantities not directly assigned by
the outlet boundary condition.
"""
function outletCompatibility(dt :: Float64, v :: Vessel)
	W1M1, W2M1 = riemannInvariants(v.M-1, v)
	W1M, W2M   = riemannInvariants(v.M, v)

	W2M += (W2M1 - W2M)*(v.u[end] + v.c[end])*dt/v.dx
	W1M = v.W1M0 - v.Rt * (W2M - v.W2M0)

	v.u[end], v.c[end] = inverseRiemannInvariants(W1M, W2M)
	v.Q[end] = v.A[end]*v.u[end]
end


"""
    function threeElementWindkessel(dt :: Float64, v :: Vessel)

The three element windkessel simulates the perfusion of downstream vessels.
This 0D model is coupled by `wk3` function to 1D model terminal branches via the
solution of a Riemann problem at the 0D/1D interface.
"""
function threeElementWindkessel(dt :: Float64, v :: Vessel)
	Pout = 0.0

	Al = v.A[end]
	ul = v.u[end]

	v.Pc += dt/v.Cc * (Al*ul - (v.Pc - Pout)/v.R2)

	As = Al

	ssAl = sqrt(sqrt(Al))
	sgamma =  2*sqrt(6*v.gamma[end])
	sA0 = sqrt(v.A0[end])
	bA0 = v.beta[end]/sA0

	fun(As) = As*v.R1*(ul + sgamma * (ssAl - sqrt(sqrt(As)))) -
		(v.Pext + bA0*(sqrt(As) - sA0)) + v.Pc

	dfun(As) = v.R1*(ul + sgamma * (ssAl - 1.25*sqrt(sqrt(As)))) - bA0*0.5/sqrt(As)

    try
        As = newtonSolver(fun, dfun, As)
    catch e
        vlab = v.label
        println("\nNewton solver doesn't converge at $vlab outlet!")
        throw(e)
    end

	us = (pressure(As, v.A0[end], v.beta[end], v.Pext) - Pout)/(As*v.R1)

	v.A[end] = As
	v.u[end] = us
end


"""
    function newtonSolver(f, df, x0 :: Float64)

Solve windkessel equation by means of Newton method.
"""
function newtonSolver(f, df, x0 :: Float64)
	xn = x0 - f(x0)/df(x0)
	if abs(xn-x0)<= 1e-5
		return xn
	else
		newtonSolver(f, df, xn)
	end
end


"""
    updateGhostCells(v :: Vessel)

Ghost cells are updated by assuming both ends as transmissive boundaries.
"""
function updateGhostCells(v :: Vessel)
	v.U00A = v.A[1]
	v.U00Q = v.Q[1]
	v.U01A = v.A[2]
	v.U01Q = v.Q[2]

	v.UM1A = v.A[v.M]
	v.UM1Q = v.Q[v.M]
	v.UM2A = v.A[v.M-1]
	v.UM2Q = v.Q[v.M-1]
end


function updateGhostCells(vessels :: Array{Vessel,1})
	for vessel in vessels
		updateGhostCells(vessel)
	end
end
