# Two-phase flow thermohydraulic oscillations
Continued work on [my master's thesis](https://zir.nsk.hr/islandora/object/fsb%3A9362).


# Governing Equations

## Mass Conservation Law
∂ρ/∂t + ∂(ρu)/∂x = 0

## Momentum Conservation Law
∂(ρu)/∂t + ∂(ρu²)/∂x + ∂p/∂x + (∂p/∂x)|_friction + ρg sin(β) = 0

## Energy Conservation Law
∂/∂t[ρ(h + u²/2)] + ∂/∂x[ρu(h + u²/2)] + ρu g = ∂p/∂t + q_w

## Heat Conduction in the Pipe Wall
ρc ∂θ_s/∂t = -∂q_x/∂x + φ_V

Rewriting in terms of heat flux:
∂θ_s/∂t = (1 / ρc) ∂/∂x(λ ∂θ_s/∂x) + φ_V / ρc

Simplified form:
∂θ_s/∂t = (λ / ρc) ∂²θ_s/∂x² + φ_V / ρc

