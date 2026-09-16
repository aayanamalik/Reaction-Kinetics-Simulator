# Reaction Kinetics Simulator

Python models of chemical reaction kinetics, from analytically solvable first-order
decay through to a non-isothermal reactor whose ignition transient defeats fixed-step
explicit integration.

![Haber process ignition](docs/haber.png)

## Models

| # | Model | Method |
|---|-------|--------|
| 1 | First-order decay, A → B | Analytic |
| 2 | Rate constant comparison | Analytic |
| 3 | Temperature dependence (Arrhenius) | Analytic |
| 4 | Consecutive reactions, A → B → C | Analytic |
| 5 | Reversible equilibrium, A ⇌ B | Analytic |
| 6 | Haber process, non-isothermal | Heun (explicit, fixed step) |
| 7 | Haber process, non-isothermal | SciPy Radau (implicit, adaptive) |

Models 6 and 7 integrate the same system from the same parameters, so the two
solvers can be compared directly.

## The system

Four coupled ODEs in `y = [N₂, H₂, NH₃, T]`. A mass-action rate law
(forward `k_f[N₂][H₂]³`, reverse `k_r[NH₃]²`) is coupled to an energy balance over a
jacketed reactor, with both rate constants set by the Arrhenius equation. The
coupling runs in a loop: heat release raises the temperature, temperature raises the
rate constants exponentially, and higher rate constants raise heat release further.

### Default parameters

```python
# Chemistry
A0, B0, C0 = 1.0, 3.0, 0.0   # initial N2, H2, NH3 (mol dm^-3)
T0 = 450                     # initial temperature (K)
forward_Ea = 110000          # J/mol, iron catalyst
reverse_Ea = 180000          # J/mol, iron catalyst
A_factor   = 5e8             # pre-exponential factor
reactant_orders = [1, 3, 0]  # [N2]^1 [H2]^3
product_orders  = [0, 0, 2]  # [NH3]^2
coeffs = [-1, -3, 2]         # stoichiometry per unit of reaction

# Reactor and thermodynamics (200 atm)
delta_H = -92400             # J per mol of N2
Cp      = 2000               # J/kg/K, constant
density = 60.0               # kg/m3, constant
volume  = 1                  # dm3
surroundings_temp = 750      # K, jacket
UA = 5000.0                  # W/K
end_time = 100               # s
```

These values are chosen to make the system stiff, not fitted to plant data. See
[Limitations](#limitations).

## Why two solvers

Under the defaults above the reactor ignites at t ≈ 0.024 s, spiking to 1279 K before
the jacket pulls it back to 750 K. The whole event lasts a few hundredths of a second
inside a 0.5 second integration, and that mismatch is what breaks fixed-step explicit
integration.

Measured against that prediction:

| Solver | Steps | Function evals | Result |
|--------|------:|---------------:|--------|
| Heun, dt = 5 × 10⁻³ | — | — | diverges at t = 0.025 s |
| Heun, dt = 5 × 10⁻⁴ | — | — | diverges at t = 0.024 s |
| Heun, dt = 2 × 10⁻⁴ | 500,000 | 1,000,000 | completes |
| RK45 (explicit, adaptive) | — | 9,374 | completes |
| Radau (implicit, adaptive) | 253 | 1,979 | completes |

Radau needs about 250× fewer steps than the smallest working fixed step. The comparison with RK45 separates the
two effects: adaptivity alone accounts for most of the saving, and implicitness saves a further factor of 4.7. The Heun's implementation used has a fixed timestep and therefore requires a huge number of calculations to complete successfully.

## Results

**Numerics.** See the table above. Stiffness here is driven by the ignition
transient, not by a wide separation of steady timescales.

**Physics.** The image at the top uses a milder parameter set (`A_factor = 5e4`,
`UA = 2.0`, `surroundings_temp = 650`, `end_time = 400`) chosen to show thermal
runaway on a visible timescale. Holding everything else fixed and varying only the
cooling:

| UA (W/K) | Peak T (K) | Excess over jacket (K) |
|---------:|-----------:|-----------------------:|
| 10 | 671.5 | 21.5 |
| 2 | 1146.9 | 496.9 |

The same reactor, one parameter apart: adequate cooling gives smooth operation
peaking 21 K above the jacket, while reducing UA fivefold produces runaway to 1147 K
after 181 seconds. This is parametric sensitivity — the classic reason jacketed
exothermic reactors are designed with margin rather than to nominal duty.

## Running

Requires [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/aayanamalik/Reaction-Kinetics-Simulator
cd Reaction-Kinetics-Simulator
uv run src/reactions/simulator.py
```

## Limitations

- Cp and mixture density are treated as constants, so there is no upper bound on
  temperature from dissociation, phase change, or variable heat capacity. Runaway
  cases are therefore only qualitatively meaningful.
- Forward and reverse rate constants share a pre-exponential factor. The implied
  equilibrium constant is not thermodynamically consistent with ΔH and ΔS.
- The default `UA = 5000 W/K` implies far more heat transfer area than a 1 dm³
  reactor could have, and `A_factor = 5e8` is not fitted to kinetic data. Both are
  chosen to produce stiff behaviour for the solver comparison.
- Conditions are assumed isobaric and the reactor perfectly mixed.
