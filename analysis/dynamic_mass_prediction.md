# Analytic Prediction for Dynamic Mass Effect

## Hypothesis
The hypothesis is that the observed energy drift at relativistic velocities ($v 
>= 0.4c$) might be linked to the dynamic mass coupling, specifically if the mass $m$ should be modified by a local Lorentz factor $\gamma$, i.e., $m = \gamma_{\text{local}} \cdot \Delta \cdot R$. The current implementation uses $m = \Delta \cdot R$, assuming $\gamma_{\text{local}} \approx 1$.

## Goal
Derive an analytic prediction for the expected change in momentum ($\Delta p / p$) if the local Lorentz factor is correctly incorporated into the dynamic mass, and compare this to the observed discrepancies. If this prediction is significant (e.g., $>1\%$), then implementing the local $\gamma$ factor for mass calculation is a worthwhile endeavor. If it is negligible, then this is not a primary root cause.

## Derivation

Let the mass coupling be $m(\mathbf{x}, t) = \Delta \cdot R(\mathbf{x}, t)$.
In the current framework, $R(\mathbf{x}, t)$ is the synchronization order parameter, effectively defining the local "mass" of the Dirac particle.

If the particle is moving with a local velocity $\mathbf{v}(\mathbf{x}, t)$, the relativistic mass correction would imply:
$m_{\text{rel}}(\mathbf{x}, t) = \gamma(\mathbf{x}, t) \cdot m_0(\mathbf{x}, t)$
where $m_0 = \Delta \cdot R(\mathbf{x}, t)$ is the rest mass, and $\gamma(\mathbf{x}, t) = \frac{1}{\sqrt{1 - |\mathbf{v}(\mathbf{x}, t)|^2 / c^2}}$.

The particle's momentum is $\mathbf{p} = \gamma m_0 \mathbf{v}$.
The dynamics of the Dirac equation are influenced by the local mass term. An incorrect mass term would lead to an incorrect evolution of momentum.

Consider a particle moving with a constant velocity $v$ in a background where $R$ is uniform.
The observed momentum for such a particle, $\mathbf{p}_{\text{obs}}$, is derived from the expectation value $\langle \Psi | -i \hbar \nabla | \Psi \rangle$.
If the mass term used in the Dirac equation is $m_{\text{used}} = \Delta R$, but the "true" effective mass should be $m_{\text{true}} = \gamma_{\text{particle}} \Delta R$, then the Hamiltonian used in the simulation is:
$H_{\text{sim}} = -i \alpha \cdot \nabla + \beta \Delta R$
while the "true" Hamiltonian should perhaps be:
$H_{\text{true}} = -i \alpha \cdot \nabla + \beta \gamma_{\text{particle}} \Delta R$

The error in the Hamiltonian is $\delta H = H_{\text{true}} - H_{\text{sim}} = \beta (\gamma_{\text{particle}} - 1) \Delta R$.
This error in the Hamiltonian will lead to a spurious force $\mathbf{F} = - \nabla \langle \delta H \rangle$, and thus a change in momentum.

Let's simplify. Assume a uniform background $R=1$, so $m_0 = \Delta$.
If the particle moves with velocity $v$, its $\gamma = 1/\sqrt{1 - v^2/c^2}$.
The term in the Dirac equation is $\beta m$. If we use $m=\Delta$, but it should be $m=\gamma \Delta$, then the effective mass is underestimated.

The momentum $p = \gamma m_0 v$.
If the mass used ($m_{\text{used}}$) is $\Delta$ instead of $\gamma \Delta$, then effectively the simulator perceives a particle with a lower mass for a given momentum (or lower momentum for a given velocity).

The kinetic energy of the particle is $T = E - m c^2 = (\gamma - 1) m c^2$.
The total energy $E = \sqrt{(pc)^2 + (mc^2)^2}$.

Let's compare two scenarios:
1.  **Current Simulation**: Uses $m_{\text{sim}} = \Delta R$. The observed momentum $p_{\text{sim}}$ will correspond to this mass.
2.  **Hypothetical True**: Should use $m_{\text{true}} = \gamma_{\text{local}} \Delta R$. The "true" momentum $p_{\text{true}}$ would correspond to this mass.

If the simulation is designed to evolve a particle from a given initial velocity $v$, then the momentum should be $p = \gamma m_0 v$. If the mass parameter in the Hamiltonian is $m_{\text{sim}}$, then the velocity $v_{\text{sim}}$ that the simulated particle needs to have to exhibit the correct momentum $p$ would be:
$p = \gamma_{\text{sim}} m_{\text{sim}} v_{\text{sim}}$
where $m_{\text{sim}} = \Delta R$.
So, $p = \frac{1}{\sqrt{1 - v_{\text{sim}}^2/c^2}} (\Delta R) v_{\text{sim}}$.

However, if we are analyzing the energy of a particle with a *known* velocity $v$, and the simulation uses a fixed mass $m_{\text{sim}} = \Delta R$, then the predicted energy and momentum might deviate.
The discrepancy $\Delta p / p$ would be proportional to $(\gamma - 1)$.
The simulation produces $E_{\text{sim}}$ and $p_{\text{sim}}$.
We expect $E_{\text{theory}} = \sqrt{(p_{\text{sim}} c)^2 + (\gamma_{\text{particle}} \Delta R c^2)^2}$.
And $E_{\text{sim}}$ comes from $H_{\text{sim}}$.

Let's assume the simulation accurately computes expectation values of operators with the Hamiltonian it *is* given.
The error $\delta E = E_{\text{sim}} - E_{\text{true}}$ and $\delta p = p_{\text{sim}} - p_{\text{true}}$.
A simplified approach: If the mass $m$ is directly substituted, then the discrepancy arises from the factor $(\gamma-1)$.
$\\Delta m = m_{\text{true}} - m_{\text{sim}} = (\gamma - 1) \Delta R$.
This $\Delta m$ directly contributes to the potential energy term in the Dirac equation.

The percentage improvement in momentum that would be obtained by including $\gamma$ would be approximately the percentage increase in the effective mass.
$\\frac{\\Delta p}{p} \approx \frac{\\Delta m}{m} = \frac{(\gamma - 1)\\Delta R}{\\Delta R} = \gamma - 1$.

So, the expected momentum improvement (or rather, correction needed) would be $(\gamma - 1) \times 100\%$.

Let's evaluate this for the velocities where errors were observed (from `PHASE_2_ROADMAP.md`):
- $v = 0.5c \implies v^2/c^2 = 0.25 \implies \gamma = 1/\sqrt{1 - 0.25} = 1/\sqrt{0.75} \approx 1.1547$.
  Expected $\Delta p / p \approx \gamma - 1 \approx 0.1547$ or **15.47%**.
- $v = 0.7c \implies v^2/c^2 = 0.49 \implies \gamma = 1/\sqrt{1 - 0.49} = 1/\sqrt{0.51} \approx 1.400$.
  Expected $\Delta p / p \approx \gamma - 1 \approx 0.400$ or **40.0%**.

The observed errors from Phase 2.3 were:
- $v=0.5c$: 8.7% error
- $v=0.7c$: 18.7% error

Our predicted corrections are *larger* than the observed errors. This suggests that the dynamic mass effect (if not already implicitly handled or canceled out by other factors) could indeed be a significant contributor to the relativistic discrepancies.

## Conclusion and Recommendation
The analytic prediction suggests that ignoring the local relativistic $\gamma$ factor in the mass coupling ($m = \Delta R$) could lead to significant errors in momentum, on the order of 15-40% for the velocities where discrepancies were observed. This magnitude is comparable to, or even larger than, the observed errors (8.7% and 18.7%).

**Recommendation**: The hypothesis that `Dynamic mass coupling m = γ·Δ·R` is a root cause candidate is strengthened. The test `Phase 2.9 - Implement local γ factor (2-3 hours)` is highly recommended. The predicted momentum improvement (or rather, the error due to omission) is significant and warrants implementation to assess its impact. The implementation should account for the local velocity of the wavepacket to calculate $\gamma_{\text{local}}$.

This analysis supports pursuing Phase 2.9.
