from scipy.integrate import solve_ivp, odeint
import numpy as np
from math import exp


def solve_PRmodel(t_dur, g_c, I_stim, stim_start, stim_end, stim_loc='soma'):
    """
    Solve the Pinsky-Rinzel model for a two-compartment neuron

    Parameters:
    -----------
    t_dur : float
        Duration of the simulation in ms
    g_c : float
        Coupling conductance between soma and dendrite in mS/cm^2
    I_stim : float
        Stimulus current amplitude in μA/cm^2
    stim_start : float
        Time when stimulus starts in ms
    stim_end : float
        Time when stimulus ends in ms
    stim_loc : str, optional
        Location where stimulus is applied, either 'soma' or 'dendrite'
        Default is 'soma'

    Returns:
    --------
    sol : OdeResult
        Solution from solve_ivp containing time points and state variables
    """
    # Check that stimulus location is valid
    if stim_loc not in ['soma', 'dendrite']:
        raise ValueError("stim_loc must be either 'soma' or 'dendrite'")
    # Parameters from the paper
    C_m = 2.0     # membrane capacitance [uF cm**-2]
    p = 0.5       # proportion of the membrane area taken up by the soma

    # Conductances
    g_C = 0       # Calcium-dependent K+ current [mS cm**-2]
    g_L = 2.0     # leak conductance [mS cm**-2]
    g_Na = 20.0   # sodium conductance [mS cm**-2]
    g_DR = 20     # Delayed rectifier K+ current [mS cm**-2]
    g_AHP = 0     # After-hyperpolarization K+ current [mS cm**-2]
    g_Ca = 10     # [mS cm**-2]

    E_L = -70.0   # leak reversal potential [mV]
    E_Na = 50.0   # sodium reversal potential [mV]
    E_K = -100.0  # potassium reversal potential [mV]
    E_Ca = 120.0  # calcium reversal potential [mV]

    def alpha_m(Vs):
        V1 = Vs + 46.9
        alpha = - 0.32 * V1 / (exp(-V1 / 4.) - 1.)
        return alpha

    def beta_m(Vs):
        V2 = Vs + 19.9
        beta = 0.28 * V2 / (exp(V2 / 5.) - 1.)
        return beta

    def alpha_h(Vs):
        alpha = 0.128 * exp((-43. - Vs) / 18.)
        return alpha

    def beta_h(Vs):
        V5 = Vs + 20.
        beta = 4. / (1 + exp(-V5 / 5.))
        return beta

    def alpha_n(Vs):
        V3 = Vs + 24.9
        alpha = - 0.016 * V3 / (exp(-V3 / 5.) - 1)
        return alpha

    def beta_n(Vs):
        V4 = Vs + 40.
        beta = 0.25 * exp(-V4 / 40.)
        return beta

    def alpha_s(Vd):
        alpha = 1.6 / (1 + exp(-0.072 * (Vd-5.)))
        return alpha

    def beta_s(Vd):
        V6 = Vd + 8.9
        beta = 0.02 * V6 / (exp(V6 / 5.) - 1.)
        return beta

    def alpha_c(Vd):
        V7 = Vd + 53.5
        V8 = Vd + 50.
        if Vd <= -10:
            alpha = 0.0527 * exp(V8/11. - V7/27.)
        else:
            alpha = 2 * exp(-V7 / 27.)
        return alpha

    def beta_c(Vd):
        V7 = Vd + 53.5
        if Vd <= -10:
            beta = 2. * exp(-V7 / 27.) - alpha_c(Vd)
        else:
            beta = 0.
        return beta

    def alpha_q(Ca):
        alpha = min(0.00002*Ca, 0.01)
        return alpha

    def beta_q(Ca):
        return 0.001

    def chi(Ca):
        return min(Ca/250., 1.)

    def m_inf(Vs):
        return alpha_m(Vs) / (alpha_m(Vs) + beta_m(Vs))

    def dVdt(t, V):

        Vs, Vd, n, h, s, c, q, Ca = V

        I_leak_s = g_L*(Vs - E_L)
        I_leak_d = g_L*(Vd - E_L)
        I_Na = g_Na * m_inf(Vs)**2 * h * (Vs - E_Na)
        I_DR = g_DR * n * (Vs - E_K)
        I_ds = g_c * (Vd - Vs)

        I_Ca = g_Ca * s**2 * (Vd - E_Ca)
        I_AHP = g_AHP * q * (Vd - E_K)
        I_C = g_C * c * chi(Ca) * (Vd - E_K)
        I_sd = -I_ds

        # Determine if stimulus is active
        stim_active = t > stim_start and t < stim_end

        # Apply stimulus based on location setting
        if stim_loc == 'soma':
            # Stimulus applied to soma
            dVsdt = (1./C_m) * (-I_leak_s - I_Na - I_DR +
                                I_ds/p + (I_stim/p if stim_active else 0))
            dVddt = (1./C_m) * (-I_leak_d - I_Ca - I_AHP - I_C + I_sd/(1-p))
        else:  # stim_loc == 'dendrite'
            # Stimulus applied to dendrite
            dVsdt = (1./C_m) * (-I_leak_s - I_Na - I_DR + I_ds/p)
            dVddt = (1./C_m) * (-I_leak_d - I_Ca - I_AHP - I_C +
                                I_sd/(1-p) + (I_stim/(1-p) if stim_active else 0))

        dhdt = alpha_h(Vs)*(1-h) - beta_h(Vs)*h
        dndt = alpha_n(Vs)*(1-n) - beta_n(Vs)*n
        dsdt = alpha_s(Vd)*(1-s) - beta_s(Vd)*s
        dcdt = alpha_c(Vd)*(1-c) - beta_c(Vd)*c
        dqdt = alpha_q(Ca)*(1-q) - beta_q(Ca)*q
        dCadt = -0.13*I_Ca - 0.075*Ca

        return dVsdt, dVddt, dndt, dhdt, dsdt, dcdt, dqdt, dCadt

    t_span = (0, t_dur)

    Vs0 = -68.
    Vd0 = -68.
    n0 = 0.001
    h0 = 0.999
    s0 = 0.009
    c0 = 0.007
    q0 = 0.01
    Ca0 = 0.2
    V0 = [Vs0, Vd0, n0, h0, s0, c0, q0, Ca0]

    sol = solve_ivp(dVdt, t_span, V0, max_step=0.05)

    return sol
