function sol = neuron_model(t_dur, g_c, I_soma, I_dend, stim_start, stim_end)
    % Parameters from Yi et al., 2017
    C = 3;  % membrane capacitance [uF/cm^2]
    p = 0.5;
    % Maximal conductances [mS/cm^2]
    g_Na = 20;
    g_K  = 50;
    g_Ca = 0;
    g_L_s  = 2;
    g_L_d = 2;

    % Reversal potentials [mV]
    E_Na = 50;
    E_K  = -100;
    E_Ca = 120;
    E_L  = -70;

    % Gating kinetics constants
    phi_w = 0.15;
    beta_w = 0;
    gamma_w = 10;
    beta_m = -12;
    gamma_m = 18;

    n_zw = 1;
    tau_n = 15/n_zw;
    tau_h = 80/n_zw;

    % Initial conditions: [Vs, w, Vd, n, h]
    Y0 = [-65; 0.01; -65; 0.01; 1];
    tspan = [0, t_dur];
    opts = odeset('MaxStep', 0.01);

    [t, Y] = ode23(@(t, Y) dYdt(t, Y), tspan, Y0, opts);
    sol.t = t;
    sol.y = Y;

    function dY = dYdt(t, Y)
        Vs = Y(1);
        w  = Y(2);
        Vd = Y(3);
        n  = Y(4);
        h  = Y(5);
        % Activation functions
        m_inf = 0.5 * (1 + tanh((Vs - beta_m) / gamma_m));
        w_inf = 0.5 * (1 + tanh((Vs - beta_w) / gamma_w));
        tau_w = 1 / cosh((Vs - beta_w) / (2 * gamma_w));

        n_inf = 1 / (1 + exp(-(Vd + 9) / 0.5));
        h_inf = 1 / (1 + exp((Vd + 21) / 0.5));

        % Currents
        I_Na = g_Na * m_inf * (Vs - E_Na);
        I_K  = g_K  * w      * (Vs - E_K);
        I_Ls = g_L_s * (Vs - E_L);

        I_Ca = g_Ca * n * h * (Vd - E_Ca);
        I_Ld = g_L_d * (Vd - E_L);

        I_ds = g_c * (Vd - Vs);

        % Input currents
        Is = 0; Id = 0;
        if t >= stim_start && t <= stim_end
            Is = I_soma;
            Id = I_dend;
        end

        % ODEs
        dVs = (1/C) * (-I_Na - I_K - I_Ls + I_ds/p+ Is/p);
        dw  = phi_w * (w_inf - w) / tau_w;

        dVd = (1/C) * (-I_Ca - I_Ld - I_ds/(1-p) + Id/(1-p));
        dn  = (n_inf - n) / tau_n;
        dh  = (h_inf - h) / tau_h;

        dY = [dVs; dw; dVd; dn; dh];
    end
end
