import numpy as np
import matplotlib.pyplot as plt
from solve_PRmodel import solve_PRmodel

# Simulation parameters
t_dur = 400       # simulation duration in ms
g_c = 40          # coupling conductance between compartments (mS/cm^2)
I_stim = 1.5      # stimulation current (μA/cm^2)
stim_start = 50   # when stimulation starts (ms)
stim_end = 350    # when stimulation ends (ms)

# Solve the model
solution = solve_PRmodel(t_dur, g_c, I_stim, stim_start, stim_end)

# Extract the solution
t = solution.t
Vs = solution.y[0]  # Somatic voltage
Vd = solution.y[1]  # Dendritic voltage
Ca = solution.y[7]  # Calcium concentration

# Create figure with multiple subplots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

# Plot somatic voltage
ax1.plot(t, Vs, 'b', linewidth=1.5)
ax1.set_ylabel('Somatic Voltage (mV)')
ax1.set_title('Pinsky-Rinzel Model Simulation')
ax1.grid(True)

# Plot dendritic voltage
ax2.plot(t, Vd, 'r', linewidth=1.5)
ax2.set_ylabel('Dendritic Voltage (mV)')
ax2.grid(True)

# Plot calcium concentration
ax3.plot(t, Ca, 'g', linewidth=1.5)
ax3.set_xlabel('Time (ms)')
ax3.set_ylabel('Ca Concentration')
ax3.grid(True)

# Add a shaded region to indicate when the stimulus is on
for ax in [ax1, ax2, ax3]:
    ax.axvspan(stim_start, stim_end, alpha=0.2, color='gray')
    ax.text(stim_start + 5, ax.get_ylim()[1] * 0.9, 'Stimulus On', fontsize=10)

# Add annotation for parameters
param_text = f'Parameters: $g_c$ = {g_c} mS/cm$^2$, $I_{{stim}}$ = {I_stim} μA/cm$^2$'
fig.text(0.5, 0.01, param_text, ha='center', fontsize=12)

plt.tight_layout()

# Save the figure
plt.savefig('pinsky_rinzel_simulation.png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

# Let's also create a phase plot of Vs vs Vd
plt.figure(figsize=(8, 8))
plt.plot(Vs, Vd, 'k', linewidth=0.8)
plt.xlabel('Somatic Voltage (mV)')
plt.ylabel('Dendritic Voltage (mV)')
plt.title('Phase Plot: Somatic vs Dendritic Voltage')
plt.grid(True)
plt.savefig('phase_plot.png', dpi=300, bbox_inches='tight')
plt.show()

# Create a zoomed view of a specific action potential
# Let's find where action potentials occur
ap_threshold = -20  # threshold to detect action potentials
ap_indices = np.where(Vs > ap_threshold)[0]

if len(ap_indices) > 0:
    # Find the first action potential after stimulus onset
    first_ap_after_stim = next(
        (i for i in ap_indices if t[i] > stim_start), None)

    if first_ap_after_stim is not None:
        # Create a window around this action potential
        window_start = max(0, first_ap_after_stim - 100)
        window_end = min(len(t) - 1, first_ap_after_stim + 200)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        ax1.plot(t[window_start:window_end],
                 Vs[window_start:window_end], 'b', linewidth=1.5)
        ax1.set_ylabel('Somatic Voltage (mV)')
        ax1.set_title('Zoomed View of Action Potential')
        ax1.grid(True)

        ax2.plot(t[window_start:window_end],
                 Vd[window_start:window_end], 'r', linewidth=1.5)
        ax2.set_ylabel('Dendritic Voltage (mV)')
        ax2.set_xlabel('Time (ms)')
        ax2.grid(True)

        plt.tight_layout()
        plt.savefig('zoomed_action_potential.png',
                    dpi=300, bbox_inches='tight')
        plt.show()
