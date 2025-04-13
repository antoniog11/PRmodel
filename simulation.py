from solve_PRmodel import solve_PRmodel
import numpy as np
import matplotlib.pyplot as plt
import os

# Create the simulation_plots folder if it doesn't exist
os.makedirs('simulation_plots', exist_ok=True)

# Simulation parameters
t_dur = 400             # simulation duration in ms
g_c = 0                # coupling conductance between compartments (mS/cm^2)
I_input_somatic = 45    # Somatic input stimulation current (μA/cm^2)
I_input_dendrite = 40   # Dendrite input stimulation current (μA/cm^2)
stim_start = 50   # when stimulation starts (ms)
stim_end = 350    # when stimulation ends (ms)


# Run simulations with stimulation in different locations
soma_stim = solve_PRmodel(
    t_dur, g_c, I_input_somatic, stim_start, stim_end, stim_loc='soma')
dend_stim = solve_PRmodel(t_dur, g_c, I_input_dendrite, stim_start,
                          stim_end, stim_loc='dendrite')

# Create figure for comparison
fig, axs = plt.subplots(3, 2, figsize=(14, 12), sharex=True)

# Time arrays for both simulations
t_soma = soma_stim.t
t_dend = dend_stim.t

# Plot somatic voltage for both stimulation locations
axs[0, 0].plot(t_soma, soma_stim.y[0], 'b', linewidth=1.5)
axs[0, 0].set_ylabel('Somatic Voltage (mV)')
axs[0, 0].set_title('Somatic Stimulation: Somatic Voltage')
axs[0, 0].grid(True)
axs[0, 0].axvspan(stim_start, stim_end, alpha=0.2, color='gray')

axs[0, 1].plot(t_dend, dend_stim.y[0], 'b', linewidth=1.5)
axs[0, 1].set_ylabel('Somatic Voltage (mV)')
axs[0, 1].set_title('Dendritic Stimulation: Somatic Voltage')
axs[0, 1].grid(True)
axs[0, 1].axvspan(stim_start, stim_end, alpha=0.2, color='gray')

# Plot dendritic voltage for both stimulation locations
axs[1, 0].plot(t_soma, soma_stim.y[1], 'r', linewidth=1.5)
axs[1, 0].set_ylabel('Dendritic Voltage (mV)')
axs[1, 0].set_title('Somatic Stimulation: Dendritic Voltage')
axs[1, 0].grid(True)
axs[1, 0].axvspan(stim_start, stim_end, alpha=0.2, color='gray')

axs[1, 1].plot(t_dend, dend_stim.y[1], 'r', linewidth=1.5)
axs[1, 1].set_ylabel('Dendritic Voltage (mV)')
axs[1, 1].set_title('Dendritic Stimulation: Dendritic Voltage')
axs[1, 1].grid(True)
axs[1, 1].axvspan(stim_start, stim_end, alpha=0.2, color='gray')

# Plot input currents for both stimulation locations
axs[2, 0].plot(t_soma, [0 if t < stim_start or t >
               stim_end else I_input_somatic for t in t_soma], 'k', linewidth=1.5)
axs[2, 0].set_ylabel('Input Current (μA/cm²)')
axs[2, 0].set_xlabel('Time (ms)')
axs[2, 0].set_title('Somatic Input Current')
axs[2, 0].grid(True)
axs[2, 0].set_ylim(-5, I_input_somatic*1.1)

axs[2, 1].plot(t_dend, [0 if t < stim_start or t >
               stim_end else I_input_dendrite for t in t_dend], 'k', linewidth=1.5)
axs[2, 1].set_ylabel('Input Current (μA/cm²)')
axs[2, 1].set_xlabel('Time (ms)')
axs[2, 1].set_title('Dendritic Input Current')
axs[2, 1].grid(True)
axs[2, 1].set_ylim(-5, I_input_dendrite*1.1)

plt.tight_layout()

# Add annotation for parameters
param_text = f'Parameters: $g_c$ = {g_c} mS/cm$^2$, $I_{{soma}}$ = {I_input_somatic} μA/cm$^2$, $I_{{dend}}$ = {I_input_dendrite} μA/cm$^2$'
fig.text(0.5, 0.01, param_text, ha='center', fontsize=12)

# Create filename with parameters and save to simulation_plots folder
filename = os.path.join(
    'simulation_plots', f'pinsky_rinzel_comparison_gc{g_c}_Is{I_input_somatic}_Id{I_input_dendrite}.png')
plt.savefig(filename, dpi=300, bbox_inches='tight')
plt.show()

# Calculate and plot calcium concentration
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(t_soma, soma_stim.y[7], 'g', linewidth=1.5)
plt.ylabel('Ca Concentration')
plt.xlabel('Time (ms)')
plt.title('Somatic Stimulation: Ca Concentration')
plt.grid(True)
plt.axvspan(stim_start, stim_end, alpha=0.2, color='gray')

plt.subplot(1, 2, 2)
plt.plot(t_dend, dend_stim.y[7], 'g', linewidth=1.5)
plt.ylabel('Ca Concentration')
plt.xlabel('Time (ms)')
plt.title('Dendritic Stimulation: Ca Concentration')
plt.grid(True)
plt.axvspan(stim_start, stim_end, alpha=0.2, color='gray')
plt.tight_layout()

# Create calcium filename with parameters and save to simulation_plots folder
calcium_filename = os.path.join(
    'simulation_plots', f'calcium_comparison_gc{g_c}_Is{I_input_somatic}_Id{I_input_dendrite}.png')
plt.savefig(calcium_filename, dpi=300, bbox_inches='tight')
plt.show()

# Create phase plots (Vs vs Vd)
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(soma_stim.y[0], soma_stim.y[1], 'k', linewidth=0.8)
plt.xlabel('Somatic Voltage (mV)')
plt.ylabel('Dendritic Voltage (mV)')
plt.title('Phase Plot: Somatic Stimulation')
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(dend_stim.y[0], dend_stim.y[1], 'k', linewidth=0.8)
plt.xlabel('Somatic Voltage (mV)')
plt.ylabel('Dendritic Voltage (mV)')
plt.title('Phase Plot: Dendritic Stimulation')
plt.grid(True)
plt.tight_layout()

# Create phase plot filename with parameters and save to simulation_plots folder
phase_filename = os.path.join(
    'simulation_plots', f'phase_plot_comparison_gc{g_c}_Is{I_input_somatic}_Id{I_input_dendrite}.png')
plt.savefig(phase_filename, dpi=300, bbox_inches='tight')
plt.show()
