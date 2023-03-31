import QRE
import ClassicalSims
import numpy as np
import matplotlib.pyplot as plt

# Load molecules and basis sets
molecules, basis_sets = QRE.getModels()

# Compute ground state energies using QRE
qre_energies = []
for molecule in molecules:
    for basis_set in basis_sets:
        qre_energies.append(QRE.computeEnergy(molecule, basis_set))

# Compute ground state energies using ClassicalSims
classical_energies = []
for molecule in molecules:
    for basis_set in basis_sets:
        classical_energies.append(ClassicalSims.computeEnergy(molecule, basis_set))

# Compute RMSE
rmse = np.sqrt(np.mean((np.array(qre_energies) - np.array(classical_energies))**2))

def plot_rmse(qre_rmse, cs_rmse):
    # Set up the figure
    fig, ax = plt.subplots()

    # Set the x-axis label and tick marks
    ax.set_xlabel('Molecule')
    ax.set_xticks(np.arange(len(qre_rmse)))
    ax.set_xticklabels([f'Molecule {i+1}' for i in range(len(qre_rmse))])

    # Set the y-axis label and tick marks
    ax.set_ylabel('RMSE')
    ax.set_ylim([0, max(max(qre_rmse), max(cs_rmse))])
    ax.yaxis.set_ticks(np.arange(0, max(max(qre_rmse), max(cs_rmse))+0.1, 0.1))

    # Plot the RMSE values for the QRE and ClassicalSims models
    ax.plot(qre_rmse, label='QRE')
    ax.plot(cs_rmse, label='ClassicalSims')

    # Add a legend
    ax.legend()

    # Show the plot
    plt.show()

# Generate models using the QRE and ClassicalSims programs
qre_models = QRE.getModels()
cs_models = ClassicalSims.getModels()

# Calculate the RMSE values for the QRE and ClassicalSims models
qre_rmse = QRE.calculateRMSE(qre_models)
cs_rmse = ClassicalSims.calculateRMSE(cs_models)

# Plot the RMSE values
plot_rmse(qre_rmse, cs_rmse)
