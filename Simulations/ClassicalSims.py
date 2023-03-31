import re
import numpy as np
import scipy.linalg as la
from pyscf import gto, scf, dft
import psi4
import matplotlib.pyplot as plt


def get_energy_calculation(molecule, basis, method):
    # Define the molecule
    mol = gto.M(atom=molecule, basis=basis)
    if method == 'HF':
        mf = scf.RHF(mol)
    else:
        mf = dft.RKS(mol)
        mf.xc = method
    mf.kernel()
    energy = mf.energy_tot()
    # Print the total energy
    return energy

def get_wavefunction(molecule, basis, method):
    """
    Calculates the wavefunction of a molecule using the given basis set and method.

    Parameters:
    molecule (str): A string that specifies the molecule in XYZ format.
    basis (str): The name of the basis set to use.
    method (str): The name of the method to use (e.g., 'HF', 'B3LYP', etc.).

    Returns:
    wavefunction (ndarray): A numpy array representing the wavefunction.
    """
    mol = psi4.geometry(molecule)
    psi4.set_options({'basis': basis, 'reference': 'rhf', 'dft_functional': method})
    wavefunction = psi4.properties('SCF DIPOLE', return_wfn=True)[1]
    return wavefunction

def get_density_matrix(molecule, basis, method):
    """
    Calculates the density matrix of a molecule using the given basis set and method.

    Parameters:
    molecule (str): A string that specifies the molecule in XYZ format.
    basis (str): The name of the basis set to use.
    method (str): The name of the method to use (e.g., 'HF', 'B3LYP', etc.).

    Returns:
    density_matrix (ndarray): A numpy array representing the density matrix.
    """
    mol = gto.M(atom=molecule, basis=basis)
    if method == 'HF':
        mf = scf.RHF(mol)
    else:
        mf = dft.RKS(mol)
        mf.xc = method
    mf.kernel()
    density_matrix = mf.make_rdm1()
    return density_matrix

def get_molecular_orbitals(molecule, basis, method):
    """
    Calculates the molecular orbitals of a molecule using the given basis set and method.

    Parameters:
    molecule (str): A string that specifies the molecule in XYZ format.
    basis (str): The name of the basis set to use.
    method (str): The name of the method to use (e.g., 'HF', 'B3LYP', etc.).

    Returns:
    molecular_orbitals (ndarray): A numpy array representing the molecular orbitals.
    """
    mol = gto.M(atom=molecule, basis=basis)
    if method == 'HF':
        mf = scf.RHF(mol)
    else:
        mf = dft.RKS(mol)
        mf.xc = method
    mf.kernel()
    mo_coeff = mf.mo_coeff
    molecular_orbitals = mo_coeff[:,:mol.nelec]
    return molecular_orbitals

def get_spectroscopic_data(molecule, basis, method):
    """
    Calculates the spectroscopic data of a molecule using the given basis set and method.

    Parameters:
    molecule (str): A string that specifies the molecule in XYZ format.
    basis (str): The name of the basis set to use.
    method (str): The name of the method to use (e.g., 'HF', 'B3LYP', etc.).

    Returns:
    spectroscopic_data (dict): A dictionary containing the spectroscopic data.
    """
    mol = gto.M(atom=molecule, basis=basis)
    if method == 'HF':
        mf = scf.RHF(mol)
    else:
        mf = dft.RKS(mol)
        mf.xc = method
    mf.kernel()
    mo_energy = mf.mo_energy
    eigvals, eigvecs = eigh(mf.get_hcore() + mf.get_jk()[0])
    spectroscopic_data = {
        'energies': mo_energy,
        'oscillator_strengths': 2 / 3 * h * c * np.abs(mf.dip_moment()[0])**2 / eigvals**3,
        'transition_dipoles': eigvals,
'molecular_orbitals': mf.mo_coeff,}
return spectroscopic_data

def getModels():
    # Define the molecules to simulate
    molecules = [
        ('H2', 'sto-3g', 'HF'),
        ('LiH', 'sto-3g', 'HF'),
        ('H2O', 'sto-3g', 'HF'),
        ('H2O', '6-31g', 'B3LYP'),
        ('CO2', 'sto-3g', 'HF'),
        ('NH3', 'sto-3g', 'HF')
    ]

    models = []

    # Calculate the electronic structure of each molecule
    for molecule in molecules:
        # Get the energy
        energy = get_energy_calculation(*molecule)
        # Get the wavefunction
        wavefunction = get_wavefunction(*molecule)
        # Get the density matrix
        density_matrix = get_density_matrix(*molecule)
        # Get the molecular orbitals
        molecular_orbitals = get_molecular_orbitals(*molecule)
        # Get the spectroscopic data
        spectroscopic_data = get_spectroscopic_data(*molecule)

        # Store the data in a dictionary
        model = {
            'molecule': molecule[0],
            'basis': molecule[1],
            'method': molecule[2],
            'energy': energy,
            'wavefunction': wavefunction,
            'density_matrix': density_matrix,
            'molecular_orbitals': molecular_orbitals,
            'spectroscopic_data': spectroscopic_data
        }

        models.append(model)

    return models

def energy_visualization(models):
    ""Visualizes the energies of the molecules in the given models.""
    for model in models:
        plt.figure()
        plt.title(f"{model['molecule']} Energy")
        plt.bar(['HF', 'B3LYP'], model['energy'])
        plt.xlabel("Method")
        plt.ylabel("Energy (Hartrees)")
        plt.savefig(f"{model['molecule']}_energy.png")
        plt.close()

def molecular_orbital_visualization(models):
    ""Visualizes the molecular orbitals of the molecules in the given models.""
    for model in models:
        num_orbitals = len(model['molecular_orbitals'])
        num_rows = int(np.ceil(num_orbitals/2))
        fig, axs = plt.subplots(num_rows, 2, figsize=(8, num_rows*3))
        fig.suptitle(f"{model['molecule']} Molecular Orbitals")
        for i, orbital in enumerate(model['molecular_orbitals']):
            ax = axs[i//2, i%2]
            ax.set_title(f"Orbital {i+1}")
            ax.contourf(orbital, cmap='coolwarm')
            ax.set_aspect('equal')
            ax.axis('off')
        plt.savefig(f"{model['molecule']}_molecular_orbitals.png")
        plt.close()

def density_matrix_visualization(models):
    ""Visualizes the density matrices of the molecules in the given models.""
    for model in models:
        plt.figure()
        plt.title(f"{model['molecule']} Density Matrix")
        plt.imshow(model['density_matrix'], cmap='coolwarm')
        plt.colorbar()
        plt.xlabel("Atom Index")
        plt.ylabel("Atom Index")
        plt.savefig(f"{model['molecule']}_density_matrix.png")
        plt.close()

def spectroscopic_data_visualization(models):
    """Visualizes the spectroscopic data of the molecules in the given models."""
    for model in models:
        fig, axs = plt.subplots(2, 2, figsize=(8, 6))
        fig.suptitle(f"{model['molecule']} Spectroscopic Data")
        axs[0, 0].set_title("IR Spectrum")
        axs[0, 0].plot(*model['spectroscopic_data']['IR'])
        axs[0, 0].set_xlabel("Wavenumber (cm$^{-1}$)")
        axs[0, 0].set_ylabel("Transmittance")
        axs[0, 1].set_title("UV-Vis Spectrum")
        axs[0, 1].plot(*model['spectroscopic_data']['UV-Vis'])
        axs[0, 1].set_xlabel("Wavelength (nm)")
        axs[0, 1].set_ylabel("Absorbance")
        axs[1, 0].set_title("Ramam Spectrum")
        axs[1, 0].plot(*model['spectroscopic_data']['Raman'])
        axs[1, 0].set_xlabel("Wavenumber (cm$^{-1}$)")
        axs[1, 0].set_ylabel("Intensity")
        axs[1, 1].set_title("NMR Spectrum")
        axs[1, 1].plot(*model['spectroscopic_data']['NMR'])
        axs[1, 1].set_xlabel("Chemical Shift (ppm)")
        axs[1, 1].set_ylabel("Signal Intensity")
        plt.tight_layout()
        plt.savefig(f"{model['molecule']}_spectroscopic_data.png")
        plt.close()

}
