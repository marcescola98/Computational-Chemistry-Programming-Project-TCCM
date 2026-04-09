import numpy as np
import time
import os
import argparse
from input_reader import read_infile, read_extended
import hf_kernels as hf # type: ignore
import hf_scf as hf_scf # type: ignore
import hf_gradients as grad # type: ignore

def main():

    # --------------------------------------------------------------#
    # Parse command-line arguments and determine input/output files #
    # --------------------------------------------------------------#

    parser = argparse.ArgumentParser(
        prog = 'HF Gradient Calculator',
        usage='%(prog)s [options]',
        description='A HF gradient calculator program for the Computational Chemistry Programming Project, a subject from the TCCM MsC.' \
        'This program uses Python for structure and Fortran90 for calculations in order to calculate the HF gradient from different simple molecules. '
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input file name. Required!'
    )
    parser.add_argument(
        '-o','--output',
        required=False,
        help='Output file name. Not required. In case the user doesn\'t provide an output flag, the program will automatically use the base name of the input file to create the output\'s file name'
    )

    args = parser.parse_args()
    base_name, _ = os.path.splitext(args.input)
    
    if not args.output:
        args.output = base_name + '.out'

    input_file = args.input
    output_file = args.output
    initial_time = time.time()


    #-----------------------------------------------------------------------------------------------#
    # VARIABLES THAT WILL BE USED THROUGHOUT THE PROGRAM. IF NEEDED, THEY CAN BE EDITED BY THE USER #
    #-----------------------------------------------------------------------------------------------#

    old_energy = 0.0
    convergence_deltaE = 1e-8
    convergence_deltaP = 1e-6
    max_iter = 50
    converged = False



    with open(output_file, 'w') as f_out:
        f_out.write("-" * 26 +"\n")
        f_out.write(f"| HF GRADIENT CALCULATOR |\n")
        f_out.write("-" * 26 +"\n")
        
        if 'extended' in input_file.lower():

            # --------------------------------------------------------------#
            # Extended-input mode:                                          #
            # use precomputed integrals for SCF and precomputed derivative  #
            # tensors for gradient assembly                                 #
            # --------------------------------------------------------------#

            data = read_extended(input_file)

            n_electrons = np.sum(data['Z']) - data['charge']
            nocc = n_electrons // 2


            f_out.write(f"Data will be extracted from the next input file: {input_file}\n\n")
            f_out.write(f"System information:\n")
            f_out.write(f"- {data['natoms']} atoms\n")
            f_out.write(f"- {n_electrons} electrons\n")
            f_out.write(f"- {nocc} occupied orbitals\n")
            f_out.write(f"- {data['nbasis']} basis functions\n")
            f_out.write(f"- Overall charge: {data['charge']}\n\n")
            f_out.write(f"Coordinates (in Bohr):\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{data['atom_coords'][i,j]:>12.8f}")
                f_out.write("\n")
            f_out.write("\n")


            enuc = hf.hf_kernels.nuclear_repulsion(
                data['natoms'],
                data['Z'],
                data['atom_coords']
            )


            f_out.write(f'With a nuclear repulsion energy of: {enuc:.8f} Hartree\n\n')
            f_out.write(f"Pre-computed matrices and integrals read from input file. Starting SCF... \n\n")


            H_core = data['T_ext'] + data['V_ext']
            X = hf_scf.hf_scf.orthogonalizer(data['nbasis'], data['S_ext'])
            P = np.zeros((data['nbasis'], data['nbasis']))


            f_out.write(f"{'Iter':<5} | {'Electronic Energy (Ha)':<25} | {'Max delta P':<16}| {'Delta E':<15}\n")
            f_out.write("-" * 50 + "\n")

            # -------------------------------------------------------------- #
            # Self-consistent field (SCF) procedure                          #
            # Start from P = 0, build Fock matrix, solve Roothaan equations, #
            # update density, and iterate until both energy and density      #
            # changes are below threshold                                    #
            # -------------------------------------------------------------- #

            for i in range(max_iter):

            # Step 1: Build the G matrix using the current density matrix P and the computed ERI

                G = hf_scf.hf_scf.compute_g_matrix(data['nbasis'], P, data['ERI_ext'])
            
            # Step 2: Form the Fock Matrix

                F = H_core + G
            
            # Step 3: Solve Fock Equations (Diagonalize and Back-Transform)

                C, eps = hf_scf.hf_scf.solve_fock_equations(data['nbasis'], F, X)
                
            # Step 4: Build New Density Matrix

                P_new = hf_scf.hf_scf.compute_density_matrix(data['nbasis'], nocc, C)
                
            # Step 5: Calculate Energy

                energy = hf_scf.hf_scf.calculate_electronic_energy(data['nbasis'], P_new, H_core, F)
                
            # Step 6: Check Convergence (Energy and Density Matrix)

                delta_e = abs(energy - old_energy)
                delta_p = np.max(np.abs(P_new - P))

                f_out.write(f"{i+1:<5} | {energy:<25.8f} | {delta_p:<15.8f} | {delta_e:<15.10f}\n")

            # If both energy and density changes are below the thresholds, we consider SCF converged and break the loop

                if delta_e < convergence_deltaE and delta_p < convergence_deltaP:
                    total_energy = energy + enuc
                    converged = True
                    P = P_new
                    f_out.write(f"\nSCF Converged in {i+1} iterations. \n\n")
                    break
                
                P = P_new
                old_energy = energy
            
            if not converged:
            
                total_energy = energy + enuc
                f_out.write(f"\nSCF did not converge in {max_iter} iterations. \n\n")
            


            f_out.write(f"Final SCF electronic energy = {energy:.8f} Hartree \n")
            f_out.write(f"Total energy (electronic + nuclear repulsion) = {total_energy:.8f} Hartree \n")
            f_out.write("Final SCF density matrix P:\n")
            for i in range(data['nbasis']):
                f_out.write(f"{i+1:<3}   ")
                for j in range(data['nbasis']):
                    f_out.write(f"{P[i,j]:>14.8f}        ")
                f_out.write("\n")
            f_out.write("\n")
            f_out.write("Final SCF Fock matrix F:\n")
            for i in range(data['nbasis']):
                f_out.write(f"{i+1:<3}   ")
                for j in range(data['nbasis']):
                    f_out.write(f"{F[i,j]:>14.8f}        ")
                f_out.write("\n")
            f_out.write("\n")
            f_out.write("\nFinal orbital energies (Hartree):\n")
            for i in range(data['nbasis']):
                f_out.write(f"Orbital {i+1:<3} {eps[i]: .8f}\n")
            f_out.write(f"\nFinal orbital coefficients:\n")
            for i in range(data['nbasis']):
                f_out.write(f"Orbital {i+1:<3}   ")
                for j in range(data['nbasis']):
                    f_out.write(f"{C[j,i]:>14.8f}        ")
                f_out.write("\n")
            
            #------------------------------#
            # Gradient calculation section #
            #------------------------------#

            f_out.write(f"Calculating gradients... \n\n")
            
            # First, we compute the W matrix, which is needed for the overlap gradient contribution.

            W = grad.hf_gradients.compute_w_matrix(
                data['nbasis'],
                nocc,
                eps,
                C
            )

            # Then, we compute the gradient contribution from nuclear repulsion, which is independent of the SCF procedure and can be computed directly from the atomic coordinates and charges.

            vnn_grad = grad.hf_gradients.compute_vnn_gradient(
                data['natoms'],
                data['Z'],
                data['atom_coords']
            )

            f_out.write(f"Gradient contribution from nuclear repulsion: \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{vnn_grad[i,j]:>14.8f}        ")
                f_out.write("\n")
            
            f_out.write("\n")

            # Then, we compute the gradient contribution from the overlap matrix.

            overlap_grad = np.zeros((data['natoms'], 3))
            for A in range(data['natoms']):
                for xyz in range(3):
                    overlap_grad[A, xyz] = -np.sum(W * data['dS_ext'][:, :, A, xyz])


            f_out.write(f"Gradient contribution from overlap (Lagrangian): \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{overlap_grad[i,j]:>14.8f}        ")
                f_out.write("\n")

            # Next, we compute the gradient contribution from the kinetic energy matrix.

            kinetic_grad = np.zeros((data['natoms'], 3))
            for A in range(data['natoms']):
                for xyz in range(3):
                    kinetic_grad[A, xyz] = np.sum(P_new * data['dT_ext'][:, :, A, xyz])


            f_out.write(f"\nGradient contribution from kinetic energy: \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{kinetic_grad[i,j]:>14.8f}        ")
                f_out.write("\n")

            # Next, we compute the gradient contribution from the nuclear attraction matrix.

            vne_grad = np.zeros((data['natoms'], 3))
            for A in range(data['natoms']):
                for xyz in range(3):
                    vne_grad[A, xyz] = np.sum(P_new * data['dV_ext'][:, :, A, xyz])


            f_out.write(f"\nGradient contribution from nuclear attraction: \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{vne_grad[i,j]:>14.8f}        ")
                f_out.write("\n")

            # Finally, we compute the gradient contribution from the electron repulsion integrals (ERI). This is the most computationally intensive part, as it involves a 4-index summation over the basis functions and the density matrix.

            eri_grad = np.zeros((data['natoms'], 3))
            nb = data['nbasis']

            for A in range(data['natoms']):
                for xyz in range(3):
                    val = 0.0
                    for mu in range(nb):
                        for nu in range(nb):
                            for la in range(nb):
                                for si in range(nb):
                                    d = data['dERI_ext'][mu, nu, la, si, A, xyz]
                                    val += 0.5  * P_new[mu, nu] * P_new[la, si] * d
                                    val -= 0.25 * P_new[mu, la] * P_new[nu, si] * d
                    eri_grad[A, xyz] = val


            f_out.write(f"\nGradient contribution from electron repulsion integrals (ERI): \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{eri_grad[i,j]:>14.8f}        ")
                f_out.write("\n")

            # After computing all individual contributions, we sum them up to get the total gradient for each atom and check for translational invariance by summing all atomic gradients.
            # The sum of all atomic gradients should be close to zero.

            total_grad = overlap_grad + kinetic_grad + vne_grad + vnn_grad + eri_grad
            grad_sum = np.sum(total_grad, axis=0)


            f_out.write(f"\nTotal gradient: \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{total_grad[i,j]:>14.8f}        ")
                f_out.write("\n")

            f_out.write("\nTranslational invariance check (sum of all atomic gradients):\n\n")
            f_out.write(f"{grad_sum[0]:>14.8f} {grad_sum[1]:>14.8f} {grad_sum[2]:>14.8f}\n")

            time_taken = time.time() - initial_time

            f_out.write(f"\nAll calculations completed.\n")
            f_out.write(f"\nTotal time taken: {time_taken:.3f} seconds\n")

            print(f"Calculations completed. Results written to {output_file}")



        else:
            
            # ----------------------------------------------------------------- #
            # Standard input mode:                                              #
            # use only the basic molecular information and basis set data       #
            # to compute all necessary integrals for SCF and gradient assembly  #
            # ----------------------------------------------------------------- #

            data = read_infile(input_file)
            f_out.write(f"Data will be extracted from the next input file: {input_file}\n\n")

            n_electrons = np.sum(data['Z']) - data['charge']
            nocc = n_electrons // 2

            f_out.write(f"System information:\n")
            f_out.write(f"- {data['natoms']} atoms\n")
            f_out.write(f"- {n_electrons} electrons\n")
            f_out.write(f"- {nocc} occupied orbitals\n")
            f_out.write(f"- {data['nbasis']} basis functions\n")
            f_out.write(f"- Overall charge: {data['charge']}\n\n")
            f_out.write(f"Coordinates (in Bohr):\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{data['atom_coords'][i,j]:>12.8f}")
                f_out.write("\n")
            f_out.write("\n")
            

            enuc = hf.hf_kernels.nuclear_repulsion(
                data['natoms'],
                data['Z'],
                data['atom_coords']
            )

            f_out.write(f'With a nuclear repulsion energy of: {enuc:.8f} Hartree\n\n')

            f_out.write(f"Building necessary matrices for SCF...\n\n")


            # ---------------------------- #
            # Integral computation section #
            # ---------------------------- #


            S = hf.hf_kernels.compute_overlap_matrix(
                data['nbasis'],
                data['n_primitives'],
                data['exponents'],
                data['coefficients'],
                data['atom_coords'],
                data['atom_map']
            )

            f_out.write(f"Overlap matrix S computed successfully.\n\n")

            T = hf.hf_kernels.compute_kinetic_matrix(
                data['nbasis'],
                data['n_primitives'],
                data['exponents'],
                data['coefficients'],
                data['atom_coords'],
                data['atom_map']
            )

            f_out.write(f"Kinetic energy matrix T computed successfully.\n\n")

            V = hf.hf_kernels.compute_nuclear_atraction_matrix(
                data['nbasis'],
                data['n_primitives'],
                data['exponents'],
                data['coefficients'],
                data['atom_coords'],
                data['atom_map'],
                data['Z']
            )

            f_out.write(f"Nuclear attraction matrix V computed successfully.\n\n")

            eri = hf.hf_kernels.compute_eri(
                data['nbasis'],
                data['n_primitives'],
                data['exponents'],
                data['coefficients'],
                data['atom_coords'],
                data['atom_map']
            )

            f_out.write(f"Electron repulsion integrals (ERI) computed successfully.\n\n")

            f_out.write(f"All matrices correctly computed. Starting SCF... \n\n")

            # -------------------------------------------------------------- #
            # Self-consistent field (SCF) procedure                          #
            # Start from P = 0, build Fock matrix, solve Roothaan equations, #
            # update density, and iterate until both energy and density      #
            # changes are below threshold                                    #
            # -------------------------------------------------------------- #

            H_core = T + V

            n_electrons = np.sum(data['Z']) - data['charge']
            nocc = n_electrons // 2

            X = hf_scf.hf_scf.orthogonalizer(data['nbasis'], S)
            P = np.zeros((data['nbasis'], data['nbasis']))

            f_out.write(f"{'Iter':<5} | {'Electronic Energy (Ha)':<25} | {'Max delta P':<16}| {'Delta E':<15}\n")
            f_out.write("-" * 70 + "\n")

            for i in range(max_iter):

            # Step 1: Build the G matrix using the current density matrix P and the computed ERI

                G = hf_scf.hf_scf.compute_g_matrix(data['nbasis'], P, eri)
            
            # Step 2: Form the Fock Matrix

                F = H_core + G
            
            # Step 3: Solve Fock Equations (Diagonalize and Back-Transform)

                C, eps = hf_scf.hf_scf.solve_fock_equations(data['nbasis'], F, X)
                
            # Step 4: Build New Density Matrix

                P_new = hf_scf.hf_scf.compute_density_matrix(data['nbasis'], nocc, C)
                
            # Step 5: Calculate Energy

                energy = hf_scf.hf_scf.calculate_electronic_energy(data['nbasis'], P_new, H_core, F)
                
            # Step 6: Check Convergence (Energy and Density Matrix)

                delta_e = abs(energy - old_energy)
                delta_p=np.max(np.abs(P_new - P))

                f_out.write(f"{i+1:<5} | {energy:<25.8f} | {delta_p:<15.8f} | {delta_e:<15.10f}\n")

            # If both energy and density changes are below the thresholds, we consider SCF converged and break the loop
               
                if delta_e < convergence_deltaE and delta_p < convergence_deltaP:
                    total_energy = energy + enuc
                    P = P_new
                    converged = True
                    f_out.write(f"\nSCF Converged in {i+1} iterations. \n\n")
                    break
                
                P = P_new
                old_energy = energy

            if not converged:

                total_energy = energy + enuc
                f_out.write(f"\nSCF did not converge in {max_iter} iterations. \n\n")
            

            f_out.write(f"Final SCF electronic energy = {energy:.8f} Hartree \n")
            f_out.write(f"Total energy (electronic + nuclear repulsion) = {total_energy:.8f} Hartree \n")
            f_out.write("Final SCF density matrix P:\n")
            for i in range(data['nbasis']):
                f_out.write(f"{i+1:<3}   ")
                for j in range(data['nbasis']):
                    f_out.write(f"{P[i,j]:>14.8f}        ")
                f_out.write("\n")
            f_out.write("\n")
            f_out.write("Final SCF Fock matrix F:\n")
            for i in range(data['nbasis']):
                f_out.write(f"{i+1:<3}   ")
                for j in range(data['nbasis']):
                    f_out.write(f"{F[i,j]:>14.8f}        ")
                f_out.write("\n")
            f_out.write("\n")
            f_out.write("\nFinal orbital energies (Hartree):\n")
            for i in range(data['nbasis']):
                f_out.write(f"Orbital {i+1:<3} {eps[i]: .8f}\n")
            f_out.write(f"\nFinal orbital coefficients:\n")
            for i in range(data['nbasis']):
                f_out.write(f"Orbital {i+1:<3}   ")
                for j in range(data['nbasis']):
                    f_out.write(f"{C[j,i]:>14.8f}        ")
                f_out.write("\n")

            
            f_out.write(f"\nCalculating gradients... \n\n")

            
            #------------------------------#
            # Gradient calculation section #
            #------------------------------#


            # First we compute the W matrix, which is needed for the overlap gradient contribution.

            W = grad.hf_gradients.compute_w_matrix(
                data['nbasis'],
                nocc,
                eps,
                C
            )

            # Second, we compute the gradient contribution from nuclear repulsion.

            vnn_grad = grad.hf_gradients.compute_vnn_gradient(
                data['natoms'],
                data['Z'],
                data['atom_coords']
            )

            f_out.write(f"Gradient contribution from nuclear repulsion: \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{vnn_grad[i,j]:>14.8f}        ")
                f_out.write("\n")
            
            f_out.write("\n")

            # Then, we compute the gradient contribution from the overlap matrix.

            overlap_grad = grad.hf_gradients.compute_overlap_grad_contrib(
                data['nbasis'], 
                data['natoms'], 
                data['n_primitives'], 
                data['exponents'], 
                data['coefficients'], 
                data['atom_coords'], 
                data['atom_map'], 
                W
            )

            f_out.write(f"Gradient contribution from overlap (Lagrangian): \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{overlap_grad[i,j]:>14.8f}        ")
                f_out.write("\n")

            # Then, we compute the gradient contribution from the kinetic energy matrix.

            kinetic_grad = grad.hf_gradients.compute_kinetic_grad_contrib(
                data['nbasis'],
                data['natoms'],
                data['n_primitives'],
                data['exponents'],
                data['coefficients'],
                data['atom_coords'],
                data['atom_map'],
                P_new
            )

            f_out.write(f"\nGradient contribution from kinetic energy: \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{kinetic_grad[i,j]:>14.8f}        ")
                f_out.write("\n")

            # Next, we compute the gradient contribution from the nuclear attraction matrix.

            vne_grad = grad.hf_gradients.compute_nuclear_attraction_grad_contrib(
                data['nbasis'],
                data['natoms'],
                data['n_primitives'],
                data['exponents'],
                data['coefficients'],
                data['atom_coords'],
                data['atom_map'],
                data['Z'],
                P_new
            )

            f_out.write(f"\nGradient contribution from nuclear attraction: \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{vne_grad[i,j]:>14.8f}        ")
                f_out.write("\n")

            # Finally, we compute the gradient contribution from the electron repulsion integrals (ERI).

            eri_grad = grad.hf_gradients.compute_eri_grad_contrib(
                data['nbasis'],
                data['natoms'],
                data['n_primitives'],
                data['exponents'],
                data['coefficients'],
                data['atom_coords'],
                data['atom_map'],
                P_new
            )

            f_out.write(f"\nGradient contribution from electron repulsion integrals (ERI): \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{eri_grad[i,j]:>14.8f}        ")
                f_out.write("\n")

            # After computing all individual contributions, we sum them up to get the total gradient for each atom and check for translational invariance by summing all atomic gradients.
            # The sum of all atomic gradients should be close to zero.

            total_grad = overlap_grad + kinetic_grad + vne_grad + vnn_grad + eri_grad
            grad_sum = np.sum(total_grad, axis=0)


            f_out.write(f"\nTotal gradient: \n\n")
            for i in range(data['natoms']):
                f_out.write(f"{data['labels'][i]}   ")
                for j in range (0,3):
                    f_out.write(f"{total_grad[i,j]:>14.8f}        ")
                f_out.write("\n")

            f_out.write("\nTranslational invariance check (sum of all atomic gradients):\n\n")
            f_out.write(f"{grad_sum[0]:>14.8f} {grad_sum[1]:>14.8f} {grad_sum[2]:>14.8f}\n")


            time_taken = time.time() - initial_time

            f_out.write(f"\nAll calculations completed.\n")
            f_out.write(f"\nTotal time taken: {time_taken:.3f} seconds\n")

            print(f"Calculations completed. Results written to {output_file}")


if __name__ == '__main__':
    main()