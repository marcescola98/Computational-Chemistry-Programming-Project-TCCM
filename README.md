# Hartree-Fock Energy and Gradient Calculator

This project implements a simple Hartree-Fock code for small closed-shell molecules using contracted s-type Gaussian basis functions.  
It can perform SCF energy calculations and evaluate analytical energy gradients.  
The project includes two modes:
- a standard mode, where integrals are computed internally
- an extended mode, where precomputed integrals and derivative tensors are read from the input file

## Files

- `main.py` – main program driver
- `input_reader.py` – reads standard and extended input files
- `hf_kernels.f90` – builds one- and two-electron integrals
- `hf_scf.f90` – SCF routines (orthogonalization, Fock matrix, density, energy)
- `hf_gradients.f90` – gradient routines

## Requirements

- Python 3.x
- NumPy
- A Fortran compiler
- `f2py` (available through NumPy)

## Compilation

Compile the Fortran modules with:

```bash
python -m numpy.f2py -c hf_kernels.f90 -m hf_kernels
python -m numpy.f2py -c hf_scf.f90 -m hf_scf -llapack
python -m numpy.f2py -c hf_gradients.f90 -m hf_gradients
```
or

```bash
f2py -c hf_kernels.f90 -m hf_kernels
f2py -c hf_scf.f90 -m hf_scf -llapack
f2py -c hf_gradients.f90 -m hf_gradients
```

### Running the program
```md``
## Usage

Run the code with:

```bash
python main.py -i example.input
```
Alternatively you can also specify a name for the output file:

```bash
python main.py -i example.input -o result.out
```


### Input modes
```md``
## Input modes

- **Standard input**: reads molecular and basis-set data and computes all integrals internally.
- **Extended input**: reads precomputed integrals and integral derivatives from the input file.

## Notes

- Coordinates are converted from angstrom to bohr in the Python input reader.
- The code is intended for small closed-shell systems with s-type basis functions only.
- Compiled `.so` files are platform-dependent and should be rebuilt locally.

