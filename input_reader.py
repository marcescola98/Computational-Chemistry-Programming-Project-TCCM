import numpy as np

def read_infile(input_name, angstrom_to_bohr = 1.8897259886, return_index=False):

    """This function reads the input file containing the molecular information and basis set data.
      It parses the file line by line, looking for specific keywords to identify different sections of the input.
      The data is stored in a dictionary for easy access. If return_index is True, it also returns the list of
      lines and the current index in the parsing process, which can be useful for further parsing in the read_extended function."""
    
    with open(input_name, 'r') as f:
        
        #We filter out empty lines and strip leading and trailing whitespace from each line
        
        lines = [line.strip() for line in f.readlines() if line.strip()]

    data = {}

    i = 0
    while i < len(lines):
        
        line = lines[i]

        #---------- Number of atoms ---------

        if 'number of atoms' in line.lower():
                data['natoms'] = int(lines[i+1])
                i += 2
        
        #--------- Atom labels, atomic numbers, and coordinates ---------

        elif 'atom labels' in line.lower():
            i += 1
            
            natoms = data['natoms']
            data['labels'] = []
            data['Z'] = np.zeros((natoms), dtype=np.int32)
            data['atom_coords'] = np.zeros((natoms, 3), dtype=np.float64)

            for j in range(natoms):
                parts = lines[i].split()
                data['labels'].append(parts[0])
                data['Z'][j] = int(parts[1])
                data['atom_coords'][j] = [float(k) * angstrom_to_bohr for k in parts[2:5] ]
                i += 1
                
        #--------- Overall charge ---------

        elif 'overall charge' in line.lower():
            data['charge'] = int(lines[i+1])
            i += 2
            
        #--------- Basis set information ---------

        elif 'basis funcs' in line.lower():
            data['nbasis'] = int(lines[i+1])
            i += 2

        #--------- Maximum number of primitives per basis function ---------
           
        elif 'maximum number of primitives' in line.lower():
            data['max_prim'] = int(lines[i+1])
            i += 2

        #--------- Basis set data ---------

        elif 'basis set' in line.lower():
            i += 1
            
            nbasis = data['nbasis']
            max_prim = data['max_prim']

            data['atom_map'] = np.zeros((nbasis), dtype=np.int32)
            data['n_primitives'] = np.zeros((nbasis), dtype=np.int32)
            data['exponents'] = np.zeros((nbasis, max_prim), dtype=np.float64)
            data['coefficients'] = np.zeros((nbasis, max_prim), dtype=np.float64)

            for j in range(nbasis):
                header = lines[i].split()

                # Store the atom index to which this basis function is associated
                data['atom_map'][j] = int(header[2])
                i += 1

                n_p = int(lines[i])
                data['n_primitives'][j] = n_p
                i += 1

                for k in range(n_p):
                    zeta, djk = map(float, lines[i].split())
                    data['exponents'][j, k] = zeta
                    data['coefficients'][j, k] = djk
                    i += 1

            if return_index:
                return data, lines, i
            return data
        
        else:
            i += 1

    if return_index:
        return data, lines, i
    return data

def read_extended(input_name,angstrom_to_bohr = 1.8897259886):
    
    """This function reads the extended input file which contains not only the basic molecular information and basis set data,
    but also the integrals and their derivatives. It assumes that the basic information has already been read by the read_infile function,
    and it continues parsing from where read_infile left off."""
    
    data, lines, i = read_infile(input_name, angstrom_to_bohr, return_index=True)
    
    nbasis = data['nbasis']
    natoms = data['natoms']

    # ---------- Overlap integrals ----------

    while i < len(lines) and 'overlap integrals' not in lines[i].lower():
        i += 1
    i += 1
    nS = int(lines[i])
    i += 1

    data['S_ext'] = np.zeros((nbasis, nbasis), dtype=np.float64)

    for _ in range(nS):
        mu, nu, val = lines[i].split()
        mu = int(mu) - 1
        nu = int(nu) - 1
        val = float(val)

        data['S_ext'][mu, nu] = val
        data['S_ext'][nu, mu] = val
        i += 1

    # ---------- Kinetic integrals ----------

    while i < len(lines) and 'kinetic integrals' not in lines[i].lower():
        i += 1
    i += 1
    nT = int(lines[i])
    i += 1

    data['T_ext'] = np.zeros((nbasis, nbasis), dtype=np.float64)

    for _ in range(nT):
        mu, nu, val = lines[i].split()
        mu = int(mu) - 1
        nu = int(nu) - 1
        val = float(val)

        data['T_ext'][mu, nu] = val
        data['T_ext'][nu, mu] = val
        i += 1

    # ---------- Nuclear attraction integrals ----------

    while i < len(lines) and 'nuclear attraction integrals' not in lines[i].lower():
        i += 1
    i += 1
    nV = int(lines[i])
    i += 1

    data['V_ext'] = np.zeros((nbasis, nbasis), dtype=np.float64)

    for _ in range(nV):
        mu, nu, val = lines[i].split()
        mu = int(mu) - 1
        nu = int(nu) - 1
        val = float(val)

        data['V_ext'][mu, nu] = val
        data['V_ext'][nu, mu] = val
        i += 1

    # ---------- Two-electron integrals ----------

    while i < len(lines) and 'two-electron integrals' not in lines[i].lower():
        i += 1
    i += 1
    nERI = int(lines[i])
    i += 1

    data['ERI_ext'] = np.zeros((nbasis, nbasis, nbasis, nbasis), dtype=np.float64)

    for _ in range(nERI):
        mu, nu, la, si, val = lines[i].split()
        mu = int(mu) - 1
        nu = int(nu) - 1
        la = int(la) - 1
        si = int(si) - 1
        val = float(val)

        # Fill all equivalent permutations related by ERI symmetry

        perms = [
            (mu, nu, la, si),
            (mu, nu, si, la),
            (nu, mu, la, si),
            (nu, mu, si, la),
            (la, si, mu, nu),
            (la, si, nu, mu),
            (si, la, mu, nu),
            (si, la, nu, mu),
        ]
        for p in perms:
            data['ERI_ext'][p] = val

        i += 1

    # ----------  Overlap derivative vectors ----------

    while i < len(lines) and 'derivatives of overlap integrals' not in lines[i].lower():
        i += 1
    i += 1
    ndS = int(lines[i])
    i += 3

    data['dS_ext'] = np.zeros((nbasis, nbasis, natoms, 3), dtype=np.float64)

    for _ in range(ndS):
        mu, nu, dx, dy, dz = lines[i].split()
        mu = int(mu) - 1
        nu = int(nu) - 1
        at_mu = data['atom_map'][mu] - 1
        at_nu = data['atom_map'][nu] - 1

        vec = np.array([float(dx), float(dy), float(dz)], dtype=np.float64)

        data['dS_ext'][mu, nu, at_mu, :] = vec
        data['dS_ext'][nu, mu, at_mu, :] = vec
        data['dS_ext'][mu, nu, at_nu, :] = -vec
        data['dS_ext'][nu, mu, at_nu, :] = -vec
        i += 1

    # ---------- Kinetic derivative vectors ----------

    while i < len(lines) and 'derivatives of kinetic energy integrals' not in lines[i].lower():
        i += 1
    i += 1
    ndT = int(lines[i])
    i += 3

    data['dT_ext'] = np.zeros((nbasis, nbasis, natoms, 3), dtype=np.float64)

    for _ in range(ndT):
        mu, nu, dx, dy, dz = lines[i].split()
        mu = int(mu) - 1
        nu = int(nu) - 1
        at_mu = data['atom_map'][mu] - 1
        at_nu = data['atom_map'][nu] - 1

        vec = np.array([float(dx), float(dy), float(dz)], dtype=np.float64)

        data['dT_ext'][mu, nu, at_mu, :] = vec
        data['dT_ext'][nu, mu, at_mu, :] = vec
        data['dT_ext'][mu, nu, at_nu, :] = -vec
        data['dT_ext'][nu, mu, at_nu, :] = -vec
        i += 1

    # ---------- Nuclear-attraction derivative vectors ----------

    while i < len(lines) and 'derivatives of nucleus-electron energy integrals' not in lines[i].lower():
        i += 1
    i += 1
    ndV = int(lines[i])
    i += 3

    data['dV_ext'] = np.zeros((nbasis, nbasis, natoms, 3), dtype=np.float64)

    for _ in range(ndV):
        mu, nu, at, dx, dy, dz = lines[i].split()
        mu = int(mu) - 1
        nu = int(nu) - 1
        at = int(at) - 1
        vec = np.array([float(dx), float(dy), float(dz)], dtype=np.float64)

        data['dV_ext'][mu, nu, at, :] = vec
        data['dV_ext'][nu, mu, at, :] = vec
        i += 1

    # ---------- ERI derivative vectors ----------

    while i < len(lines) and 'derivatives of two-electron integrals' not in lines[i].lower():
        i += 1
    i += 1
    ndERI = int(lines[i])
    i += 2

    data['dERI_ext'] = np.zeros((nbasis, nbasis, nbasis, nbasis, natoms, 3), dtype=np.float64)

    for _ in range(ndERI):
        mu, nu, la, si, at, dx, dy, dz = lines[i].split()
        mu = int(mu) - 1
        nu = int(nu) - 1
        la = int(la) - 1
        si = int(si) - 1
        at = int(at) - 1
        vec = np.array([float(dx), float(dy), float(dz)], dtype=np.float64)

        # Fill all equivalent permutations related by ERI symmetry
        
        perms = [
            (mu, nu, la, si),
            (mu, nu, si, la),
            (nu, mu, la, si),
            (nu, mu, si, la),
            (la, si, mu, nu),
            (la, si, nu, mu),
            (si, la, mu, nu),
            (si, la, nu, mu),
        ]
        for p in perms:
            data['dERI_ext'][p[0], p[1], p[2], p[3], at, :] = vec

        i += 1

    return data