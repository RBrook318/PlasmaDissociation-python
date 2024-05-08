import numpy as np

def convert_to_bohr(coordinates):
    angstrom_to_bohr = 1.88973
    return [(atom, x * angstrom_to_bohr, y * angstrom_to_bohr, z * angstrom_to_bohr) for atom, x, y, z in coordinates]

def organise_modes(modes):
    numeric_modes = modes[:, 1:].astype(float)
    # Combine the columns
    column1 = np.concatenate((numeric_modes[:, 0], numeric_modes[:, 3], numeric_modes[:, 6]))
    column2 = np.concatenate((numeric_modes[:, 1], numeric_modes[:, 4], numeric_modes[:, 7]))
    column3 = np.concatenate((numeric_modes[:, 2], numeric_modes[:, 5], numeric_modes[:, 8]))
    # Stack the columns horizontally to create a new array
    reshaped_data = np.column_stack((column1, column2, column3))
    return reshaped_data

def create_geom(n,nmod,T,modes,m,mom_num):
    Ax = modes[:, 0]
    Ay = modes[:,1]
    Az = modes[:,2]
    Ax = Ax.reshape(n, nmod, order = 'F')
    Ay = Ay.reshape(n, nmod, order = 'F')
    Az = Az.reshape(n, nmod, order = 'F')
    rn = np.random.randn(nmod, mom_num)  # Use np.random.randn for standard normal distribution
    # Initialize arrays for random
    Meff = np.zeros(nmod)
    rv = np.zeros((nmod, mom_num))
    for i in range(nmod):
        for j in range(n):
            Meff[i] = Meff[i]+np.sum(((Ax[j, i]**2) + (Ay[j, i]**2) + (Az[j, i]**2)) * m[j])
        rv[i, :] = rn[i, :] * np.sqrt(2 * T / Meff[i])
    # Calculate the velocity by applying it through the tranformation matrix of normal modes.
    Vx = np.dot(Ax, rv)
    Vy = np.dot(Ay, rv)
    Vz = np.dot(Az, rv)
    Px = np.zeros((n,mom_num))
    Py = np.zeros((n,mom_num))
    Pz = np.zeros((n,mom_num))
    for i in range(n):
        Px[i,:] = Vx[i,:]*m[i]
        Py[i,:] = Vy[i,:]*m[i]
        Pz[i,:] = Vz[i,:]*m[i]
    
    return Px, Py, Pz
    
   
