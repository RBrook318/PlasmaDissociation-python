import numpy as np

# Load the text file as a NumPy array
A = np.loadtxt('a.txt')
m = np.loadtxt('m.txt')
# geom = np.loadtxt('opt_geom.txt')
m = m * 1822.887  # Multiply m by 1836
geometry = """C   2.9929992543  -3.9946887565  -0.4491417187
O   0.6121754810  -2.7846264943  -0.8315142595
C   0.0506538168  -0.6636895674  0.4514345584
C   1.8776150005  1.5576166831  -0.1124543892
F   4.0994600695  1.1209082135  0.9776846052
F   2.2523968838  1.7788567366  -2.5790964978
F   1.0035042413  3.7331141883  0.7577522955
C   -2.7329998278  -0.0449445324  -0.1811134382
F   -4.1508599710  -2.0659541518  0.1987047447
F   -3.5955157717  1.7996000678  1.2735573933
F   -2.9572846560  0.6596801199  -2.5762348141
F   0.1595519692  -1.0108722785  3.0068571971
H   2.7659927041  -5.8646892817  -1.2396595215
H   4.4806048647  -2.9948181308  -1.4416256266
H   3.4251142338  -4.1272396563  1.5480248344
"""
T = 0.015834
n = 15
nmod = 3 * n - 6
Ax = A[:, 0]
Ay = A[:,1]
Az = A[:,2]

Ax = Ax.reshape(n, nmod, order = 'F')
Ay = Ay.reshape(n, nmod, order = 'F')
Az = Az.reshape(n, nmod, order = 'F')

rn = np.random.randn(nmod, 500)  # Use np.random.randn for standard normal distribution


# Initialize arrays for random
Meff = np.zeros(nmod)
rv = np.zeros((nmod, 500))


for i in range(nmod):
    for j in range(n):
        Meff[i] = Meff[i]+np.sum(((Ax[j, i]**2) + (Ay[j, i]**2) + (Az[j, i]**2)) * m[j])
    rv[i, :] = rn[i, :] * np.sqrt(2 * T / Meff[i])

# Saves effective masses for checking. 
# np.savetxt('Meff.txt', Meff, fmt='%f', delimiter='\t')

# Calculate the velocity by applying it through the tranformation matrix of normal modes.
Vx = np.dot(Ax, rv)
Vy = np.dot(Ay, rv)
Vz = np.dot(Az, rv)


Px = np.zeros((n,500))
Py = np.zeros((n,500))
Pz = np.zeros((n,500))
for i in range(n):
    Px[i,:] = Vx[i,:]*m[i]
    Py[i,:] = Vy[i,:]*m[i]
    Pz[i,:] = Vz[i,:]*m[i]
print(Px.shape)

# Save each momentum to separate files
for j in range(500):
    with open(f'Geometry.{j + 1}', 'w') as file:
        # Write the "momentum" line
        file.write(geometry)
        file.write("momentum\n")
        # Write Px, Py, and Pz for each atom on the same line
        for atom in range(n):
            # Access the Px, Py, and Pz values using the corresponding indices
            px_value = Px[atom, j]
            py_value = Py[atom, j]
            pz_value = Pz[atom, j]
            file.write(f'{px_value}  {py_value}  {pz_value}\n')
