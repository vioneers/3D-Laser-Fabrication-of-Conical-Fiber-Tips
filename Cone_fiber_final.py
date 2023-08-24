# Importing the libraries
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # Import Axes3D explicitly  from numpy import pi
from math import pi, sin, cos, sqrt
# ===============================

# Opening the output files
big_base = open("big_base.gwl", "w")
exterior = open("exterior.gwl", "w")
interior = open("interior.gwl", "w")

# ===============================

Rmax = 30  # large radius
Rmin = 0.1 # small radius

# ===============================

# Center
xc = 50
yc = 50
zc = 0

# Starting point
x0 = xc
y0 = yc - Rmax # start from the margin not from the center
z0 = zc

# ===============================

dS = 0.5 # maximum increment for point-to-point path

# Vectors initialization
x = []
y = []
z = []

# ===============================

# Big base of the cone
R_max_spiral = Rmax # max spiral radius
R_min_spiral = Rmin # min spiral radius
deltaR = 1          # radius increment between two consecutive spirals
N_spiral_loops = int((R_max_spiral - R_min_spiral) / deltaR)

for i in range (0, N_spiral_loops):
    Ri = R_min_spiral + deltaR * i # initial radius of spiral i
    N_seg = int(2. * pi * Ri / dS) # number of segments in the spiral loop
    if N_seg < 6:                  # too few when reaching the center
        N_seg = 12
    dR_phi = deltaR / N_seg        # gradual increment between segments of the same spiral loop
    phi_steps = np.linspace (0, 2. * pi, int(N_seg), endpoint = True) # angular steps array for constant distance increment
    
    zi = 0.
    for phi in phi_steps:
        xi = xc + Ri * sin(phi)
        yi = yc + Ri * cos(phi)
        
        x.append(xi)
        y.append(yi)
        z.append(zi)
        
        string = '{0:5.2f}\t{1:5.2f}\t{2:5.2f}\n'
        Ri = Ri + dR_phi           # gradual increment between segments of the same spiral loop
        big_base.write(string.format(xi,yi,zi))

big_base.write('write')            # Write function at the end of each file for the Nanoscribe
big_base.close()

#Plotting the big base
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.set_title('Cone Base')
ax1.set_xlabel('X-axis')
ax1.set_ylabel('Y-axis')
ax1.plot(x, y, 'b')
ax1.set_aspect('equal', 'box')
# plt.show()

# ===============================

# Exterior of the cone
x = [] # reinitialize the vectors
y = []
z = []
h = 100                 # cone height
dZ = 0.5                # step size in Z direction
Nspirals_Z = int (h/dZ) # number of spirals in Z direction
deltaR_Z = (Rmax - Rmin) / Nspirals_Z

for i in range (0, Nspirals_Z):
    Ri = Rmax - deltaR_Z * i # initial radius of spiral i
    N_seg = int(2. * pi * Ri / dS) # number of segments in the spiral loop
    dR_phi = deltaR / N_seg        # gradual increment between segments of the same spiral loop
    phi_steps = np.linspace (0, 2. * pi, int(N_seg), endpoint = True) # angular steps array for constant distance increment
    
    zi = dZ * i         # z increases by dZ every spiral
    dz_phi = dZ / N_seg # z increase evenly distributed from segment to segment
    
    j = 0
    for phi in phi_steps:
        xi = xc + Ri * sin(phi)
        yi = yc + Ri * cos(phi)
        
        x.append(xi)
        y.append(yi)
        z.append(zi)
        
        string = '{0:5.2f}\t{1:5.2f}\t{2:5.2f}\n'
        Ri = Ri - dR_phi # radius decreases gradually 
        zi = zi + dz_phi # z increases gradually
        exterior.write(string.format(xi,yi,zi))
        j = j + 1 

exterior.write('write')  # Write function at the end of each file for the Nanoscribe
exterior.close() 

#Plotting the exterior
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111, projection='3d')
ax2.set_title('Cone Exterior')
ax2.set_xlabel('X-axis')
ax2.set_ylabel('Y-axis')
ax2.set_zlabel('Z-axis')
ax2.plot(x, y, z, 'b')
ax2.view_init(elev=45., azim=45) #elev -> XY, azim -> Z (see code documentation)
# plt.show()

# ===============================

# Interior of the cone
x = [] # reinitialize the vectors
y = []
z = []
h = 100                 # cone height
dZ = 0.5                # step size in Z direction -> 0.5 to test use 20
Nspirals_Z = int (h/dZ) # number of spirals in Z direction
deltaR_Z = (Rmax - Rmin) / Nspirals_Z
zi = 0.                 # initial z

for i in range (0, Nspirals_Z): #to test use 1 instead of Nspirals_Z 
    Ri = Rmax - deltaR_Z * i # initial radius of spiral i
    
    # Spiral
    R_max_spiral = Ri  # max spiral radius
    R_min_spiral = 0.1 # min spiral radius
    deltaR = 1         # spiral radius increment -> 1 to test use 5
    N_spiral_loops = int((R_max_spiral - R_min_spiral) / deltaR)
    
    for j in range (0, N_spiral_loops):
        Ri = R_max_spiral - deltaR * j     # initial radius of spiral i
        N_seg = int(2. * pi * Ri / dS)     # number of segments in the spiral loop
        if N_seg < 6:                      # not to be too fragmented towards the center when the radius is very small
            N_seg = 12
    
        dR_phi = deltaR / N_seg        # gradual increment between segments of the same spiral loop
        phi_steps = np.linspace (0, 2. * pi, int(N_seg), endpoint = True) # angular steps array for constant distance increment    
    
        for phi in phi_steps:
            xi = xc + Ri * sin(phi)
            yi = yc + Ri * cos(phi)
            
            x.append(xi)
            y.append(yi)
            z.append(zi)
            
            string = '{0:5.2f}\t{1:5.2f}\t{2:5.2f}\n'
            Ri = Ri - dR_phi # radius decreases gradually 
            interior.write(string.format(xi,yi,zi))
            
    zi = zi + dZ/2
    
    #Ellipse
    alpha = 0
    N_alpha = 12 # number of ellipsess overlapping in the center of the "circle"
    # RMi = Rmax - 2 * deltaR_Z * i
    # Rmi = RMi / 3
    # Ri = Rmi * 5
    
    deltaR = 1/2 * deltaR_Z # diff between the radii of two consecutive spiral levels -> HALF because the other half has already been covered with the spiral
    dRMe = deltaR / 12    # diff between the radii of two consecutive ellipses on the same level
    
    #Ri = R_max_spiral - 1 / 2 * deltaR_Z # in the middle between two consecutive spiral levels 1/2 ????? 3/4
    RMi = R_max_spiral / 2  # because the ellipse begins in the middle and the spiral on the margin
    Rmi = RMi / 3
    Ri = 5 * Rmi 
    #RMi = (Ri + Rmi) / 2
    
    for i_alpha in range(0,N_alpha):
        alpha = alpha + 2 * pi / N_alpha # gradual increase in angle
        dtheta = dS / Ri
        N_seg = int(2. * pi * Ri / dS)
        
        for j in range (0, N_seg):
            xi_0 = RMi - Rmi + RMi * sin (j * dtheta)
            yi_0 = Rmi * cos (j * dtheta)
            
            xi = xc + xi_0 * cos(alpha) + yi_0 * sin(alpha)
            yi = yc - xi_0 * sin(alpha) + yi_0 * cos(alpha)
            
            x.append(xi)
            y.append(yi)
            z.append(zi)
            
            string = '{0:5.2f}\t{1:5.2f}\t{2:5.2f}\n'
            interior.write(string.format(xi,yi,zi))
            
        RMi = RMi - dRMe
        Rmi = RMi / 3
        Ri = 5 * Rmi 
            
    zi = zi + dZ/2
    
interior.write('write')  # Write function at the end of each file for the Nanoscribe
interior.close() 

#Plotting the interior of the cone
fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111,  projection='3d')
ax3.set_xlabel('X-axis')
ax3.set_ylabel('Y-axis') 
ax3.set_zlabel('Z-axis')
ax3.set_title('Cone Interior')
# ax3.scatter(x, y, z, color='grey')
ax3.plot(x, y, z, 'b')
ax3.view_init(elev=45., azim=45)
# plt.show()

# Plotting all figures simultaneously 
plt.show()