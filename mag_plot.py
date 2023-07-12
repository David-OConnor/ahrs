# http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html

# https://matplotlib.org/stable/gallery/mplot3d/scatter3d.html#sphx-glr-gallery-mplot3d-scatter3d-py
import matplotlib.pyplot as plt
import numpy as np

points = [
(0.10657064, -0.126957, -0.5290689),
(0.19409773, -0.0040284432, -0.59315777),
(0.074221015, 0.14331493, 0.21411787),
(0.19409773, 0.027466659, 0.21997742),
(0.4947661, -0.07910398, -0.037354656),
(0.53028965, 0.1640675, -0.24610126),
(0.5594653, 0.036378063, -0.24011964),
(0.15881832, 0.17981505, 0.21265298),
(0.16357921, 0.015991699, -0.58156073),
(0.27369, 0.10681479, -0.5726493),
(0.26062807, 0.47462386, -0.16052736),
(-0.09253212, -0.11328471, 0.12073123),
(0.45582446, 0.23816645, -0.022583697),
(0.560564, 0.0253914, -0.17285684),
(0.50648516, 0.16028321, -0.31934568),
(0.045533616, 0.0888699, 0.21497238),
(0.4001587, 0.08911405, -0.5240638),
(0.4268929, 0.15601063, -0.46253854),
(0.24097416, 0.33533737, 0.09753716),
(0.37818536, 0.016479995, 0.15710929),
(-0.07702872, 0.15930662, -0.51979125),
(-0.048219245, 0.20325327, -0.5222327),
(0.50160223, 0.20459609, -0.08471938),
(0.40992463, 0.3348491, -0.06884976),
(-0.02844325, -0.08264413, -0.56215096),
(0.059694204, -0.2106998, -0.44508195),
(0.3170263, 0.016235847, 0.1739555),
(0.47340313, 0.31629384, -0.21594897),
(0.06372265, 0.34119692, -0.49903867),
(-0.13782158, -0.033936583, 0.14380322),
(0.3503525, 0.13830988, -0.5317545),
(0.3526719, 0.09814753, -0.55384994),
(0.33436078, 0.13330485, -0.53199863),
(0.28321177, 0.10937834, 0.18298899),
(0.16357921, 0.015991699, -0.58156073),
(0.19409773, -0.0040284432, -0.59315777),
(0.11890011, 0.1297647, -0.58900726),
(0.10364086, 0.16113773, -0.5783868),
(0.074221015, 0.14331493, 0.21411787),
(0.045533616, 0.0888699, 0.21497238),
(0.15881832, 0.17981505, 0.21265298),
(0.06762902, 0.3685415, 0.10364086),
(0.10779138, 0.056764428, -0.57850885),
(0.40992463, 0.3348491, -0.06884976),
(0.10657064, -0.126957, -0.5290689),
(0.27369, 0.10681479, -0.5726493),
(0.2542802, -0.16211432, 0.12671286),
(0.059694204, -0.2106998, -0.44508195),
(0.10180975, 0.08972442, -0.57704395),
(0.19409773, 0.027466659, 0.21997742),
(0.27405623, -0.011108737, -0.5710623),
(0.3526719, 0.09814753, -0.55384994),
(0.39332256, 0.40382093, -0.18689536),
(-0.048219245, 0.20325327, -0.5222327),
(0.51564074, 0.036133915, -0.39283425),
(0.46961883, 0.14050722, -0.4123661),
(0.06433302, 0.41871396, 0.044068728),
(0.13623463, 0.23999757, -0.54249704),
(0.15454574, 0.103274636, -0.57545704),
(0.09509568, 0.036744285, -0.57008576),
(0.45619068, 0.19580676, 0.029175695),
(0.13013093, 0.08337657, -0.5746025),
(-0.10095523, 0.3319193, 0.082888275),
(0.4890286, 0.009643849, 0.0019531846),
(0.47242653, 0.05285806, 0.045167394),
(0.45582446, 0.23816645, -0.022583697),
(0.52735984, 0.0380871, -0.015259255),
(0.13415937, 0.10913419, -0.5832698),
]


def ls_ellipsoid(xx,yy,zz):

   # change xx from vector of length N to Nx1 matrix so we can use hstack
   x = xx[:,np.newaxis]
   y = yy[:,np.newaxis]
   z = zz[:,np.newaxis]

   #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz = 1
   J = np.hstack((x*x,y*y,z*z,x*y,x*z,y*z, x, y, z))
   K = np.ones_like(x) #column of ones

   # print(J.shape, "J SHAPE")

   # print("\n\nJ: ", J);

   # print("K", K)

   #np.hstack performs a loop over all samples and creates
   #a row in J for each x,y,z sample:
   # J[ix,0] = x[ix]*x[ix]
   # J[ix,1] = y[ix]*y[ix]
   # etc.

   JT=J.transpose()
   JTJ = np.dot(JT,J)
   InvJTJ=np.linalg.inv(JTJ);
   ABC= np.dot(InvJTJ, np.dot(JT,K))

   # Rearrange, move the 1 to the other side
#  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz - 1 = 0
#    or
#  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz + J = 0
#  where J = -1
   eansa=np.append(ABC,-1)

   return eansa


def polyToParams3D(vec):

   # convert the polynomial form of the 3D-ellipsoid to parameters
   # center, axes, and transformation matrix
   # vec is the vector whose elements are the polynomial
   # coefficients A..J
   # returns (center, axes, rotation matrix)

   #Algebraic form: X.T * Amat * X --> polynomial form

   print('\npolynomial\n',vec)

   Amat=np.array(
   [
   [ vec[0],     vec[3]/2.0, vec[4]/2.0, vec[6]/2.0 ],
   [ vec[3]/2.0, vec[1],     vec[5]/2.0, vec[7]/2.0 ],
   [ vec[4]/2.0, vec[5]/2.0, vec[2],     vec[8]/2.0 ],
   [ vec[6]/2.0, vec[7]/2.0, vec[8]/2.0, vec[9]     ]
   ])

   print('\nAlgebraic form of polynomial\n',Amat)

   #See B.Bartoni, Preprint SMU-HEP-10-14 Multi-dimensional Ellipsoidal Fitting
   # equation 20 for the following method for finding the center
   A3=Amat[0:3,0:3]
   A3inv=np.linalg.inv(A3)
   ofs=vec[6:9]/2.0
   center=-np.dot(A3inv,ofs)

   # Center the ellipsoid at the origin
   Tofs=np.eye(4)
   Tofs[3,0:3]=center
   R = np.dot(Tofs,np.dot(Amat,Tofs.T))
   print ('\nAlgebraic form translated to center\n',R,'\n')

   R3=R[0:3,0:3]
   R3test=R3/R3[0,0]
   print('normed \n',R3test)

   s1=-R[3, 3]
   R3S=R3/s1
   (el,ec)=np.linalg.eig(R3S)

   recip=1.0/np.abs(el)
   axes=np.sqrt(recip)

   print('\nAxes are\n',axes  ,'\n')

   inve=np.linalg.inv(ec) #inverse is actually the transpose here

   return (center,axes,inve)


# Rearrange, move the 1 to the other side
#  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz - 1 = 0
#    or
#  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz + J = 0
#  where J = -1
   eansa=np.append(ABC,-1)

   return (eansa)


xx = np.array([pt[0] for pt in points])
yy = np.array([pt[1] for pt in points])
zz = np.array([pt[2] for pt in points])

poly_pts = ls_ellipsoid(xx, yy, zz)

(center, axes, R) = polyToParams3D(poly_pts)

L = np.diag([1/axes[0],1/axes[1],1/axes[2]])
M=np.dot(R.T,np.dot(L,R))

print('\nHard iron\n',center)

print('\nSoft iron\n',M)



# This script plots magnetometer data, for visualizing offsets and
# calibration


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_title("Magnetometer readings, during a calibration maneuver.")
ax.set_xlabel('Mag X')
ax.set_ylabel('Mag Y')
ax.set_zlabel('Mag Z')

# C+P from println console from device as list of (x, y, z) tuples,
# eg (1., 0., 0.)

# ax.set_aspect('equal')

for x, y, z in points:
	ax.scatter(x, y, z)

plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_title("Magnetometer readings, with calibration applied.")
ax.set_xlabel('Mag X')
ax.set_ylabel('Mag Y')
ax.set_zlabel('Mag Z')


for x, y, z in points:
	# apply cal
	pts_cal = M @ (np.array([x, y, z]) - center)

	ax.scatter(pts_cal[0], pts_cal[1], pts_cal[2])

plt.show()


