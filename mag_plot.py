# http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html

import numpy as np

points = [
(0.1853084, -0.11059908, -0.5229652),
(0.2106998, -0.12024293, -0.5164952),
(0.23963134, -0.12048708, -0.5105136),
(0.2684408, -0.13135166, -0.50648516),
(0.3103122, -0.32032228, -0.4528947),
(-0.07593005, 0.013061922, -0.07641835),
(-0.09277627, -0.010986663, -0.0794702),
(-0.02294992, -0.026490066, -0.040406507),
(-0.08703879, -0.26404613, -0.3178808),
(-0.074098945, -0.2689291, -0.32715842),
(0.29309976, 0.021118809, -0.4166387),
(0.44837794, -0.02600177, -0.2600177),
(0.08520768, -0.17615284, -0.5050203),
(0.36597797, -0.0019531846, -0.39747307),
(0.329722, 0.012451552, -0.40943632),
(0.36097294, -0.39674062, -0.29139072),
(0.3497421, -0.4015015, -0.27845088),
(0.48316905, -0.14062929, -0.2567217),
(0.4776757, -0.19788201, -0.26966155),
(0.07959227, -0.20410779, 0.086916715),
(0.2624592, 0.14636678, -0.2656331),
(-0.104373306, -0.03869747, -0.08630635),
(0.1617481, 0.17114781, -0.19019136),
(-0.08581805, -0.2318186, -0.34192938),
(-0.13183996, -0.22510453, -0.23340556),
(0.14404736, -0.2097232, 0.100833155),
(0.14270455, 0.16284677, -0.29651785),
(-0.12134159, -0.22888882, -0.10754723),
(-0.12317271, -0.21570483, -0.109500416),
(-0.12146367, -0.20801416, -0.1208533),
(0.023682363, 0.09533983, -0.39515367),
(0.021729179, 0.08850368, -0.40003663),
(0.08215583, -0.0987579, -0.50807214),
(0.38795128, -0.2633137, -0.42872402),
(0.36085087, -0.29822686, -0.43055514),
(0.046998505, 0.07068087, -0.026367992),
(0.2318186, 0.13013093, -0.30591753),
(0.34095278, 0.07470931, -0.3272805),
(0.47523424, -0.07763909, -0.25025177),
(0.21875668, 0.15063937, -0.30152288),
(-0.098391674, -0.2830897, -0.26087222),
(0.36463514, -0.19788201, -0.46229437),
(0.0075685903, 0.097048864, -0.3389996),
(0.13440351, 0.1424604, -0.36146122),
(-0.09131138, -0.27991578, -0.29016998),
(-0.11487167, -0.23841059, -0.12683493),
(-0.10046694, -0.2624592, -0.12573627),
(-0.0045167394, -0.36378065, -0.10095523),
(-0.0700705, -0.31556138, -0.113528855),
(0.1589404, 0.07470931, -0.052125614),
(0.07580798, 0.093264565, -0.4358043),
(-0.037598804, -0.11023286, -0.45240638),
(0.44898832, -0.2948088, -0.24854274),
(0.12439345, 0.12317271, -0.39796138),
(0.11926634, 0.117435224, -0.41163367),
(0.20410779, 0.13147373, -0.081789605),
(0.006225776, -0.19385357, 0.05468917),
(0.02356029, -0.20056765, 0.06982635),
(0.040650655, -0.28662986, 0.036378063),
(0.21570483, 0.11316264, -0.061891537),
(0.017090365, -0.17358929, 0.05346843),
(-0.046143986, 0.03503525, -0.05981628),
(-0.04382458, -0.03381451, -0.42164373),
(-0.04113895, 0.05774102, -0.36976227),
(0.019531846, -0.1861629, 0.052735984),
(0.25977355, -0.40504166, -0.35767692),
(0.23584704, -0.3558458, -0.4275033),
(0.29505295, -0.15747552, -0.50343335),
(0.32850122, -0.17456588, -0.48609883),
(0.2318186, -0.35193944, -0.4537492),
(0.08667257, 0.08362072, -0.44068727),
(-0.023194067, -0.18494217, 0.034302805),
(0.023438215, -0.096926786, -0.48915067),
(0.09350871, 0.08032472, -0.44837794),
(0.0926542, 0.07715079, -0.4504532),
(0.2379223, 0.04870754, 0.013550218),
(0.37501144, -0.21997742, -0.44605854),
(0.3932005, -0.23206274, -0.44007692),
(0.2981048, 0.055421613, -0.4528947),
(0.43140966, -0.3300882, -0.22510453),
(-0.10400708, -0.26673177, -0.18408765),
(-0.082888275, -0.20752586, -0.36756492),
(-0.103152566, -0.28601947, -0.22669148),
(-0.113528855, -0.27808467, -0.21448408),
(-0.09936827, -0.27588734, -0.19580676),
(-0.11316264, -0.23621327, -0.062623985),
(0.027466659, -0.36475724, -0.390759),
(0.4593646, -0.286874, -0.19629505),
(-0.10278634, -0.25366986, -0.15784173),
(0.15295877, -0.37342447, -0.423719),
(-0.053102206, -0.011108737, -0.4231086),
(0.43739128, -0.054200873, -0.36353648),
(0.43324077, -0.049073763, -0.36378065),
(0.4326304, -0.042725913, -0.37501144),
(0.4212775, 0.037354656, -0.23535874),
(0.16821803, 0.107181005, -0.39295632),
(0.051271096, 0.08594012, -0.374279),
(0.19641712, 0.07361065, -0.42677084),
(0.32947782, 0.073000275, -0.32801294),
(0.032593768, 0.10730308, -0.38966033),
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
   print('\nCenter at:',center)

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

print('\nSoft iron\n',M)



# This script plots magnetometer data, for visualizing offsets and
# calibration

# https://matplotlib.org/stable/gallery/mplot3d/scatter3d.html#sphx-glr-gallery-mplot3d-scatter3d-py
import matplotlib.pyplot as plt


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

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.set_title("Magnetometer readings, during a calibration maneuver.")
# ax.set_xlabel('Mag X')
# ax.set_ylabel('Mag Y')
# ax.set_zlabel('Mag Z')


# for x, y, z in points:
# 	# apply cal
# 	# pts_cal = M @ (np.array([x, y, z]) - center)
# 	pts_cal = (np.array([x, y, z]) - center)

# 	ax.scatter(pts_cal[0], pts_cal[1], pts_cal[2])

# plt.show()


