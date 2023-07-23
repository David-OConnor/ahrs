# http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html

# https://matplotlib.org/stable/gallery/mplot3d/scatter3d.html#sphx-glr-gallery-mplot3d-scatter3d-py
import matplotlib.pyplot as plt
import numpy as np

points = [
   (0.3051851, -0.27564317, -0.5238197),
   (0.30909148, -0.2745445, -0.52699363),
   (0.36146122, -0.49586475, -0.02783288),
   (0.37647635, -0.46534625, -0.01147496),
   (0.34119692, -0.24683371, -0.5134434),
   (0.3150731, -0.2684408, -0.52833647),
   (0.3706168, -0.16016114, -0.0020752586),
   (0.35145116, -0.23218483, 0.03747673),
   (0.20032349, -0.3177587, -0.53566086),
   (0.20288706, -0.30494094, -0.5294351),
   (0.4471572, -0.21460617, -0.038575396),
   (0.3842891, -0.4954985, -0.35389262),
   (0.3075045, -0.4043092, -0.48365733),
   (0.26953948, -0.3737907, -0.50917083),
   (0.40101323, -0.44373912, -0.050294504),
   (0.37549976, -0.44312876, 0.0026856288),
   (0.3239845, -0.156621, -0.47181615),
   (0.37769708, -0.5294351, -0.28955963),
   (0.40174565, -0.529313, -0.22669148),
   (0.3399762, -0.044312876, -0.311655),
   (0.38380077, -0.049806207, -0.29224524),
   (0.37073886, -0.054200873, -0.16797388),
   (0.35804316, -0.509415, -0.052125614),
   (0.34632406, -0.52577287, -0.07238991),
   (0.32044435, -0.09424116, -0.41871396),
   (0.37537766, -0.105838194, -0.065309614),
   (0.4080935, -0.4034547, 0.003418073),
   (0.31336406, -0.440199, -0.46107364),
   (0.2989593, -0.18909268, -0.5055086),
   (0.039185766, -0.09778131, -0.26172674),
   (0.07751702, -0.059327982, -0.25049594),
   (0.20166631, -0.38575396, -0.5089267),
   (0.16443373, -0.08130131, -0.084109016),
   (0.11670278, -0.529313, -0.21179846),
   (0.4491104, -0.47328106, -0.11243019),
   (0.45338297, -0.4721824, -0.11035493),
   (0.3521836, -0.51564074, -0.11316264),
   (0.31299785, -0.054933317, -0.110965304),
   (0.3150731, -0.2684408, -0.52833647),
   (0.30909148, -0.2745445, -0.52699363),
   (0.20166631, -0.38575396, -0.5089267),
   (0.34815517, -0.40382093, -0.4730369),
   (0.39515367, -0.36805323, -0.4710837),
   (0.39649647, -0.3573107, -0.46632284),
   (0.24292734, -0.38550982, -0.5069735),
   (0.4326304, -0.50379956, -0.20435195),
   (0.30225533, -0.34717858, -0.50428784),
   (0.24744408, -0.34583575, -0.52577287),
   (0.18787195, -0.04504532, -0.2031312),
   (0.30225533, -0.319834, -0.532487),
   (0.37867367, -0.54628134, -0.30286568),
   (0.45594653, -0.49061555, -0.33472702),
   (0.39857173, -0.3668325, -0.003540147),
   (0.4080935, -0.4034547, 0.003418073),
   (0.4804834, -0.45765558, -0.1513718),
   (0.3856319, -0.37354657, -0.47401348),
   (0.41724905, -0.45216224, -0.06860561),
   (0.3075045, -0.4043092, -0.48365733),
   (0.36597797, -0.07373272, -0.09668264),
   (0.32044435, -0.09424116, -0.41871396),
   (0.35743278, -0.29456466, 0.045777764),
   (0.34632406, -0.52577287, -0.07238991),
   (0.28247932, -0.33582568, -0.53456223),
   (0.5215003, -0.29126865, -0.29493088),
   (0.38380077, -0.049806207, -0.29224524),
   (0.38868374, -0.46961883, -0.07483139),
   (0.2924894, -0.36414686, -0.5100253),
   (0.33179724, -0.35804316, -0.49940488),
   (0.13183996, -0.36414686, -0.50843835),
   (0.35145116, -0.23218483, 0.03747673),
   (0.18860438, -0.34644613, -0.52955717),
   (0.19678335, -0.33118686, -0.5252846),
   (0.44605854, -0.38026062, -0.4340953),
   (0.34119692, -0.24683371, -0.5134434),
   (0.27869502, -0.56776637, -0.3408307),
   (0.27503282, -0.32520524, -0.52919096),
   (0.07751702, -0.059327982, -0.25049594),
   (0.043092135, -0.18018128, -0.08899198),
]


x_offset = 0.123
y_offset = 0.071
z_offset = -0.180

for point in points:
   with_cal = np.array([point[0] - x_offset, point[1] - y_offset, point[2] - z_offset])
   mag = np.linalg.norm(with_cal)
   print(f"Mag: {mag}")


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


