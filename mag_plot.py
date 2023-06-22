# temp hijacking this
import numpy as np

points = [
(0.41688284, 0.1495407, -0.10693686),
(0.16418958, 0.22327341, -0.028931547),
(0.15002899, 0.17871639, 0.039796136),
(0.1439253, 0.12854396, 0.08252205),
(0.47181615, 0.0028077029, -0.006714072),
(0.16895047, -0.06335643, 0.1782281),
(0.2238838, -0.22205268, 0.18347728),
(0.26868495, -0.22681357, 0.17334513),
(0.35670033, -0.3977172, 0.07324442),
(0.37647635, -0.33936584, 0.10242012),
(0.29078037, -0.48133793, 0.03723258),
(0.50294507, 0.10144353, -0.23438215),
(0.3361919, -0.43433943, 0.03112888),
(0.3869747, -0.3894162, 0.048463393),
(0.41334268, -0.3173925, 0.070192575),
(0.16382337, 0.017456587, 0.15759759),
(0.40174565, 0.13830988, -0.0695822),
(0.05566576, -0.15051728, 0.1861629),
(0.4364147, 0.114627525, -0.08459731),
(0.15491195, 0.07470931, 0.120975375),
(0.13257241, 0.015869625, -0.60963774),
(0.14587848, -0.0987579, -0.63039035),
(0.1495407, -0.105838194, -0.63734853),
(0.1767632, -0.09729301, -0.6399121),
(-0.034546953, -0.22254097, -0.53028965),
(0.143437, 0.26880702, -0.14258248),
(0.31690422, 0.22168645, -0.123783074),
(0.26966155, 0.23364972, -0.11035493),
(0.15295877, 0.26319164, -0.136845),
(0.3098239, 0.21729179, -0.14539018),
(0.37623218, 0.16663106, -0.045289468),
(0.5205237, 0.017334513, -0.06946013),
(0.3846553, -0.22791223, 0.1208533),
(0.3704947, -0.24243905, 0.12683493),
(0.18543047, -0.17224647, 0.19019136),
(0.3864864, -0.25464645, 0.110965304),
(0.31299785, -0.24036378, -0.5953551),
(0.30896938, -0.227546, -0.6023133),
(0.46070743, -0.04895169, 0.04174932),
(0.2200995, 0.24414808, -0.09094516),
(0.33545947, -0.21045564, 0.159917),
(0.39234596, -0.21472824, 0.12475967),
(0.41444135, -0.18628499, 0.1076693),
(0.43836787, 0.07373272, -0.4395886),
(0.48890653, 0.1259804, -0.24732201),
(-0.14429152, -0.005981628, -0.3822138),
(-0.14172795, 0.007202368, -0.3960082),
(-0.092043824, 0.03320414, -0.4547258),
(-0.06640828, 0.034913175, 0.08130131),
(0.20532854, -0.41029084, 0.10180975),
(0.18066958, -0.42799157, 0.08545183),
(0.31153294, 0.23194067, -0.19202246),
(0.3949095, 0.16589862, -0.13061923),
(0.0926542, -0.4626606, 0.061769463),
(-0.022827845, -0.002929777, 0.12671286),
(0.09399701, -0.49220252, 0.03112888),
(0.0860622, -0.46620074, 0.05896176),
(0.074221015, -0.45899838, 0.06225776),
(0.22498246, -0.49806207, 0.033082064),
(0.2573321, -0.4865871, 0.044923246),
(0.46436965, 0.058717612, -0.04358043),
(0.37464523, 0.19629505, -0.20899075),
(0.37403485, 0.19348735, -0.18347728),
(-0.19385357, -0.3460799, -0.10058901),
(-0.11560412, -0.13183996, -0.48890653),
(-0.11816767, -0.16370128, -0.48255867),
(-0.18030335, -0.24597919, -0.32520524),
(-0.16443373, -0.23487045, -0.35975218),
(0.38074893, -0.13098544, 0.12781152),
(0.30970183, 0.08508561, 0.09155553),
(0.38929412, -0.17627491, 0.12146367),
(0.374279, -0.16919462, 0.14526811),
(-0.16870633, -0.24268319, -0.3573107),
(-0.1853084, -0.2351146, -0.31189916),
(0.20252083, -0.3582873, 0.14563432),
(-0.18091373, -0.29126865, -0.3706168),
(-0.1866512, -0.27771842, -0.35901976),
(0.58302563, -0.22559282, -0.21228676),
(0.24756615, -0.08545183, 0.17004913),
(0.41810358, -0.15918455, 0.117435224),
(0.26319164, -0.08142339, 0.17688528),
(0.44691306, -0.20874661, 0.08374279),
(0.37965026, 0.18762779, -0.12671286),
(0.2690512, 0.28760645, -0.25281534),
(0.16150396, 0.17578661, -0.5271157),
(0.5416425, 0.04358043, -0.26966155),
(-0.19226661, -0.0380871, -0.1716361),
(0.41444135, 0.07214576, 0.006591998),
(0.54713583, -0.06201361, -0.119876705),
(0.38526568, 0.070192575, 0.052125614),
(0.36548966, 0.042237617, 0.08484146),
(0.43751335, 0.08777124, -0.032105472),
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


# C+P from println console from device as list of (x, y, z) tuples,
# eg (1., 0., 0.)

ax.set_aspect('equal')

for x, y, z in points:

	# apply cal
	pts_cal = M @ (np.array([x, y, z]) - center)

	ax.scatter(pts_cal[0], pts_cal[1], pts_cal[2])

ax.set_title("Magnetometer readings, during a calibration maneuver.")
ax.set_xlabel('Mag X')
ax.set_ylabel('Mag Y')
ax.set_zlabel('Mag Z')

plt.show()


