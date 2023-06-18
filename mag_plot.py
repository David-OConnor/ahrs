# This script plots magnetometer data, for visualizing offsets and
# calibration

# https://matplotlib.org/stable/gallery/mplot3d/scatter3d.html#sphx-glr-gallery-mplot3d-scatter3d-py
import matplotlib.pyplot as plt


fig = plt.figure()
ax = fig.add_subplot(projection='3d')


# C+P from println console from device as list of (x, y, z) tuples,
# eg (1., 0., 0.)

points = [
	
(-0.15, 0.05, 0.2),
(0.10232702, -0.011037018, -0.27938473),
(0.12576523, -0.1028367, -0.25545824),
(0.15835902, -0.15813623, -0.22323067),
(0.16434065, -0.18059786, -0.2077273),
(0.18033233, -0.23638569, -0.13936584),
(0.18838924, -0.23626362, -0.1239845),
(0.20914182, -0.22771843, -0.12813501),
(0.20096284, -0.15813623, -0.14595781),
(0.19754478, -0.14067966, -0.16414686),
(0.19266182, -0.13005921, -0.17440106),
(0.16165501, -0.10247047, -0.2179815),
(0.11136052, -0.0421659, -0.25496995),
(0.13284555, -0.095512256, -0.21065708),
(0.10794243, -0.17095402, -0.170983),
(0.13516495, -0.19756615, -0.12093265),
(0.15689412, -0.22503282, -0.038532674),
(0.17068848, -0.22307964, 0.062666714),
(0.17361826, -0.19195075, 0.13823053),
(0.18484908, -0.14446394, 0.18620564),
(0.18887752, -0.099540696, 0.21953185),
(0.20193943, -0.038381603, 0.226368),
(0.22110507, 0.034984894, 0.22356029),
(0.23245797, 0.09992828, 0.18937956),
(0.23819545, 0.170365, 0.14384595),
(0.24551988, 0.19355907, 0.10929899),
(0.25089115, 0.19636677, 0.09098788),
(0.25076905, 0.21016113, 0.06608479),
(0.2477172, 0.22542039, 0.033735156),
(0.2543092, 0.22126988, -0.031208232),
(0.25968048, 0.2013718, -0.08150272),
(0.26224402, 0.17402722, -0.12300791),
(0.26724905, 0.13264413, -0.15474714),
(0.2505249, 0.07356029, -0.20235603),
(0.23624226, -0.036428418, -0.20821558),
(0.23477736, -0.09477982, -0.18575396),
(0.23685262, -0.14495224, -0.1441267),
(0.22903988, -0.18438216, -0.09786065),
(0.22183752, -0.20708792, -0.058796957),
(0.21353647, -0.22710808, 0.0030945837),
(0.19376048, -0.22942747, 0.069014564),
(0.19022033, -0.2217368, 0.12345958),
(0.19070864, -0.2166097, 0.13053986),
(0.20474714, -0.22967161, 0.08842433),
(0.2120716, -0.23052613, 0.024579614),
(0.21146122, -0.22478865, -0.04292734),
(0.1918073, -0.21026187, -0.10823695),
(0.18252969, -0.21575518, -0.11617176),
(0.14444259, -0.22124852, -0.12569354),
(0.06094393, -0.23504288, -0.15059663),
(0.0028366894, -0.23138066, -0.15792109),
(-0.035006262, -0.22381206, -0.16085087),
(-0.0621067, -0.22625355, -0.14644612),
(-0.11411024, -0.23382215, -0.092245236),
(-0.18003023, -0.20244913, -0.035480812),
(-0.21860561, -0.16460617, -0.000201419),
(-0.23008057, -0.13726158, -0.009845272),
(-0.22690666, -0.13347729, 0.030927464),
(-0.2322779, -0.12175818, 0.056318864),
(-0.22605214, -0.10625477, 0.103317365),
(-0.23679465, -0.08245033, 0.09403974),
(-0.24326457, -0.04631642, 0.11674551),
(-0.21482132, -0.07573627, 0.1675283),
(-0.19016236, -0.0698767, 0.20720237),
(-0.17649007, -0.06743522, 0.23784295),
(-0.13925749, -0.063650936, 0.2583514),
(-0.12436446, -0.06877804, 0.27837154),
(-0.062594995, -0.08159582, 0.30571613),
(0.032500684, -0.110405296, 0.30693686),
(0.1423673, -0.108085886, 0.29326457),
(0.220983, -0.1221244, 0.23833126),
(0.2750618, -0.06719108, 0.21648),
(0.31107363, -0.04118931, 0.15800653),
(0.34159213, -0.0515656, 0.08915678),
(0.35123602, -0.04131138, -0.018268377),
(0.31534624, -0.02739494, -0.13606982),
(0.27115542, -0.020192575, -0.1950316),
(0.25845972, 0.021312602, -0.21114536),
(0.26578417, -0.006886501, -0.20455335),
(0.23196965, -0.116753146, -0.18428908),
(0.15957975, -0.2179525, -0.15145116),
(0.088410586, -0.25457472, -0.11104466),
(-0.019014567, -0.26043427, -0.09688406),
(-0.091648616, -0.23199104, -0.08601947),
(-0.17990814, -0.17779016, -0.07246925),
(-0.22568591, -0.1353084, -0.014728233),
(-0.25412917, -0.06499374, -0.021076083),
(-0.27012086, 0.016063418, -0.0043519437),
(-0.26609242, 0.09504532, 0.012860507),
(-0.2577914, 0.13044679, 0.032392353),
(-0.2403348, 0.17732322, 0.025312051),
(-0.19077274, 0.23872647, -0.00068971515),
(-0.11411024, 0.30843076, -0.023761705),
(-0.052340776, 0.32979372, -0.01765801),
(0.0130909085, 0.3411466, -0.016681418),
(0.08902097, 0.33333385, -0.02388379),
(0.1321131, 0.32320172, -0.01546067),
(0.17593768, 0.3083087, -0.012897119),
(0.23074892, 0.26973328, -0.011920527),
(0.27054507, 0.2371395, -0.01546067),
]

for x, y, z in points:
	ax.scatter(x, y, z)

ax.set_title("Magnetometer readings, during a calibration maneuver")
ax.set_xlabel('Mag X')
ax.set_ylabel('Mag Y')
ax.set_zlabel('Mag Z')

plt.show()