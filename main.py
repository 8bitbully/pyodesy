from geodesy import ReferenceEllipsoid, GeodeticPoint


e = ReferenceEllipsoid("grs80")
# print(e)

# pt = GeodeticPoint(41.27737178, 39.21259856, 1000)
# pt.show().geog2geoc(e).show().geoc2geog(e).show().set_ellipsoid("hayford").geog2geoc(e).show().geoc2geog(e).show()
# pt2 = GeodeticPoint.empty().show()

# pt2 = GeodeticPoint(39.08885772, 32, 0.0, "reduced").geog2geoc(e).show()

# pt3 = GeodeticPoint(3794779.641, 3072951.963, 4089718.725, 'geocentric').geoc2geog(e).show()

# pt4 = GeodeticPoint(42.0).lat2eqdist(e)
# print(pt4)

# pt5 = GeodeticPoint(42.0).dist2eqlat(e, 4332953.916).show()

# pt6 = GeodeticPoint(26.0).longt2dist(e, [45, 36])
# print(pt6)

# pt7 = GeodeticPoint(39, 28).ellipsoidArea(e, GeodeticPoint(40, 29))
# print(pt7)

points = [GeodeticPoint(41.27737178, 39.21259856, 1000), 
GeodeticPoint(41.47737178, 39.31259856, 1000),
GeodeticPoint(41.07737178, 39.51259856, 1000)]

i = 0
for point in points:
    i += 1
    print(i)
    point.show().geog2geoc(e).show().geoc2geog(e).show()


