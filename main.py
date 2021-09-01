from geodesy import ReferenceEllipsoid, GeodeticPoint


e = ReferenceEllipsoid("grs80")
pt = GeodeticPoint(41.27737178, 39.21259856, 1000)
pt.show().geog2geoc(e).show().geoc2geog(e).show().set_ellipsoid("hayford").geog2geoc(e).show().geoc2geog(e).show()

pt2 = GeodeticPoint.empty().show()

# print(e)

# pt = geog2geoc(e, 41.27737178, 39.21259856, 1000)
# print(f"X: {pt[0]}, Y: {pt[1]}, Z: {pt[2]}")

# pt2 = geog2geoc(e, 39.08885772, 32, 0, "reduced")
# print(f"X: {pt2[0]}, Y: {pt2[1]}, Z: {pt2[2]}")

# pt3 = geoc2geog(e, 3794779.641, 3072951.963, 4089718.725)
# print(f"B: {pt3[0]}, L: {pt3[1]}, h: {pt3[2]}")

# g = lat2eqdist(e, 42)
# print(f"G: {g}")

# b = dist2eqlat(e, 4332953.916)
# print(f"B: {b}")

# sp = longt2dist(e, 26, 45, 36)
# print(f"Sp: {sp}")

# f = ellipsoidArea(e, 39, 40, 28, 29)
# print(f"F: {f}")
# # fi = pi / 180
# # sinf = lambda rad: sin(rad  *  fi)
# # print(sinf(41.27737178))
# # print(sin(41.27737178))


