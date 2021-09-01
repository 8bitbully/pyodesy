from math import sqrt
from geodesy.utils import *

class NivoEllipsoid:

    SecondEccentricity = 0
    Flattening = 0
    ThirdFlattening = 0
    _a = 0
    _b = 0
    _invf = 0
    _f = 0
    _ecc = 0
    _ecc2 = 0

    @property
    def SemimajorAxis(self):
        return self._a
    @property
    def SemiminorAxis(self):
        return self._b
    @property
    def InverseFlattening(self):
        return self._invf
    @property
    def Eccentricity(self):
        return self._ecc
    @property
    def MeanRadius(self):
        return (2*self._a + self._b) / 3
    
    @SemimajorAxis.setter
    def SemimajorAxis(self, a):
        self._a = a
        self._b = (1 - self.Flattening) * a

    @SemiminorAxis.setter
    def SemiminorAxis(self, b):
        a_ = self._a
        self._b = b
        self._ecc = sqrt(a_**2 - b**2) / a_
        self._ecc2 = sqrt(a_**2 - b**2) / b
        f = (a_- b) / a_
        self._invf = 1 / f
        self.Flattening = f
        self.ThirdFlattening = f / (2 - f)
        self.SecondEccentricity = self._ecc2

    @InverseFlattening.setter
    def InverseFlattening(self, invf):
        f = 1 / invf
        self._ecc  = sqrt((2 - f) * f)
        self._ecc2 = self._ecc / (1 - f)
        self._invf = invf
        self.Flattening = f
        self.ThirdFlattening = f / (2 - f)
        self._b = (1 - self.Flattening) * self._a
        self.SecondEccentricity = self._ecc2

    @Eccentricity.setter
    def Eccentricity(self, ecc):
         self._ecc = ecc
         e2 = ecc ** 2
         f = e2 / (1 + sqrt(1 - e2))
         self._b = (1 - f) * self.a
         self._ecc2 = e2 - (1 - e2)
         self._invf = 1 / f
         self.Flattening = f
         self.ThirdFlattening = f / (2 - f)
         self.SecondEccentricity = self._ecc2
    ### Methods:
    def geog2geoc(self, ellipsoid):

        B = self.x
        L = self.y
        h = self.z

        a_ = ellipsoid.SemimajorAxis
        b_ = ellipsoid.SemiminorAxis
        e2 = ellipsoid.Eccentricity ** 2
        eu2= ellipsoid.SecondEccentricity ** 2

        type = {
            "geographic": B,
            "reduced": atanf((a_/b_)*(tanf(B)/fi)) / fi,
            "geocentric": atanf((a_**2/b_**2)*(tanf(B)/fi)) / fi
        }

        B = type[self.type]

        n2 = eu2 * cosf(B) ** 2
        V = sqrt(1 + n2)
        W = sqrt(1 - e2*sinf(B)**2)

        xv = (a_**2 / (b_*V)+h) * cosf(B) * cosf(L)
        yv = (a_**2 / (b_*V)+h) * cosf(B) * sinf(L)
        zv = (b_ / V+h) * sinf(B)

        xw = (a_ / W+h) * cosf(B) * cosf(L)
        yw = (a_ / W+h) * cosf(B) * sinf(L)
        zw = (b_**2 / (a_*W)+h) * sinf(B)

        x = round((xv + xw) / 2, 3)
        y = round((yv + yw) / 2, 3)
        z = round((zv + zw) / 2, 3)

        return GeodeticPoint(x, y, z, 'geocentric')


    # pollard iteration
    def geoc2geog(self, ellipsoid):
        a_ = ellipsoid.SemimajorAxis
        b_ = ellipsoid.SemiminorAxis
        eu2= ellipsoid.SecondEccentricity**2

        x_, y_, z_ = self.x, self.y, self.z
        z__ = (b_*z_) / sqrt(x_**2 + y_**2 + z_**2)
        k_ = (a_/b_)**2

        while True:
            z_old = z__
            k = sqrt(x_**2+y_**2+(z_ + eu2*z__)**2)
            l_, m, n = x_/k, y_/k, (z_ + eu2*z__)/k
            r = l_**2 + m**2 + k_*n**2
            s = l_*x_ + m*y_ + k_*n*z_
            t = x_**2 + y_**2 + k_*z_**2 - a_**2
            h = (s + sqrt(s**2 - r*t)) / r
            root = sqrt(s**2 - r*t)

            if h < 0:
                z_new = z_ - n*h
            else:
                z_new = z_ - n*((s - root) /r)

            z__ = z_new
            if abs(z_new - z_old) < 1e-7: #0.0000001:
                break

        if h < 0:
            h = (s + root) / r
        else:
            h = (s - root) / r

        b = atanf(((z_ + eu2*z__) / sqrt(x_**2 + y_**2)) / fi) / fi
        l = atan2f(y_, x_) / fi

        return GeodeticPoint(round(b,5), round(l,5), round(h,3), 'geographic')

# def lat2eqdist(ellipsoid, B):
#     a_  = ellipsoid.SemimajorAxis
#     b_  = ellipsoid.SemiminorAxis
#     c   = a_**2 / b_
#     eu2 = ellipsoid.SecondEccentricity**2

#     A_ = c * (1 - 3/4*eu2 + 45/64*eu2**2 - 175/256*eu2**3 + 11025/16384*eu2**4 - 43659/65536*eu2**5) * fi
#     B_ = c * (-3/8*eu2 + 15/32*eu2**2 - 525/1024*eu2**3 + 2205/4096*eu2**4 - 72765/131072*eu2**5)
#     C_ = c * (15/256*eu2**2 - 105/1024*eu2**3 + 2205/16384*eu2**4 - 10395/65536*eu2**5)
#     D_ = c * (-35/3072*eu2**3 + 315/12288*eu2**4 - 31185/786432*eu2**5)
#     E_ = c * (315/131072*eu2**4 - 3465/524288*eu2**5 )
#     F_ = c * (-693/1310720*eu2**5)

#     G = A_*B + B_*sinf(2*B) + C_*sinf(4*B) + D_*sinf(6*B) + E_*sinf(8*B) + F_*sinf(10*B)

#     return round(G, 5)

# def dist2eqlat(ellipsoid, G):
#     a_  = ellipsoid.SemimajorAxis
#     b_  = ellipsoid.SemiminorAxis
#     c   = a_**2 / b_
#     eu2 = ellipsoid.SecondEccentricity**2

#     A_ = c * (1 - 3/4*eu2 + 45/64*eu2**2 - 175/256*eu2**3 + 11025/16384*eu2**4 - 43659/65536*eu2**5) * fi #... * invro
#     B__= (3/8*eu2 - 3/16*eu2**2 + 213/2048*eu2**3 - 255/4096*eu2**4 + 20861/524288*eu2**5) / fi
#     C__= (21/256*eu2**2 - 21/256*eu2**3 + 533/8192*eu2**4 - 197/4096*eu2**5) / fi
#     D__= (151/6144*eu2**3 - 453/12288*eu2**4 + 5019/131072*eu2**5) / fi
#     E__= (1097/131072*eu2**4 - 1097/65536*eu2**5) / fi
#     F__= (8011/2621440*eu2**5) / fi

#     sigma = G / A_

#     B = sigma + B__*sinf(2*sigma) + C__*sinf(4*sigma) + D__*sinf(6*sigma) + E__*sinf(8*sigma) + F__*sinf(10*sigma)

#     return round(B, 5)

# def longt2dist(ellipsoid, L1, L2, B):
#     a_  = ellipsoid.SemimajorAxis
#     b_  = ellipsoid.SemiminorAxis
#     c   = a_**2 / b_
#     e2  = ellipsoid.Eccentricity**2
#     eu2 = ellipsoid.SecondEccentricity**2

#     n2 = eu2 * cosf(B)**2
#     V  = sqrt(1 + n2)
#     W  = sqrt(1 - e2*sinf(B)**2)

#     N_ = c / V
#     N__= a_ / W
#     N  = (N_ + N__) / 2

#     l = L2 - L1

#     Sp = N * l * cosf(B) * fi
    
#     return round(Sp,3)

# def ellipsoidArea(ellipsoid, B1, B2, L1, L2):
#     b_  = ellipsoid.SemiminorAxis
#     e2  = ellipsoid.Eccentricity**2
            
#     Fa_ = 1 + 1/2*e2 + 3/8*e2**2 + 5/16*e2**3 + 35/128*e2**4
#     Fb_ = 1/6*e2 + 3/16*e2**2 + 3/16*e2**3 + 35/192*e2**4
#     Fc_ = 3/80*e2**2 + 1/16*e2**3 + 5/64*e2**4
#     Fd_ = 1/112*e2**3 + 5/256*e2**4
#     Fe_ = 5/2304*e2**4

#     dB = B2 - B1
#     dL = L2 - L1
#     Bm = (B1 + B2) / 2

#     # Z = 2*pi*b^2 INTEGRAL(B1 -> B2) cos(B) / (1 - e^2*sin^2*B)^2
#     # --> dB
#     Z_ = 4 * pi * b_**2 * (Fa_*sinf(dB/2)*cosf(Bm) - 
#                            Fb_*sinf(3*dB/2)*cosf(3*Bm) + 
#                            Fc_*sinf(5*dB/2)*cosf(5*Bm) - 
#                            Fd_*sinf(7*dB/2)*cosf(7*Bm) + 
#                            Fe_*sinf(9*dB/2)*cosf(9*Bm)) # m^2

#     Z = Z_ / 1e+6 #km**2
#     F = dL * Z / (2*pi) * fi

#     return round(F, 3)

class ReferenceEllipsoid(NivoEllipsoid):

    def __init__(self, ellipsoid: str):
        self.ELLIPSOID = self._ellipsoid(ellipsoid.upper())
        self.Name = self.ELLIPSOID["Name"]
        self.SemimajorAxis = self.ELLIPSOID["SemimajorAxis"]
        self.SemiminorAxis = self.ELLIPSOID["SemiminorAxis"]
        self.InverseFlattening = self.ELLIPSOID["InverseFlattening"]

    def _ellipsoid(self, ell):
        ellipsoids = {
            "GRS80":{
                "Name":"Geodetic Reference System 1980",
                "SemiminorAxis":6356752.314140347,
                "SemimajorAxis":6378137,
                "InverseFlattening":298.2572221008827 
            },
            "HAYFORD":{
                "Name":"Hayford 1909",
                "SemiminorAxis":6356911.94613,
                "SemimajorAxis":6378388,
                "InverseFlattening":297 
            }
        }
        return ellipsoids[ell]

    def __str__(self):
        result = f"Name: {self.Name} \nSemimajorAxis: {self.SemimajorAxis} \nSemiminorAxis: {self.SemiminorAxis}\nEccentricity: {self.Eccentricity} "
        return result

class GeodeticPoint(ReferenceEllipsoid):

    def __init__(self, x: float, y: float, z: float = None, type: str = 'geographic') -> None:
        self._x = x
        self._y = y
        self._z = z
        self._type = type

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    @property
    def type(self):
        return self._type

    @staticmethod
    def empty():
        return GeodeticPoint(None, None, None, None)

    @x.setter
    def x(self, c_x):
        self._x = c_x
    
    @y.setter
    def y(self, c_y):
        self._y = c_y

    @z.setter
    def z(self, c_z):
        self._z = c_z
    
    @type.setter
    def type(self, c_type):
        self._type = c_type

    def show(self):
        print(f"({self._x}, {self._y}, {self._z}), type= [{self._type}]")
        return GeodeticPoint(self._x, self._y, self._z, self._type)

    def set_ellipsoid(self, name: str):
        self.ELLIPSOID = self._ellipsoid(name.upper())
        self.Name = self.ELLIPSOID["Name"]
        self.SemimajorAxis = self.ELLIPSOID["SemimajorAxis"]
        self.SemiminorAxis = self.ELLIPSOID["SemiminorAxis"]
        self.InverseFlattening = self.ELLIPSOID["InverseFlattening"]
        return GeodeticPoint(self._x, self._y, self._z, self._type)
