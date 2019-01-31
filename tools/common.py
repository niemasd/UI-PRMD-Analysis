#!/usr/bin/env python3
from gzip import open as gopen
from math import cos,degrees,inf,radians,sin,sqrt
from numpy import mat
from os.path import isfile
from statistics import mean
INF = inf
ANGLE_UNITS = ['degrees','radians']

# convert an Eulerian angle to quaternion
def eulerian_to_quaternion(yaw_z, pitch_y, roll_x, unit='degrees'):
    if unit not in ANGLE_UNITS:
        raise ValueError("Invalid unit: %s (valid options: %s)" % (unit, ','.join(ANGLE_UNITS)))
    if unit == 'degrees':
        yaw_z = radians(yaw_z)
        pitch_y = radians(pitch_y)
        roll_x = radians(roll_x)
    cy = cos(yaw_z * 0.5)
    sy = sin(yaw_z * 0.5)
    cp = cos(pitch_y * 0.5)
    sp = sin(pitch_y * 0.5)
    cr = cos(roll_x * 0.5)
    sr = sin(roll_x * 0.5)
    q_w = cy * cp * cr + sy * sp * sr
    q_x = cy * cp * sr - sy * sp * cr
    q_y = sy * cp * sr + cy * sp * cr
    q_z = sy * cp * cr - cy * sp * sr
    return mat([[q_w],[q_x],[q_y],[q_z]])

# distance functions
def d_squared_euclidean(a,b): # squared Euclidean distance
    assert len(a) == len(b), "Can't compute Euclidean distance on different-dimensional points"
    return sum((a[i]-b[i])**2 for i in range(len(a)))

def d_euclidean(a,b): # Euclidean distance
    return sqrt(d_squared_euclidean(a,b))

DISTANCE = {
    'angles': {
        'euclidean': d_euclidean,
        'squared_euclidean': d_squared_euclidean,
    },
    'positions': {
        'euclidean': d_euclidean,
        'squared_euclidean': d_squared_euclidean,
    },
}

def distance(a,b,f):
    assert len(a) == len(b), "Can't compute distance on different number of points"
    if isinstance(a[0],tuple) or isinstance(a[0],list):
        if not isinstance(b[0],tuple) and not isinstance(b[0],list):
            raise TypeError("Can't compute distance between single point and multiple points")
        return sum(f(a[i],b[i]) for i in range(len(a))) # just sum individual distances (maybe expand in the future)
    else:
        return f(a,b)

# class to load time-series data
class TimeSeries:
    # constructor
    def __init__(self, data, data_type):
        if not isinstance(data_type,str):
            raise TypeError("Type must be a string")
        data_type = data_type.lower()
        if data_type not in DISTANCE:
            raise ValueError("Invalid type: %s (valid options: %s)" % (data_type, ','.join(sorted(DISTANCE.keys()))))
        self.type = data_type
        if isinstance(data, list):
            self.values = data
        else:
            if isfile(data):
                if data.lower().endswith('.gz'):
                    data = gopen(data).read().decode()
                else:
                    data = open(data).read()
            self.values = list()
            for l in data.strip().splitlines():
                self.values.append(list())
                cols = [float(e) for e in l.strip().split(',')]
                for i in range(0,len(cols),3):
                    self.values[-1].append(tuple(cols[i:i+3]))
            if len(self.values) == 0:
                raise ValueError("No time points!")

    # getters
    def num_rows(self): # each row is a time point
        return len(self.values)
    def num_cols(self): # each column is a 3D point/angle
        return len(self.values[0])
    def get_centers(self): # get the centers of this TimeSeries object
        return [tuple(mean(self.values[i][j][v] for i in range(len(self.values))) for v in range(3)) for j in range(len(self.values[0]))]
    def __len__(self):
        return self.num_rows() # length = number of rows

    # check equality
    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        if other.type != self.type:
            return False
        if other.num_rows() != self.num_rows():
            return False
        if other.num_cols() != self.num_cols():
            return False
        for i in range(self.num_rows()):
            if self.values[i] != other.values[i]:
                return False
        return True
    def __ne__(self, other):
        return not self == other

    # return a copy of this TimeSeries where positions are centered around given values
    def get_centered(self, centers):
        if self.type != 'positions':
            raise TypeError("Centering a TimeSeries object only makes sense for positions")
        if len(centers) != len(self.values[0]):
            raise ValueError("The number of columns in the centers must be the same as in the TimeSeries")
        my_centers = self.get_centers()
        deltas = [tuple(my_centers[j][v]-centers[j][v] for v in range(3)) for j in range(len(centers))] # subtract by these
        return TimeSeries([[tuple(row[j][v]-deltas[j][v] for v in range(3)) for j in range(len(centers))] for row in self.values], self.type)

    # global alignment of two TimeSeries objects (needs affine gap penalties)
    @staticmethod
    def align(x, y, distance_metric='euclidean'):
        raise RuntimeError("NOT IMPLEMENTED! Need to find a way to compute gap costs")
        if not isinstance(x, TimeSeries) or not isinstance(y, TimeSeries):
            raise TypeError("Both objects must be TimeSeries objects")
        if x.type != y.type:
            raise TypeError("Both objects must be same type of TimeSeries")
        distance_metric = distance_metric.lower()
        if distance_metric not in DISTANCE[x.type]:
            raise ValueError("Invalid distance metric: %s (options: %s)" % (distance_metric, ','.join(sorted(DISTANCE[x.type].keys()))))
        dist_func = DISTANCE[x.type][distance_metric]
        y = y.get_centered(x.get_centers())
        M = [[(inf,None) for j in range(len(y)+1)] for i in range(len(x)+1)] # (score,backtrack) tuples, where 0 = left, 1 = up, 2 = diag
        X = [[(inf,None) for j in range(len(y)+1)] for i in range(len(x)+1)]
        Y = [[(inf,None) for j in range(len(y)+1)] for i in range(len(x)+1)]

if __name__ == "__main__":
    pass
