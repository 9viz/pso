import numpy as np

data = np.loadtxt("OUT.methane")
for i in range(len(data)):
    data[i][0::2] *= 2 * np.pi / 360 # DEG to RAD.

def rot(axis, ang, point):
    """AXIS may be 'x', 'y', 'z'."""
    sin = np.sin(ang)
    cos = np.cos(ang)
    rot = []
    if axis == "x":
        rot = [
            [1, 0, 0],
            [0,  cos, sin],
            [0, -sin, cos]
        ]
    elif axis == "y":
        rot = [
            [cos, 0, sin],
            [0, 1, 0],
            [-sin, 0, cos]
        ]
    elif axis == "z":
        rot = [
            [cos, sin, 0],
            [-sin, cos, 0],
            [0, 0, 1]
        ]
    rot = np.asarray(rot, dtype=np.float64)
    return np.matmul(rot, point)

# TODO: This conversion is not perfect, not all bond angles are equal?
# Other way:
#     print("C", "0", "0", "0", sep="\t")
#     h = [0,0,dist]
#     print("H", *h, sep="\t")
#     h = rot("y", ang, h)
#     print("H", *h, sep="\t")
#     h = rot("z", ang, h)
#     print("H", *h, sep="\t")
#     h = rot("z", ang, h)
#     print("H", *h, sep="\t")
def frame(f, ang, dist):
    print("C", "0", "0", "0", sep="\t", file=f)
    h = [0,0,dist]
    print("H", *rot("y", ang/2, h), sep="\t", file=f)
    print("H", *rot("y", -ang/2, h), sep="\t", file=f)
    h = [0,0,-dist]
    print("H", *rot("x", ang/2, h), sep="\t", file=f)
    print("H", *rot("x", -ang/2, h), sep="\t", file=f)

fs = [None] * 5
for i in range(5):
    # Global best.
    if i == 4:
        fs[i] = open("meth-global.xyz", "w")
    else:
        fs[i] = open(f"meth-ind-{i+1}.xyz", "w")

for n,i in enumerate(data):
    i = np.split(i,len(i)/2)
    for j in range(len(i)):
        f = fs[j]
        print("5", file=f)
        print(f"Generation {n}", file=f)
        frame(f, i[j][0], i[j][1])

for i in fs: i.close()
