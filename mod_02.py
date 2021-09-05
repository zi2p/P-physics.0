import math
import matplotlib.pyplot as plt
import numpy as np

q = 1.6 * 10 ** -19
m = 9.31 * 10 ** -31

f = open("input.txt", "r")
r = float(f.readline()) / 100
R = float(f.readline()) / 100
V0 = float(f.readline())
L = float(f.readline()) / 100
f.close()

global x, y, Vy, U, V, t

x = 0
y = (R + r) / 2
V = V0
Vx = V0
Vy = 0
U = 0
t = 0
dt = (L / V0) / 1000


def now(t, x, y, Vy, U):
    while True:
        t += dt
        Vy += a(U, y) * dt

        y += Vy * dt + (a(U, y) * dt ** 2) / 2
        x += Vx * dt
        V = (Vx ** 2 + Vy ** 2) ** (1 / 2)
        if x > L and y > r:
            U += 0.0000001
            t = 0
            x = 0
            y = (R + r) / 2
        elif y <= r:
            f1 = open("output.txt", "w")
            print(V, " ", t, " ", U)
            f1.write(str(V)+'\n')
            f1.write(str(t)+'\n')
            f1.write(str(U)+'\n')
            f1.close()  # дальше не надо продолжать, мы все нашли
            return


def a(U, y):
    return - q * U / (m * y * math.log(R / r))


def grafs(U, tl):
    y = (R + r) / 2
    Vy = 0
    x = 0
    linspace = np.linspace(0, tl, num = 100000)
    dt = linspace[1] - linspace[0]
    xarr, yarr, Varr, aarr = [], [], [], []
    for t in linspace:
        Vy += a(U, y) * dt
        y += Vy * dt + (a(U, y) * dt ** 2) / 2      # зависимость от времени
        yarr.append(y)
        x += Vx * dt           # зависимость от времени
        xarr.append(x)
        # V = (Vx ** 2 + Vy ** 2) ** (1 / 2)  # зависимость от времени
        Varr.append(Vy)
        aarr.append(a(U, y))
    return xarr, yarr, Varr, aarr, linspace

now(t, x, y, Vy, U)

y = (R + r) / 2
t = 0
Vy = 0
x = 0

with open('output.txt') as f:
    _V = float(f.readline().split()[0])
    _t = float(f.readline().split()[0])
    _U = float(f.readline().split()[0])

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 6))
xarr, yarr, Varr, aarr, linspace = grafs(_U, _t)

#ax1.set_title('Зависимость x от времени')
ax1.plot(xarr, yarr)
ax1.set(ylabel='y', xlabel='x')

#ax2.set_title('Зависимость y от времени')
ax2.plot(linspace, yarr)
ax2.set(ylabel='y', xlabel='t')

#ax3.set_title('Зависимость V от времени')
ax3.plot(linspace, Varr)
ax3.set(ylabel='Vу', xlabel='t')

#ax4.set_title('Зависимость a от времени')
ax4.plot(linspace, aarr)
ax4.set(ylabel='a', xlabel='t')
plt.tight_layout()
plt.show()
