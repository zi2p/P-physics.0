import math
import matplotlib.pyplot as plt

f = open('input.txt')

L = float(f.readline()) * 1000
orbit = float(f.readline()) * 1000
f.close()

radius = 1738000
g = 1.62
G = 6.6743015151515152 * (10 ** -11)
M_moon = (g * 1738000 ** 2) / G

array_x = []
array_y = []
array_t = []
array_h = []
array_a = []

def height(x, y):
    return math.sqrt(x ** 2 + y ** 2) - radius

def g(h):
    return G * M_moon / (radius + h) ** 2

def moon():
    dx = 10
    x = -radius
    y = 0
    array_x.append(x)
    array_y.append(y)

    while x < radius:
        x += dx
        y = math.sqrt(radius ** 2 - x ** 2)
        array_x.append(x)
        array_y.append(y)

    while x > -radius:
        x -= dx
        y = -math.sqrt(radius ** 2 - x ** 2)
        array_x.append(x)
        array_y.append(y)

    while x < 0:
        x += dx
        y = math.sqrt(radius ** 2 - x ** 2)
        array_x.append(x)
        array_y.append(y)

def clear():
    array_x.clear()
    array_y.clear()
    array_t.clear()
    array_h.clear()
    array_a.clear()

def append(x, y, t, h, ay):
    array_x.append(x)
    array_y.append(y)
    array_t.append(t)
    array_h.append(h)
    array_a.append(ay)

def print_graphics(stage):
    plt.title("Траектория полета")
    plt.plot(array_x, array_y, color="blue")
    plt.savefig("#" + str(stage) + " Орбита.png")
    plt.show()

def full_a(ax, ay):
    return math.sqrt(ax ** 2 + ay ** 2)

def full_V(Vx, Vy):
    return math.sqrt(Vx ** 2 + Vy ** 2)

def tangent_angle_cos(x, y):
    return x / math.sqrt(x ** 2 + y ** 2)

def perpendicular_angle_cos(x, y):
    return y / math.sqrt(x ** 2 + y ** 2)

def run(dt):
    clear()
    m_spacecraft = 2200
    t = 0
    Vx = 0
    Vy = 0
    ax = 0
    ay = 0
    m_fuel = 4000
    y = radius
    x = 0

    while height(x, y) < orbit / 3 and ay < 29.43 and m_fuel > 0 and (height(x, y) > 0 or t == 0):
        t += dt
        F = 3660 * 8 * dt
        F_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt
        ay = F - F_mg
        Vy += ay * dt
        y += Vy * dt
        m_fuel -= 8 * dt
        append(x, y, t, height(x, y), ay)
    print_graphics(1)

    while height(x, y) < orbit and full_a(ax, ay) < 29.43 and m_fuel > 0 and height(x, y) > 0:
        t += dt
        Fx = 3660 * 6 * dt * math.cos(math.radians(30))
        Fy = 3660 * 6 * dt * math.sin(math.radians(30))
        Fx_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * (tangent_angle_cos(x, y))
        Fy_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * math.sqrt(1 - tangent_angle_cos(x, y) ** 2)
        ax = Fx - Fx_mg
        ay = Fy - Fy_mg
        Vx += ax * dt
        Vy += ay * dt
        x += Vx * dt
        y += Vy * dt
        m_fuel -= 6 * dt
        append(x, y, t, height(x, y), full_a(ax, ay))
    print_graphics(2)

    # торможение
    while Vy > 5.5 and full_a(ax, ay) < 29.43 and m_fuel > 0 and height(x, y) > 100:
        t += dt
        Fx = 3660 * 4 * dt * math.cos(math.radians(math.degrees(math.acos(perpendicular_angle_cos(x, y)))) - 90)
        Fy = 3660 * 4 * dt * math.sin(math.radians(math.degrees(math.acos(perpendicular_angle_cos(x, y)))) - 90)
        Fx_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * tangent_angle_cos(x, y)
        Fy_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * math.sqrt(1 - tangent_angle_cos(x, y) ** 2)
        ax = Fx - Fx_mg
        ay = Fy - Fy_mg
        Vx += ax * dt
        Vy += ay * dt
        x += Vx * dt
        y += Vy * dt
        m_fuel -= 4 * dt
        append(x, y, t, height(x, y), full_a(ax, ay))
    print_graphics(3)

    # разгон до первой космической по орбите
    while full_V(Vx, Vy) < math.sqrt(G * M_moon / (radius + orbit)) and full_a(ax, ay) < 29.43:
        t += dt
        Fx = 3660 * 5.5 * dt * math.cos(math.radians(math.degrees(math.acos(tangent_angle_cos(x, y)))) - 90)
        Fy = 3660 * 5.5 * dt * math.sin(math.radians(math.degrees(math.acos(tangent_angle_cos(x, y)))) - 90)
        Fx_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * tangent_angle_cos(x, y)
        Fy_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * math.sqrt(1 - tangent_angle_cos(x, y) ** 2)
        ax = Fx - Fx_mg
        ay = Fy - Fy_mg
        Vx += ax * dt
        Vy += ay * dt
        x += Vx * dt
        y += Vy * dt
        m_fuel -= 5.5 * dt
        append(x, y, t, height(x, y), full_a(ax, ay))
    print_graphics(4)

    # расчёт времени ождания
    t_os = 0
    x_os = 0
    y_os = radius + orbit
    while math.fabs(75945 - x_os) > 100 or math.fabs(1793390 - y_os) > 100 and t_os < 1000:
        t_os += dt
        x_os += full_V(Vx, Vy) * dt * math.cos(full_V(Vx, Vy) * t_os / (radius + orbit))
        y_os -= full_V(Vx, Vy) * dt * math.sin(full_V(Vx, Vy) * t_os / (radius + orbit))
    t_wait = t_os
    t_os = 0
    x_os = 0
    y_os = radius + orbit
    while (math.fabs(x - x_os) > 100 or math.fabs(y - y_os) > 100) and t_os < 1000:
        t_os += dt
        x_os += full_V(Vx, Vy) * dt * math.cos(full_V(Vx, Vy) * t_os / (radius + orbit))
        y_os -= full_V(Vx, Vy) * dt * math.sin(full_V(Vx, Vy) * t_os / (radius + orbit))
    m_spacecraft -= 200
    print("До взлёта нужно подождать " + str(round(math.fabs(t - t_os + t_wait), 0)) + " секунд.")
    clear()
    moon()

    # позиция торможения
    while math.fabs(x + 144027) > 100 or t < 4000:
        t += dt
        t_os += dt
        x += full_V(Vx, Vy) * dt * math.cos(full_V(Vx, Vy) * t_os / (radius + orbit))
        y += full_V(Vx, Vy) * dt * -math.sin(full_V(Vx, Vy) * t_os / (radius + orbit))
        append(x, y, t, orbit, g)
    print_graphics(5)
    Vx = math.sqrt(G * M_moon / (radius + orbit)) * math.cos(
        math.radians(math.degrees(math.acos(tangent_angle_cos(x, y)))) - 90)
    Vy = math.sqrt(G * M_moon / (radius + orbit)) * math.sin(
        math.radians(math.degrees(math.acos(tangent_angle_cos(x, y)))) - 90)
    clear()

    # торможение по горизорнтали
    while height(x, y) > 0 and full_a(ax, ay) < 29.43 and Vx > 100 and Vy < 10:
        t += dt
        Fx = -3660 * 5 * dt * math.cos(math.radians(0))
        Fy = 3660 * 5 * dt * math.sin(math.radians(0))
        Fx_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * tangent_angle_cos(x, y)
        Fy_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * math.sqrt(1 - tangent_angle_cos(x, y) ** 2)
        ax = Fx - Fx_mg
        ay = Fy - Fy_mg
        Vx += ax * dt
        Vy += ay * dt
        x += Vx * dt
        y += Vy * dt
        m_fuel -= 5 * dt
        append(x, y, t, height(x, y), full_a(ax, ay))
    print_graphics(6)

    # торможение по вертикали
    while height(x, y) > 0 and full_a(ax, ay) < 29.43 and Vy < 0:
        t += dt
        Fx = -3660 * 7 * dt * math.cos(math.radians(90))
        Fy = 3660 * 7 * dt * math.sin(math.radians(90))
        Fx_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * tangent_angle_cos(x, y)
        Fy_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * math.sqrt(1 - tangent_angle_cos(x, y) ** 2)
        ax = Fx - Fx_mg
        ay = Fy - Fy_mg
        Vx += ax * dt
        Vy += ay * dt
        x += Vx * dt
        y += Vy * dt
        m_fuel -= 7 * dt
        append(x, y, t, height(x, y), full_a(ax, ay))
    print_graphics(7)

    # торможение по горизорнтали
    while height(x, y) > 0 and full_a(ax, ay) < 29.43 and Vx > 0.6:
        t += dt
        Fx = -3660 * 7 * dt * math.cos(math.radians(0))
        Fy = 3660 * 7 * dt * math.sin(math.radians(0))
        Fx_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * tangent_angle_cos(x, y)
        Fy_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * math.sqrt(1 - tangent_angle_cos(x, y) ** 2)
        ax = Fx - Fx_mg
        ay = Fy - Fy_mg
        Vx += ax * dt
        Vy += ay * dt
        x += Vx * dt
        y += Vy * dt
        m_fuel -= 7 * dt
        append(x, y, t, height(x, y), full_a(ax, ay))
    print_graphics(8)

    while height(x, y) > 4460 and full_a(ax, ay) < 29.43:
        t += dt
        Fx = 0
        Fy = 0
        Fx_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * tangent_angle_cos(x, y)
        Fy_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * math.sqrt(1 - tangent_angle_cos(x, y) ** 2)
        ax = Fx - Fx_mg
        ay = Fy - Fy_mg
        Vx += ax * dt
        Vy += ay * dt
        x += Vx * dt
        y += Vy * dt
        append(x, y, t, height(x, y), full_a(ax, ay))
    print_graphics(9)

    # торможение по вертикали
    while height(x, y) > 0 and full_a(ax, ay) < 29.43 and Vy < 0:
        t += dt
        Fx = -3660 * 2.7 * dt * math.cos(math.radians(90))
        Fy = 3660 * 2.7 * dt * math.sin(math.radians(90))
        Fx_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * tangent_angle_cos(x, y)
        Fy_mg = g(height(x, y)) * (m_fuel + m_spacecraft) * dt * math.sqrt(1 - tangent_angle_cos(x, y) ** 2)
        ax = Fx - Fx_mg
        ay = Fy - Fy_mg
        Vx += ax * dt
        Vy += ay * dt
        x += Vx * dt
        y += Vy * dt
        m_fuel -= 2.7 * dt
        append(x, y, t, height(x, y), full_a(ax, ay))
    print_graphics(10)

run(0.001)

s=2
for i in range(2,27):
    s *= 28-i
print(s)
