
import math
import matplotlib.pyplot as plt

start_y = 5.308

start_x = 0
step_x = 0.05
end_x = 300

x_range = [start_x]

while x_range[-1] < end_x:
    next_x = x_range[-1] + step_x
    x_range.append(next_x)

dy = lambda x, y: 0.000167716*y*(188.121 - y)
y = lambda x: (188.121 * 5.308) / ((188.121 - 5.308) * math.exp(-188.121 * 0.000167716 * x) + 5.308)

def true_function(function, x_range):
    results = []
    for x in x_range:
        y = function(x)
        results.append(y)
    return results

def euler(function, y_start, x_range):
    results = []
    y = y_start;
    for x in x_range:
        y_step = function(x, y)

        y += y_step * step_x

        results.append(y)

    return results

def improved_euler(function, y_start, x_range):
    results = []
    y = y_start;
    for x in x_range:
        k1 = function(x, y)
        u = y + step_x * k1
        k2 = function(x + step_x, u)

        y += step_x * 0.5 * (k1 + k2)

        results.append(y)

    return results

def roupe_kutta(function, y_start, x_range):
    results = []
    y = y_start;
    for x in x_range:
        k1 = function(x, y)
        k2 = function(x + step_x/2, y + k1*step_x/2)
        k3 = function(x + step_x/2, y + k2*step_x/2)
        k4 = function(x, y + k3 * step_x)

        y += (step_x/6)*(k1 + 2*k2 + 2*k3 + k4)

        results.append(y)

    return results

true_data = true_function(y, x_range)
euler_data = euler(dy, start_y, x_range)
improved_euler_data = improved_euler(dy, start_y, x_range)
roupe_kutta_data = roupe_kutta(dy, start_y, x_range)

print("Final true: " + str(true_data[-1]))
print("Final euler: " + str(euler_data[-1]))
print("Final improved euler: " + str(improved_euler_data[-1]))
print("Final roupe kutta: " + str(roupe_kutta_data[-1]))

fig, ax = plt.subplots()
ax.plot(x_range, true_data, 'b,-')
ax.plot(x_range, euler_data, 'r,-')
ax.plot(x_range, improved_euler_data, 'g,-')
ax.plot(x_range, roupe_kutta_data, 'y,-')
ax.legend(["True function", "Euler approximation", "Improved Euler approximation", "Roupe Kutta approximation"])
fig.suptitle("Approximated Population vs Time")
plt.xlabel("Time")
plt.ylabel("Population (millions)")
plt.show()

