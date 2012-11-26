# -*- coding: utf-8 *-*
import numpy as np
import matplotlib.pyplot as plt

multiple = (20, .1)  # example B, b, i.e., B, beta for multiple ss

resolution = 1000
Daf = lambda x, B, b: x / (1 - x) * np.exp(-B * x / (1 + b))
xmin = 0.01
xmax = .98
dx = (xmax - xmin) / resolution
xl = list(np.arange(xmin, xmax, dx))
#extendage = list(np.arange(xmin / 100., xmax / 100., dx / 100.))
#xl.extend(extendage)
B = multiple[0]
b = multiple[1]

Dal = [Daf(x, B, b) for x in xl]
yf = lambda x, Da: np.log(x / Da / (1 - x))
yl = [yf(x, Da) for (x, Da) in zip(xl, Dal)]

fig1 = plt.figure(1, figsize=(11, 8.5))
ax11 = fig1.add_subplot(1, 1, 1)
ax11.plot(Dal, xl, label=r'$x_s$')
ax11.plot(Dal, yl, label=r'$y_s$')
ax11.set_xlabel(r'$Da$')
ax11.set_ylabel(r'$x_s$ or $y_s$')
ax11.legend()
#ax11.set_yscale('log')

minB = 0.01
maxB = 30
minb = 0.01
maxb = 10
dB = (maxB - minB) / resolution
db = (maxb - minb) / resolution
Bl = np.arange(minB, maxB, dB)
bl = np.arange(minb, maxb, db)
(B, b) = np.meshgrid(Bl, bl)

bifurcations = B > 3 + b + 2 * np.lib.scimath.sqrt(2 + b)
unique = B < 4*(1 + b)

uniqueUnstable = unique & bifurcations
multipleStable = np.logical_not(unique) & np.logical_not(bifurcations)
fig2 = plt.figure(2)
ax21 = fig2.add_subplot(1, 1, 1)
ax21.imshow(uniqueUnstable, origin='lower', cmap='Greys')
B = multiple[0]
b = multiple[1]
fig2.suptitle('Unique steady states that can demonstrate oscilatory instability \n'\
 + r'(e.g., $B=%.1f$, $\beta=%.1f$ is a combination that shows multiple, stable SSs.)' % (B, b))
numticks = 5

actualtofigure = lambda actual, maxval: float(actual) * resolution / maxval
figuretoactual = lambda figure, maxval: float(figure) / resolution * maxval
a2f = actualtofigure
f2a = figuretoactual
locations = [int(x) for x in np.arange(0, resolution, resolution / numticks)]
Blabels = [f2a(x, maxB) for x in locations]
blabels = [f2a(x, maxb) for x in locations]
ax21.set_xticks(locations)
ax21.set_xticklabels(Blabels)
ax21.set_yticks(locations)
ax21.set_yticklabels(blabels)
ax21.set_xlabel(r'$B$')
ax21.set_ylabel(r'$\beta$')

ax21.scatter([a2f(multiple[0], maxB)], [a2f(multiple[1], maxb)], color='w')

print max(bl)
ax21.set_xlim((0, a2f(max(Bl), maxB)))
ax21.set_ylim((0, a2f(max(bl), maxb)))

fig3 = plt.figure(3, figsize=(11, 8.5))
ax31 = fig3.add_subplot(1, 1, 1)
tmin = 0.01
tmax = 44
timeresolution = 10000
dt = (tmax - tmin) / timeresolution
tl = list(np.arange(tmin, tmax, dt))
Da = 4
B = multiple[0]
b = multiple[1]
xp = lambda x, y: -x + Da * (1 - x) * np.exp(y)
yp = lambda x, y: -(1 + b) * y + B * Da * (1 - x) * np.exp(y)

initials = [(0.01, 0.01)]
initials = [(4, 4)]
for initial in initials:
    xold = initial[0]
    yold = initial[1]
    xl = []
    yl = []
    for t in tl:
        x = xold + dt * xp(xold, yold)
        y = yold + dt * yp(xold, yold)
        xl.append(x)
        yl.append(y)
        xold = x
        yold = y
    ax31.plot(tl, xl)
    ax31.plot(tl, yl)
    #ax31.plot(xl, yl)

plt.show()
