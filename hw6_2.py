# -*- coding: utf-8 *-*
import numpy as np
import matplotlib.pyplot as plt
import sys
write = sys.stdout.write
flush = sys.stdout.flush

multiple = (8, .4)  # example (B, b) for multiple SS's
single = (18, 4)  # example (B, b) for single, Hopf-able SS's

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

#find bifurcation points
xs1 = .5 + np.lib.scimath.sqrt(B ** 2 - 4 * B * (1 + b)) / 2 / B
xs2 = .5 - np.lib.scimath.sqrt(B ** 2 - 4 * B * (1 + b)) / 2 / B
Das1 = Daf(xs1, B, b)
Das2 = Daf(xs2, B, b)
ax11 = fig1.add_subplot(1, 1, 1)
ax11.plot(Dal, xl, color='k', label=r'$x_s$')
ax11.scatter([Das1, Das2], [xs1, xs2], color='k')
#label them
annotation1 = r'$Da=%.3f$, $x_s=%.3f$' % (Das1, xs1)
ax11.annotate(annotation1, xy=(Das1, xs1), xytext=(Das1+.01, xs1))
annotation2 = r'$Da=%.3f$, $x_s=%.3f$' % (Das2, xs2)
ax11.annotate(annotation2, xy=(Das2, xs2), xytext=(Das2+.01, xs2))

ax11.set_xlabel(r'$Da$')
ax11.set_ylabel(r'$x_s$')
#ax11.legend()
fig1.suptitle(r'Problem 2: Multiple steady states are possible for $B=%.1f$, $\beta=%.1f$' % (B, b))
#ax11.set_yscale('log')
#plt.show()

minB = 0.01
maxB = 30
minb = 0.01
maxb = 10
dB = (maxB - minB) / resolution
db = (maxb - minb) / resolution
Bl = np.arange(minB, maxB, dB)
bl = np.arange(minb, maxb, db)
(B, b) = np.meshgrid(Bl, bl)

bifurcations = B < 3 + b + 2 * np.lib.scimath.sqrt(2 + b)
unique = B > 4*(1 + b)

uniqueUnstable = unique | bifurcations
multipleStable = np.logical_not(unique) & np.logical_not(bifurcations)
#plt.show()
filename = 'hw6_2_f1.pdf'
print 'saving', filename
fig1.savefig(filename)

fig2 = plt.figure(2, figsize=(11, 8.5))
ax21 = fig2.add_subplot(1, 1, 1)
ax21.imshow(uniqueUnstable, origin='lower', cmap='Greys')

B = multiple[0]
b = multiple[1]
fig2.suptitle('Problem 2: Unique steady states that can demonstrate oscilatory instability \n'\
 + r'(e.g., $B=%.1f$, $\beta=%.1f$ satisfies this, but ' % (single[0], single[1]) \
 + r"$B=%.1f$, $\beta=%.1f$ shows multiple SS's.)" % (multiple[0], multiple[1]))

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
ax21.scatter([a2f(single[0], maxB)], [a2f(single[1], maxb)], color='k')
ax21.set_xlim((0, a2f(max(Bl), maxB)))
ax21.set_ylim((0, a2f(max(bl), maxb)))
#plt.show()
filename = 'hw6_2_f2.pdf'
print 'saving', filename
fig2.savefig(filename)

fig3 = plt.figure(3, figsize=(11, 8.5))
ax31 = fig3.add_subplot(1, 1, 1)
tmin = 0.01
tmax = 44
timeresolution = 10000
dt = (tmax - tmin) / timeresolution
tl = list(np.arange(tmin, tmax, dt))
Da = (0.080 + 0.041) / 2
B = single[0]
b = single[1]
xpf = lambda x, y: -x + Da * (1 - x) * np.exp(y)
ypf = lambda x, y: -(1 + b) * y + B * Da * (1 - x) * np.exp(y)

xmax = 1
ymax = 2
xs = np.arange(.01, xmax, 0.098)
ys = np.arange(.01, ymax, 0.097)
initials = []
for x0 in xs:
    initials.append((x0, 0.01))
    initials.append((x0, 1.9))
for y0 in ys:
    initials.append((0.01, y0))
    initials.append((1, y0))

#initials = [(.01, .1), (.01, .5), (.5, 2.8), (.94, 1)]
#patterns = ['k-', 'k--', 'k:', 'k.']
#for (initial, pattern) in zip(initials, patterns):
for (initial, count) in zip(initials, range(len(initials))):
    pattern = 'k-'
    xold = x0 = initial[0]
    yold = y0 = initial[1]
    xl = []
    yl = []
    maxslope = 0
    for t in tl:
        xp = xpf(xold, yold)
        yp = ypf(xold, yold)
        if abs(xp) > maxslope:
            maxslope = abs(xp)
        if abs(yp) > maxslope:
            maxslope = abs(yp)
        x = xold + dt * xp
        y = yold + dt * yp
        y = max([0, y])
        x = max([0, x])
        xl.append(x)
        yl.append(y)
        xold = x
        yold = y
    write('\rIC '+ str(count+1) + ' of ' + str(len(initials)) + ': x=%.3f, y=%.3f' % initial)
    flush()
    if maxslope > 30:
        print 'bigslope'
        pattern = 'r-'
    #ax31.plot(tl, xl, pattern, markersize=3, label=r'$B=%.2f$, $\beta=%.2f$' % (B, b))
    #ax31.plot(tl, yl, pattern, label=r'$B=%.2f$, $\beta=%.2f$' % (B, b))
    ax31.plot(xl, yl, pattern)
print ''
#ax31.set_ylim((0, 2))  # set if doing x vs t  or  y vst
ax31.set_xlim((0, xmax))  # set if doing phase portrait
ax31.set_ylim((0, ymax))  # set if doing phase portrait
ax31.set_xlabel(r'dimensionless conversion, $x$')
ax31.set_ylabel(r'dimensionless temperature, $y$')
fig3.suptitle('Problem 2: Phase portrait\n' +
r'$Da=%.4f$, $B=%.1f$, $\beta=%.1f$' % (Da, B, b))
#ax31.legend()
filename = 'hw6_2_f3.pdf'
print 'saving', filename
fig3.savefig(filename)
#plt.show()
