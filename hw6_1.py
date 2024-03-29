# https://github.com/tsbertalan/504hw6
# -*- coding: utf-8 *-*
import numpy as np
import matplotlib.pyplot as plt
from sys import exit
lc = (.2, .5)  # example limit cycle parameter combinaiton
sf = (.4, .8)  # example stable focus parameter combinaiton
#so = (.4, .8)
#so = (.8, .8)
minmu = .01
maxmu = 1
mina = .01
maxa = 2
resolution = 1000
dmu = (maxmu - minmu) / resolution
da = (maxa - mina) / resolution
mul = np.arange(minmu, maxmu, dmu)
al = np.arange(mina, maxa, da)
(mu, a) = np.meshgrid(mul, al)

def eigs(mu, a):
    '''Calculate two eigenvalues'''
    A = 1
    B = (mu ** 2 + a ** 3 - a ** 2) / a ** 2
    C = mu ** 2 / a
    eig1 = (-B + np.lib.scimath.sqrt(B ** 2 - 4 * A * C)) / 2 / A
    eig2 = (-B - np.lib.scimath.sqrt(B ** 2 - 4 * A * C)) / 2 / A
    return (eig1, eig2)

(eig1, eig2) = eigs(mu, a)

real = (eig1.imag == 0) & (eig2.imag == 0)
real = np.logical_not(real)  # clearly I have made some mistakes
imaginary = np.logical_not(real)
neg = (eig1.real <= 0) & (eig2.real <= 0)
stable = np.logical_not(neg)  # clearly I have made some mistakes
unstable = np.logical_not(stable)
stablefocus = stable | imaginary
limitcycle = unstable | imaginary

#negimag = neg | imaginary  # I don't understand why | is required here, not &.
                            # the operators seem to be switched (!)
#negreal = neg | real

# convert between figure coordinates and actual mu or a values
# assumes that the minumum value for either parameter is about 0.
actualtofigure = lambda actual, maxval: float(actual) * resolution / maxval
figuretoactual = lambda figure, maxval: float(figure) / resolution * maxval
a2f = actualtofigure
f2a = figuretoactual

fig1 = plt.figure(1, figsize=(11, 8.5))
numaxes = 3
numticks = 5
ax1 = fig1.add_subplot(2, 2, 1)
ax1.imshow(imaginary, cmap='Greys', origin='lower')
ax1.set_title('Oscilatory')
locations = [int(x) for x in np.arange(0, resolution, resolution / numticks)]
mulabels = [f2a(x, maxmu) for x in locations]
alabels  = [f2a(x, maxa)  for x in locations]
ax1.set_xticks(locations)
ax1.set_xticklabels(mulabels)
ax1.set_yticks(locations)
ax1.set_yticklabels(alabels)
ax1.set_xlabel(r'$\mu$')
ax1.set_ylabel(r'$\alpha$')

ax2 = fig1.add_subplot(2, 2, 2)
ax2.imshow(neg, cmap='Greys', origin='lower')
ax2.set_title('Stable')
locations = [int(x) for x in np.arange(0, resolution, resolution / numticks)]
alabels = [float(x) / resolution * maxa for x in locations]
mulabels = [float(x) / resolution * maxmu for x in locations]
ax2.set_xticks(locations)
ax2.set_xticklabels(mulabels)
ax2.set_yticks(locations)
ax2.set_yticklabels(alabels)
ax2.set_xlabel(r'$\mu$')
ax2.set_ylabel(r'$\alpha$')

ax3 = fig1.add_subplot(2, 2, 3)
ax3.imshow(limitcycle, cmap='Greys', origin='lower')
ax3.set_title('Limit Cycle')
ax3.scatter([a2f(lc[0], maxmu)], [a2f(lc[1], maxa)], color='k')
ax3.set_xlim([0, resolution])
ax3.set_ylim([0, resolution])
locations = [int(x) for x in np.arange(0, resolution, resolution / numticks)]
alabels = [float(x) / resolution * maxa for x in locations]
mulabels = [float(x) / resolution * maxmu for x in locations]
ax3.set_xticks(locations)
ax3.set_xticklabels(mulabels)
ax3.set_yticks(locations)
ax3.set_yticklabels(alabels)
ax3.set_xlabel(r'$\mu$')
ax3.set_ylabel(r'$\alpha$')

ax4 = fig1.add_subplot(2, 2, 4)
ax4.imshow(stablefocus, cmap='Greys', origin='lower')
ax4.set_title('Stable Focus')
ax4.scatter([a2f(sf[0], maxmu)], [a2f(sf[1], maxa)], color='k')
ax4.set_xlim([0, resolution])
ax4.set_ylim([0, resolution])
locations = [int(x) for x in np.arange(0, resolution, resolution / numticks)]
alabels = [float(x) / resolution * maxa for x in locations]
mulabels = [float(x) / resolution * maxmu for x in locations]
ax4.set_xticks(locations)
ax4.set_xticklabels(mulabels)
ax4.set_yticks(locations)
ax4.set_yticklabels(alabels)
ax4.set_xlabel(r'$\mu$')
ax4.set_ylabel(r'$\alpha$')

fig1.suptitle('Problem 1: Regimes')
plt.tight_layout()
filename = 'hw6_1_f1.pdf'
print 'saving', filename
fig1.savefig(filename)
#plt.show()
#exit()
ax1es = range(4)
ax2es = range(4)
figs = range(4)

#For plotting nullclines:
#from xdot:
y1 = lambda x, mu, a: ((1 + a)*x - mu) / x ** 2
#from ydot:
y2 = lambda x, mu, a: 1 / x
xmin = 0.01
xmax = 5
dx = (xmax - xmin) / resolution
xl_nullclines = list(np.arange(xmin, xmax, dx))

plt.tight_layout()

for (fignum, parameters, regime) in zip([2, 3], [lc, sf], ['Stable Limit Cycle', 'Stable Focus']):
    mu, a = parameters
    figs[fignum] = plt.figure(fignum, figsize=(11, 8.5))

    ax1es[fignum] = figs[fignum].add_subplot(2, 1, 1)
    tmin = 0.01
    tmax = 200
    timeresolution = 8000
    dt = (tmax - tmin) / timeresolution
    tl = list(np.arange(tmin, tmax, dt))
    xp = lambda x, y, mu, a: mu - (1 + a) * x + x ** 2 * y
    yp = lambda x, y: x - x ** 2 * y
    initials = [(2, 2.5), (0, 0)]
    patterns = ['k+', 'k.']
    for (initial, pattern) in zip(initials, patterns):
        xold = initial[0]
        yold = initial[1]
        xl = []
        yl = []
        for t in tl:
            x = xold + xp(xold, yold, mu, a) * dt
            y = yold + yp(xold, yold) * dt
            xl.append(x)
            yl.append(y)
            xold = x
            yold = y
        ax1es[fignum].plot(xl, yl, pattern, label=r'$x_0=%.2f$, $y_0=%.2f$' % initial, markersize=3)
#    nullclines:
    y1l = [y1(x, mu, a) for x in xl_nullclines]
    y2l = [y2(x, mu, a) for x in xl_nullclines]
    ax1es[fignum].plot(xl_nullclines, y1l, 'b-', label=r"nullcline from $x'$")
    ax1es[fignum].plot(xl_nullclines, y2l, 'r--', label=r"nullcline from $y'$")
    ax1es[fignum].set_ylim((0, max(yl) + 1))
    ax1es[fignum].set_xlim((0, xmax))
    ax1es[fignum].legend()
    ax1es[fignum].set_title(r'phase-plane trajectories for $\tau_{max}=%.0f$' % (tmax,))
    ax1es[fignum].set_xlabel(r'$x$')
    ax1es[fignum].set_ylabel(r'$y$')

    #predict period from eigenvalues of linearization
    (eig1, eig2) = eigs(mu, a)
    (eig1, eig2) = eigs(mu, a)
    eigstr1 = '%.3f +%.3f j' % (eig1.real, eig1.imag)
    eigstr2 = '%.3f %.3f j' % (eig2.real, eig2.imag)
    period = 2 * np.pi / eig1.imag

    ax2es[fignum] = figs[fignum].add_subplot(2, 1, 2)
    ax2es[fignum].plot(tl, xl, 'k-')
    ax2es[fignum].plot(tl, yl, 'k--')
    ax2es[fignum].legend([r'$x$', r'$y$'])
    ax2es[fignum].set_title(r'Conversion over time starting from $x_0=%.1f$, $y_0=%.1f$' % initial)
    ax2es[fignum].set_xlabel(r'$\tau$ (tick marks are separated by periods $T=2\pi/\omega=%.1f$, as predicted by the linearization)' % period)
    ax2es[fignum].set_ylabel(r'$x$ or $y$')

    title = 'Problem 1: ' + regime + \
    r', with the parameters $\mu=%.2f$, $\alpha=%.2f$' % (mu, a) + ', \n' +\
    r'and the eigenvalues ' + eigstr1 + ' ; ' + eigstr2
    figs[fignum].suptitle(title)
    filename = 'hw6_1_f%i.pdf' % fignum
    print 'saving', filename
    tickspots = list(np.arange(0, tmax, period))
    plt.xticks(tickspots, ['%.2f' % x for x in tickspots], rotation=17)
    figs[fignum].savefig(filename)
#    plt.tight_layout()

#plt.show()
