# -*- coding: utf-8 *-*
import numpy as np
import matplotlib.pyplot as plt
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
A = np.ones((resolution, resolution))
B = (mu ** 2 + a ** 3 - a ** 2) / a ** 2
C = mu ** 2 / a
eig1 = (-B + np.lib.scimath.sqrt(B ** 2 - 4 * A * C)) / 2 / a
eig2 = (-B - np.lib.scimath.sqrt(B ** 2 - 4 * A * C)) / 2 / a

real = (eig1.imag == 0) | (eig2.imag == 0)
real = np.logical_not(real)  # clearly I have made some mistakes
imaginary = np.logical_not(real)
neg = (eig1.real <= 0) | (eig2.real <= 0)
negreal = neg | real
negimag = neg | imaginary

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
ax1.set_title('oscilatory \n (Have imaginary eigenvalues.)')
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
ax2.set_title('stable \n (Eigenvalues have negative real part.)')
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
ax3.imshow(negimag, cmap='Greys', origin='lower')
ax3.set_title('stable oscilatory \n (Eigenvalues have negative real part and are imaginary.')
ax3.scatter([0.2 * resolution / maxmu], [0.4 * resolution / maxa], color='k')
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
ax4.imshow(negreal, cmap='Greys', origin='lower')
ax4.set_title('stable non-oscilatory \n (Eigenvalues are negative pure-real.)')
ax4.scatter([a2f(0.1, maxmu)], [a2f(0.4, maxa)], color='k')
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

plt.tight_layout()
#fig1.savefig('hw6_1_f1.pdf')
plt.show()

fig2 = plt.figure(2, figsize=(11, 8.5))

ax21 = fig2.add_subplot(2, 1, 1)
tmin = 0.01
tmax = 120
timeresolution = 10000
dt = (tmax - tmin) / timeresolution
tl = list(np.arange(tmin, tmax, dt))
xp = lambda x, y, mu, a: mu - (1 + a) * x + x ** 2 * y
yp = lambda x, y: x - x ** 2 * y
mu = 0.2
a  = 0.4
initials = [(0, 0), (0, 1), (2, 2.5)]
patterns = ['k--', 'k:', 'k-']
for (initial, pattern)in zip(initials, patterns):
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
    ax21.plot(xl, yl, pattern)
ax21.legend([r'$x_0=%.1f$, $y_0=%.1f$' % initial for initial in initials])
ax21.set_title(r'phase-plane trajectories for $\tau_{max}=%.0f$, with the parameters $\mu=%.2f$, $\alpha=%.2f$' % (tmax, mu, a))
ax21.set_xlabel(r'$x$')
ax21.set_ylabel(r'$y$')

ax22 = fig2.add_subplot(2, 1, 2)
ax22.plot(tl, xl, 'k-')
ax22.plot(tl, yl, 'k--')
ax22.legend([r'$x$', r'$y$'])
ax22.set_title(r'$x_0=%.1f$, $y_0=%.1f$' % initial)
ax22.set_xlabel(r'$\tau$')
ax22.set_ylabel(r'$x$ or $y$')
plt.tight_layout()
#fig2.savefig('hw6_1_f2.pdf')
#plt.show()
