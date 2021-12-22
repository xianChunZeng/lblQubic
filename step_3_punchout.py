from matplotlib import pyplot
import numpy
import sys
sys.path.append('submodules/qubic/src')
from qubic.qcvv.punchout import c_punchout


punchout=c_punchout(qubitid='vna',calirepo='submodules/qchip')
fcenter=punchout.opts['qchip'].getfreq(sys.argv[1]+'.readfreq')
fbw=10e6
fstart=fcenter-fbw/2
fstop=fcenter+fbw/2
fx=numpy.linspace(fstart,fstop,200)
punchout.run(103,fx=fx,attens=numpy.arange(-35,0.2,3.0),maxvatatten=0)

fig1=pyplot.figure(figsize=(10,10))
sub=fig1.subplots(2,2)
fig1.suptitle(sys.argv[1])
punchout.plotamp(sub[0,0])
punchout.plotang(sub[0,1])
punchout.plotampdiff(sub[1,0])
punchout.plotangdiff(sub[1,1])
fig2=pyplot.figure()
sub2=fig2.subplots()
punchout.plotmaxangdiff(sub2)
fig2.suptitle(sys.argv[1])
fig3=pyplot.figure()
ax3=fig3.subplots(2,2)
ax3[0,0].plot(punchout.s11.real.T,punchout.s11.imag.T)
ax3[1,0].plot(abs(punchout.s11.T))
ax3[1,1].plot(numpy.angle(punchout.s11.T))
pyplot.show()

