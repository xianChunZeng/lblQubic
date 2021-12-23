from matplotlib import pyplot
import numpy
from qubic.qcvv.vna import c_vna

vna=c_vna(qubitid='vna',calirepo='submodules/qchip')
#fx,fxstep=numpy.linspace(6.53e9,6.85e9,1000,retstep=True)
#	fx=numpy.linspace(6.69e9,6.72e9,1000)
lor=vna.opts['wiremap'].lor
bw=vna.opts['chassis'].fsample
#fx=numpy.linspace(6.2e9,6.7e9,2000)
fx=numpy.linspace(lor-bw/2,lor+bw/2,2000)
vna.seqs(fx,t0=692e-9,amp=0.5)
vna.run(100)

#vna.peaks(0.2)
	
fig1=pyplot.figure(figsize=(15,8))
sub1=fig1.subplots(3,2).reshape(6)
vna.iqapplot(sub1)
#fig2=pyplot.figure(figsize=(15,8))
#pyplot.plot(vna.fx[0:len(vna.p)],vna.p)
#pyplot.plot(vna.fx[vna.p2],vna.p[vna.p2],'x')
#	fig3=pyplot.figure(figsize=(8,8))
#	ax3=fig3.subplots()
#	vna.iqplot(ax3,hexbin=False,marker='.')
#pyplot.savefig('fig1.pdf')
pyplot.show()

