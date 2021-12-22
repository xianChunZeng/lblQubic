import time
import sys
sys.path.append('submodules/qubic/src')
from matplotlib import pyplot
import numpy
from qubic.qcvv.vna import c_vna

vna=c_vna(qubitid='vna',calirepo='submodules/qchip')
#fx,fxstep=numpy.linspace(6.53e9,6.85e9,1000,retstep=True)
#	fx=numpy.linspace(6.69e9,6.72e9,1000)
lor=vna.opts['wiremap'].lor
bw=vna.opts['chassis'].fsample
#fx=numpy.linspace(6.2e9,6.7e9,2000)
    #fx=numpy.linspace(lor,lor+bw/2,1000)
    #fx=numpy.linspace(6.3e9,6.8e9,1000)
    #fx=numpy.linspace(lor-bw/2,lor+bw/2,2000)
if len(sys.argv)==1:
    fcenter=6.53e9
    fbw=0.5e9
    fcenter=lor
    fbw=bw
else:
    fbw=30e6
    fcenter=vna.opts['qchip'].getfreq(sys.argv[1]+'.readfreq')
fx=numpy.linspace(fcenter-fbw/2,fcenter+fbw/2,1000)
vna.seqs(fx,amp=0.1)
#vna.run(1)

#vna.peaks(0.2)
	
#fig1=pyplot.figure(figsize=(15,8))
#vna.iqapplot(fig1)
#fig2=pyplot.figure(figsize=(15,8))
#pyplot.plot(vna.fx[0:len(vna.p)],vna.p)
#pyplot.plot(vna.fx[vna.p2],vna.p[vna.p2],'x')
#	fig3=pyplot.figure(figsize=(8,8))
#	ax3=fig3.subplots()
#	vna.iqplot(ax3,hexbin=False,marker='.')
#pyplot.savefig('fig1.pdf')
#vna.iqapplot(ax1)
fig1=pyplot.figure('vna',figsize=(15,8))
live=1
if live:
    pyplot.ion()
    while 1:
        vna.run(100,includegmm=False)
        print('adcminmax',vna.opts['chassis'].adcminmax())
        dt=vna.fit()
        print('dt estimate: delay=current','+' if dt<0 else '-', abs(dt))
        pyplot.clf()
        ax1=fig1.subplots(3,2)
#        ax1[0,1].set_ylim([0.0,0.01])#[-10,1e6])
        vna.iqapplot(ax1)
        time.sleep(0.1)
        pyplot.pause(0.1)
        pyplot.show(block=False)
                #pyplot.ion()
else:
    vna.run(1000)
    print('adcminmax',vna.opts['chassis'].adcminmax())
    ax1=fig1.subplots(3,2)
    vna.iqapplot(ax1)
    pyplot.show()
