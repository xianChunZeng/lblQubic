import time
from matplotlib import pyplot
import sys
sys.path.append('submodules/qubic/src')
from qubic.qcvv.chevron import c_chevron

t0=time.time()
#qubitid=sys.argv[1]
#if 1:
#for qubitid in ['Q5','Q6']:
#for qubitid in ['Q2']:
for qubitid in ['Q0','Q1','Q2','Q3','Q4','Q7']:
    t1=time.time()
    print('qubitid',qubitid,'start',t1)
    chevron=c_chevron(qubitid=qubitid,calirepo='submodules/qchip',gmixs=None)

    freq=chevron.opts['qchip'].getfreq(qubitid+'.freq')
    fbw=100e6
    if 1:
        freq=chevron.opts['wiremap'].loq
        fbw=chevron.opts['chassis'].fsample
    fsteps=100
    chevron.seqs(fstart=freq-fbw/2,fstop=freq+fbw/2,elementlength=41,elementstep=8e-9,fsteps=fsteps,overlapcheck=False)
    chevron.run(100)
    chevron.fit()

    fig1=pyplot.figure('fig1%s'%qubitid,figsize=((30,10)))
    ax1=fig1.subplots(2,1)
    chevron.plot(ax1[:])
    pyplot.savefig('fig1%s.png'%qubitid)

    fig2=pyplot.figure('fig2%s'%qubitid,figsize=(15,8))
    ax2=fig2.subplots()
    chevron.iqplot(ax2)
    chevron.gmmplot(ax2)
    pyplot.savefig('fig2%s.png'%qubitid)
    print('qubitid',qubitid,time.time()-t1)
#fig3=pyplot.figure('ampfit')
#ax3=fig3.subplots()
#chevron.ampfitplot(ax3)
#fig4=pyplot.figure('rsqr')
#ax4=fig4.subplots()
#chevron.rsqrplot(ax4)
pyplot.show()
