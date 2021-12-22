import sys
sys.path.append('submodules/qubic/src')
from matplotlib import pyplot
import numpy
from qubic.qcvv.vnaqubit import c_vnaqubit

qubitid=sys.argv[1]
vnaqubit=c_vnaqubit(qubitid=qubitid,calirepo='submodules/qchip')
loq=vnaqubit.opts['wiremap'].loq
bw=vnaqubit.opts['chassis'].fsample

fx=numpy.linspace(loq-bw/2,loq+bw/2,1000)
fig1=pyplot.figure(figsize=(15,8))
sub1=fig1.subplots(3,2)

vnaqubit.seqs(fx,ef=False)
vnaqubit.run(200)

vnaqubit.iqapplot(qubitid,sub1)
fig2=pyplot.figure(figsize=(15,8))
sub2=fig2.subplots()
vnaqubit.iqplot(qubitid,sub2)
#	pyplot.savefig('fig1.pdf')
pyplot.show()

