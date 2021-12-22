import sys
sys.path.append('submodules/qubic/src')
from matplotlib import pyplot
from qubic.qcvv.alignment import c_alignment


alignment=c_alignment(qubitid='alignment',calirepo='submodules/qchip',debug=False,sim=True)
tlo=540e-9
alignment.seqs(tlo=tlo,mon_sel0=0,mon_sel1=1)
fig1=pyplot.figure('alignment',figsize=(15,8))
live=1
if live:
    while live:
        alignment.run(100)
        print('adcminmax',alignment.opts['chassis'].adcminmax())
        pyplot.clf()
        alignment.plotmonwithfft(fig1)
        pyplot.pause(0.1)
        pyplot.show(block=False)
else:
    alignment.run(100)
    print('adcminmax',alignment.opts['chassis'].adcminmax())
    #alignment.plot(fig1)
    alignment.plotmonwithfft(fig1)
    pyplot.show()
