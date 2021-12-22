import sys
sys.path.append('../..')
from qubic.qcvv.rabioptimize import c_rabioptimize
from qubic.qcvv.ramsey import c_ramsey
from qubic.qcvv.repeatgate import c_repeatgate
from qubic.qcvv.dragalpha import c_dragalpha
from qubic.qubic import envset,qubicrun,heralding
from qubic.qubic.squbic import c_qchip
import json
import datetime
if __name__=="__main__":
    import datetime
    t0=datetime.datetime.now()
    qubitids=['Q%d'%i for i in range(2,3)]
    qubitids=['Q%d'%i for i in range(5,8)]
    qubitids=['Q%d'%i for i in range(0,-1,-1)]
    qubitids=['Q4','Q5']
    qubitids=['Q6','Q7','Q0']
    qubitids=['Q1']
    qubitids=['Q4']
    qubitids=['Q%d'%i for i in (1,2,3,4,5,7)]
    qubitids=['Q5']
    qubitids=['Q0']
    qubitids=['Q6']
    qubitids=['Q%d'%i for i in (1,3)]
    calirepo='../../../../qchip'
    nsample=20
#    with open('t2.json') as jfile:
#        qubitcfg=json.load(jfile)
#    qchip=c_qchip(qubitcfg)
#    stepclass=[c_repeatgate]
    t0=datetime.datetime.now()
    qchip=None
    plot=False
    stepclass=[
            #c_rabioptimize:dict(gmixs=None,nsample=nsample,plot=True,freadxtol=200e3,frabixtol=200e3,arabixtol=0.05,areadxtol=0.05),
            #c_rabioptimize:dict(gmixs=None,nsample=nsample,plot=True),
            (c_ramsey,dict(gmixs='.',nsample=nsample,plot=plot,elementstep=4e-9)),
            (c_ramsey,dict(gmixs='.',nsample=nsample,plot=plot,elementstep=20e-9)),
            (c_ramsey,dict(gmixs='.',nsample=nsample,plot=plot,elementstep=100e-9)),
            (c_repeatgate,dict(gmixs='.',nsample=nsample,plot=plot,nsteps=5)),
            (c_dragalpha,dict(gmixs='.',nsample=nsample*10,plot=plot)),
            ]
    for qubitid in qubitids:
        index=0
        for c_step,steppara in stepclass:
            obj=c_step(qubitid=qubitid,calirepo=calirepo,qchip=qchip,**steppara)
            updatedict=obj.optimize(**steppara)
            qchip=c_qchip(obj.opts['qchip'].updatecfg(updatedict,'%s_t%d.json'%(qubitid,index)))
            t1=datetime.datetime.now()
            print('optimstep ',obj.__class__,'qubitid',qubitid,'time',t1-t0,'step',index)
            index+=1



#    
#    rabioptimize=c_rabioptimize(qubitid=qubitid,calirepo=calirepo)
##    rabioptidict=rabioptimize.optimize(nsample=50,disp=3)
#    rabioptidict={}
#    qchip=c_qchip(rabioptimize.opts['qchip'].updatecfg(rabioptidict,'t1.json'))
#    ramsey=c_ramsey(qubitid=qubitid,calirepo=calirepo,qchip=qchip)
#    ramseyoptidict=ramsey.optimize(nsample=50)
#    print('ramseyoptidict',ramseyoptidict)
#    qchip=c_qchip(ramsey.opts['qchip'].updatecfg(ramseyoptidict,'t2.json'))
#    dragalpha=c_dragalpha(qubitid=qubitid,calirepo=calirepo,qchip=qchip)
#    dragoptidict=dragalpha.optimize(nsample=50)
#    print('dragoptidict',dragoptidict)
#    qchip=c_qchip(dragalpha.opts['qchip'].updatecfg(dragoptidict,'t3.json'))
#    repeatgate=c_repeatgate(qubitid=qubitid,calirepo=calirepo,qchip=qchip)
#    rgatedict=repeatgate.optimize(nsample=50)
#    print('roptidict',rgatedict)
#    repeatgate.opts['qchip'].updatecfg(rgatedict,wfilename='test.json')
#
