import json
class chipcfg():
    def __init__(self):
        self.readaligngate=[]
        pass
    def newchip(self):
        self.chipdict=dict(Qubits={},Gates={})
        self.addqubit(dict(M0={}))
        self.addqubit(dict(vna=self.freqs(readfreq=0)))
        self.addqubit(dict(alignment=self.freqs(readfreq=0)))
        self.addgate(self.markergate('M0'))
        self.addgate(self.rdrvread('vna',twidth=4e-6))
        self.addgate(self.rdrvread('alignment',twidth=1e-6))

    def singlequbitgates(self,qubits,gates=['X90','X90_ef','rabi','rabi_ef','Y_90']):
        for q in qubits:
            self.addqubit({q:self.freqs(readfreq=0,freq=0,freq_ef=0)})
            for gate in gates:
                self.addgate(getattr(self,gate)(q))
    def twoqubitsgates(self,qubitspairs,gates=['CR','cnot']):
        for ctl,tgt in qubitspairs:
            for gate in gates:
                self.addgate(getattr(self,gate)(ctl=ctl,tgt=tgt))
    def addqubit(self,qdict):
        for q,d in qdict.items():
            if q not in self.chipdict['Qubits']:
                self.chipdict['Qubits'][q]={}
            for f,v in d.items():
                self.chipdict['Qubits'][q][f]=v
    def addgate(self,gdict):
        for g,l in gdict.items():
            if g not in self.chipdict['Gates']:
                self.chipdict['Gates'][g]=[]
            if not (isinstance(l,list) or isinstance(l,tuple)):
                l=[l]
            for p in l:
                self.chipdict['Gates'][g].append(p)
    def markergate(self,mname,twidth=1.0e-6):
        gatename=mname+'mark'
        destname=mname+'.mark'
        Gates={gatename:dict(dest=destname,pcarrier=0.0,t0=0.0,fcarrier=0.0,amp=1.0,twidth=twidth,env=self.markenv())}
        return Gates
    def markenv(self):
        return [dict(envfunc="mark",paradict={})]
    def readfreq(self,name,readfreq=0):
        freq=name+'.readfreq'
        Qubits={name:dict(readfreq=readfreq)}
        return Qubits
    def rdrvread(self,name,t0=0,amp=0.5,twidth=1e-6,rdrvenv=None,readenv=None):
        gatename=name+'read'
        rdrvenv=self.envsel(rdrvenv) #self.square() if rdrvenv is None else rdrvenv
        readenv=self.envsel(readenv) #self.square() if readenv is None else readenv
        freq=name+'.readfreq'
        rdrvdict=dict(dest=name+'.rdrv',pcarrier=0,fcarrier=freq,t0=t0,amp=amp,twidth=twidth,env=[rdrvenv])
        readdict=dict(dest=name+'.read',pcarrier=0,fcarrier=freq,t0=t0,amp=1.0,twidth=twidth,env=[readenv])
        Gates={gatename:[rdrvdict,readdict]}
        self.readaligngate.append(gatename)
        return Gates
    def env_square(self,phase=0,amplitude=1,**kwargs):
        pdict=dict(env_func='square',paradict=dict(phase=0.0,amplitude=1.0))
        pdict.update(kwargs)
        return pdict
    def env_DRAG(self,alpha=0,sigmas=3,delta=-270e6,**kwargs):
        pdict=dict(env_func='DRAG',paradict=dict(alpha=alpha,sigmas=sigmas,delta=delta))
        pdict.update(kwargs)
        return pdict
    def env_cos_edge_square(self,ramp_fraction=0.25,**kwargs):
        pdict=dict(env_func='cos_edge_square',paradict=dict(ramp_fraction=ramp_fraction))
        pdict.update(kwargs)
        return pdict
    def envsel(self,env,**kwargs):
        if env is None or env=='square':
            pdict=self.env_square(**kwargs)
        elif env=='DRAG':
            pdict=self.env_DRAG(**kwargs)
        elif env=='cos_edge_square':
            pdict=self.env_cos_edge_square(**kwargs)
        else:
            psict={}
        return pdict
#    def qubit(self,name,**kwargs):
#        Qbuits=self.freqs(**kwargs)
#        pass
    def freqs(self,readfreq,freq=None,freq_predicted=None,freq_ef=None,**kwargs):
        return {k:v for k,v  in dict(freq=freq,readfreq=readfreq,freq_predicted=freq_predicted,freq_ef=freq_ef).items() if v is not None}
    def rabi_ef(self,name,t0=0,amp=0.5,twidth=32e-9,env='DRAG'):
        gatename=name+'rabi_ef'
        freq=name+'.freq_ef'
        env=self.envsel(env) #self.square() if rdrvenv is None else rdrvenv
        gatedict=dict(dest=name+'qdrv',pcarrier=0,fcarrier=freq,t0=t0,amp=amp,twidth=twidth,env=[env])
        Gates={gatename:[gatedict]}
        return Gates
    def rabi(self,name,t0=0,amp=0.5,twidth=32e-9,env='DRAG'):
        gatename=name+'rabi'
        freq=name+'.freq'
        env=self.envsel(env) #self.square() if rdrvenv is None else rdrvenv
        gatedict=dict(dest=name+'qdrv',pcarrier=0,fcarrier=freq,t0=t0,amp=amp,twidth=twidth,env=[env])
        Gates={gatename:[gatedict]}
        return Gates
    def X90(self,name,t0=0,amp=0.5,twidth=32e-9,env='DRAG'):
        gatename=name+'X90'
        freq=name+'.freq'
        env=self.envsel(env) #self.square() if rdrvenv is None else rdrvenv
        gatedict=dict(dest=name+'qdrv',pcarrier=0,fcarrier=freq,t0=t0,amp=amp,twidth=twidth,env=[env])
        Gates={gatename:[gatedict]}
        return Gates
    def X90_ef(self,name,t0=0,amp=0.5,twidth=32e-9,env='DRAG'):
        gatename=name+'X90_ef'
        freq=name+'.freq_ef'
        env=self.envsel(env) #self.square() if rdrvenv is None else rdrvenv
        gatedict=dict(dest=name+'qdrv',pcarrier=0,fcarrier=freq,t0=t0,amp=amp,twidth=twidth,env=[env])
        Gates={gatename:[gatedict]}
        return Gates
    def Z_90pulse(self,name,t0=0):
        Gates=dict(pcarrier="-numpy.pi/2.0",fcarrier=name+".freq",t0=t0)
        return Gates
    def Z90pulse(self,name,t0=0):
        Gates=dict(pcarrier="numpy.pi/2.0",fcarrier=name+".freq",t0=t0)
        return Gates
    def Y_90pulse(self,name,tX90):
        pulse=[self.Z_90pulse(name),name+"X90",self.Z90pulse(name,tX90)]
        return pulse
    def Y_90(self,name,tX90=32e-9):
        Gates={name+"Y-90":self.Y_90pulse(name,tX90)}
        return Gates
    def CR(self,ctl,tgt,twidth=0,amp=0.5,env=None):
        gatename='%s%sCR'%(ctl,tgt)
        dest='%s.qdrv'%(ctl)
        freq='%s.freq'%(tgt)
        env= self.envsel('cos_edge_square',ramp_length=32e-9) if env is None else env
        gdict={gatename:dict(dest=dest,pcarrier=0,t0=0,twidth=twidth,fcarrier=freq,amp=amp,env=[env])}
        return gdict
    def cnot(self,ctl,tgt,tcr=256e-9,pztgt=0,pzctl=0,axtgt=0.5,txtgt=32e-9,xenv='DRAG'):
        gatename='%s%scnot'%(ctl,tgt)
        CRgatename='%s%sCR'%(ctl,tgt)
        tgtfreq='%s.freq'%(tgt)
        ctlfreq='%s.freq'%(ctl)
        tgtdest='%s.qdrv'%(tgt)
        ctldest='%s.qdrv'%(ctl)
        xenv= self.envsel(xenv)
        pulses=[dict(pcarrier=pztgt,t0=0,fcarrier=tgtfreq)
                ,CRgatename
                ,dict(pcarrier=-pztgt,t0=tcr,fcarrier=tgtfreq)
                ,dict(pcarrier=-pzctl,t0=tcr,fcarrier=ctlfreq)
                ,dict(dest=tgtdest,pcarrier=0,fcarrier=tgtfreq,t0=tcr,amp=axtgt,twidth=txtgt,env=xenv)
                ]
        gdict={gatename:pulses}
        return gdict

if __name__=="__main__":
    a=chipcfg()
    a.newchip()
    qubits=['Q%i'%i for i in range(8)]
    qubitspairs=[('Q2','Q1'),('Q3','Q2')]
    a.singlequbitgates(qubits)
    a.twoqubitsgates(qubitspairs)
    fwrite=open('test.json','w')
    json.dump(a.chipdict,fwrite,indent=4)
    fwrite.close()
    json.dumps(a.chipdict,indent=4)
