from matplotlib import pyplot
from scipy.optimize import curve_fit
import numpy

def _linepara(t,a,b):
    return a*t+b
def _quadratic(t,a,b,c):
    return a*(t-b)**2+c
def _sinpara(t,a,f,p,o):
    return a*numpy.sin(2*numpy.pi*f*t+p)+o
# rabi fit
def _sinsquare(t, a, f, p, c = 0): # set default constant offset to zero
    return a * numpy.sin(2*numpy.pi*f * t + p)**2 + c
# t1 and spin-echo fit
def _exppara(t, a, tau, c = 0): # set default constant offset to zero
    return a * numpy.exp(-t/tau) + c
# ramsey fit
def _expsin(t, a, f, p, c, tau): # input, amplitude, decay parameter, angular freq, offset phase, constant offset
    return a * numpy.exp(- t/tau) * numpy.sin(2*numpy.pi*f * t + p) + c # ideally c = 0.5
def _rb_func(x,p,a): # fit this to a power law
    return a*p**x+0.5

def _Zcavity(w,R,Q,w0):
    Z=1.0*R/(1+1j*Q*(w/w0-w0/w))
    return Z
def _cavs11(f,A,R,Q,f0,t0,a,b):
    w0=2*numpy.pi*f0
    w=2*numpy.pi*f
    z=_Zcavity(w,R,Q,w0)
    s11=A*numpy.exp(-2*1j*w*t0)*(z-1)/(z+1)+a+b*1j
    sr=numpy.hstack([s11.real,s11.imag])
    return sr
def fitcavs11(f,s11,**kwargs):
    try:
        popt,pcov = curve_fit(_cavs11, f,numpy.hstack([s11.real,s11.imag]),**kwargs)
        yfit2=_cavs11(f,*popt)
        yfit=yfit2[0:int(len(yfit2)/2)]+1j*yfit2[int(len(yfit2)/2):]
        rsqr=calcrsqr(s11,yfit)
    except:
        popt,pcov,yfit,rsqr=errorret(f,7)
    return dict(popt=popt,pcov=pcov,yfit=yfit,rsqr=rsqr)
def s11est(f,s11,f0est):
    f0est=f0est
    t0est=0e-9
    Qest=60000
    Rest=1
    Aest=max(abs(s11))-min(abs(s11))
    aest=s11.real.mean()
    best=s11.imag.mean()
    p0=[Aest,Rest,Qest,f0est,t0est,aest,best]
    return p0
def s11fit5(f,s11):
    fits=[]
    l=len(s11)
#    print(s11.shape)
#    print(s11est(f,t,f[numpy.argmax(abs(s11))]))
    fits.append(fitcavs11(f,s11,p0=s11est(f,s11,f0est=f[numpy.argmax(abs(s11))])))
    fits.append(fitcavs11(f,s11,p0=s11est(f,s11,f0est=f[int(l/4)])))
    fits.append(fitcavs11(f,s11,p0=s11est(f,s11,f0est=f[int(l/2)])))
    fits.append(fitcavs11(f,s11,p0=s11est(f,s11,f0est=f[int(3/4)])))
    return fits[numpy.argmax(numpy.array([i['rsqr'] for i in fits]))]

def errorret(x,npara):
    popt=numpy.ones(npara)*numpy.inf
    pcov=numpy.ones((npara,npara))*numpy.inf
    yfit=numpy.ones(len(x))*numpy.inf
    rsqr=1.0
    return popt,pcov,yfit,rsqr

def fitline(x,y,**kwargs):
    try:
        popt, pcov = curve_fit(_linepara, x,y,**kwargs)
        yfit=_linepara(x,*popt)
        rsqr=calcrsqr(y,yfit)
    except:
        popt,pcov,yfit,rsqr=errorret(x,2)
    #return popt,pcov,yfit,rsqr
    return dict(popt=popt,pcov=pcov,yfit=yfit,rsqr=rsqr)
def fitquadratic(x,y,**kwargs):
    #pyplot.plot(x,y)
    #pyplot.show()
    numpy.savez('data.npz',numpy.array((x,y)))
    try:
        popt,pcov = curve_fit(_quadratic,x,y,**kwargs)
        yfit=_quadratic(x,*popt)
        rsqr=calcrsqr(y,yfit)
    except:
        popt,pcov,yfit,rsqr=errorret(x,3)
    #pyplot.plot(x,y)
    #pyplot.plot(x,yfit)
    #pyplot.show()
    return dict(popt=popt,pcov=pcov,yfit=yfit,rsqr=rsqr)
def fitsin(x,y,**kwargs):
    try:
        popt, pcov = curve_fit(_sinpara, x,y,**kwargs)
        if popt[0]<0:
            popt[0]=-popt[0]
            popt[2]=(popt[2]+numpy.pi)%(numpy.pi*2)
        yfit=_sinpara(x,*popt)
        #print(popt,numpy.sqrt(numpy.diag(pcov)))
        #print(kwargs['p0'])
        rsqr=calcrsqr(y,yfit)
    except:
        popt,pcov,yfit,rsqr=errorret(x,4)
    return dict(popt=popt,pcov=pcov,yfit=yfit,rsqr=rsqr)

def fitsinsquare(x,y,**kwargs):
    try:
        popt, pcov = curve_fit(_sinsquare, x,y,**kwargs)
        yfit=_sinsquare(x,*popt)
        rsqr=calcrsqr(y,yfit)
    except:
        popt,pcov,yfit,rsqr=errorret(x,4)
    return dict(popt=popt,pcov=pcov,yfit=yfit,rsqr=rsqr)
def fitexp(x,y,**kwargs):
    try:
        popt, pcov = curve_fit(_exppara, x, y, **kwargs)
        yfit=_exppara(x,*popt)
        rsqr=calcrsqr(y,yfit)
    except:
        popt,pcov,yfit,rsqr=errorret(x,3)
    #pyplot.figure()
    #pyplot.plot(x,y)
    #pyplot.plot(x,yfit)
    #return popt,pcov,yfit,rsqr
    return dict(popt=popt,pcov=pcov,yfit=yfit,rsqr=rsqr)

def fitexpsin(x,y,**kwargs):
    try:
        popt, pcov = curve_fit(_expsin,x,y,**kwargs)
        yfit=_expsin(x,*popt)
        rsqr=calcrsqr(y,yfit)
    except:
        popt,pcov,yfit,rsqr=errorret(x,5)
    return dict(popt=popt,pcov=pcov,yfit=yfit,rsqr=rsqr)
def calcrsqr(y,yfit):
    ymean=y.mean()
    ssr=numpy.sum(abs(y-yfit)**2)
    sst=numpy.sum(abs(y-ymean)**2)
    if sst==0:
        if ssr==0:
            rsqr=0
        else:
            rsqr=1
    else:
        rsqr=ssr/sst
    return rsqr
#    numpy.sum((yfit-y)**2/abs(yfit))
def sinestimate(x,y,**kwargs):
    spec=abs(numpy.fft.rfft(y-y.mean()))
    df=1.0/(x[-1]-x[0])
    indexmax=numpy.argmax(spec)
    f=indexmax*df
    a=(y.max()-y.min())/2.0
    c=y.mean()
    p= kwargs['p'] if 'p' in kwargs else 0
    return a,f,p,c


if __name__=="__main__":
    x=numpy.arange(0,320e-9,4e-9)
    y=_sinpara(x,0.3,1./128e-9,p=1.2,o=0.3)
    est=sinestimate(x,y)
    print('should be',[0.3,1./128e-9,1.2,0.3])
    popt,pcov,yfit=fitsin(x,y,p0=est)
    print('popt',popt)
    print('pcov',pcov)

    from matplotlib import pyplot
    pyplot.plot(x,y,'.r')
    pyplot.plot(x,yfit)
    pyplot.show()
