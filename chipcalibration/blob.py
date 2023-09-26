import itertools
import scipy.optimize
import numpy
import matplotlib
from matplotlib import pyplot
def scalemax(ax):    
    ax.autoscale_view()
    xlim=ax.get_xlim()
    ylim=ax.get_ylim()
    xylim=[min(min(xlim),min(ylim)),max(max(xlim),max(ylim))]
    ax.set_aspect('equal')
    ax.set_xlim(xylim)
    ax.set_ylim(xylim)
class c_blob:
    def __init__(self,mean,covar):
        self.mean=numpy.mat(mean)
        self.covar=numpy.mat(covar)
        self.calc()
    def calc(self):
        try:
            c=numpy.linalg.cholesky(self.covar)
        except:
            print('covar',self.covar)
            print('eigval',numpy.linalg.eigvals(self.covar))
            c=numpy.zeros(self.covar.shape)
            
        u,s,v=numpy.linalg.svd(c)
        self.smn=numpy.mat(numpy.diag(1/s))
        self.angle = numpy.pi+numpy.arctan2(u[0,1],u[0,0])#u[1], u[0])
        self.angledeg = self.angle * 180/numpy.pi
        self.x0=self.mean[0,0]
        self.y0=self.mean[0,1]
        self.a,self.b=s
        self.u=u*numpy.mat(numpy.array([(1,0),(0,numpy.linalg.det(u))]))#*v.T
        self.t=self.u*self.smn
    def dist(self,blob,plot=False):
        bi=blob.tonormalizedof(self)
#        bs=self.tonormalizedof(self)
#        mean=bs.unittomeanrotof(bi)
        mean=bi.unittomeanrotof()
        x0=mean[0,0]
        y0=mean[0,1]
        a=bi.a
        b=bi.b
        def func3(xx):
            x=xx[0]
            y=xx[1]
            k=xx[2]
            return [(x/(k*a))**2+(y/(k*b))**2-1
            ,(x-x0)**2+(y-y0)**2-k**2
            ,(y0-y)-a**2/b**2*y/x*(x0-x)
            ]
        x,y,k= scipy.optimize.fsolve(func3, [x0/2, y0/2,1])
        if plot:
            t=numpy.arange(0,2*numpy.pi,0.001)
            xcir=x0+k*numpy.cos(t)
            ycir=y0+k*numpy.sin(t)
            xe=k*a*numpy.cos(t)
            ye=k*b*numpy.sin(t)
            fig=pyplot.figure(figsize=(10,10))
            ax=fig.subplots()
            ax.plot(xcir,ycir)
            ax.plot(xe,ye)
            ax.set_aspect('equal')
        return k

    def ellipse(self,k=1,**kwargs):
        opts=dict(edgecolor='b',facecolor='none',linewidth=2,alpha=0.5)
        opts.update({kk:v for kk,v in kwargs.items() if kk in opts})
        ell=matplotlib.patches.Ellipse((self.x0,self.y0),2*self.a*k,2*self.b*k,self.angledeg)
        ell.set_edgecolor(opts['edgecolor'])
        ell.set_facecolor(opts['facecolor'])
        ell.set_linewidth(opts['linewidth'])
        ell.set_alpha(opts['alpha'])
        return ell
    def ellipsexy(self,k=1):
        t=numpy.arange(0,numpy.pi*2,0.001)
        x=self.x0+self.a*k*numpy.cos(t)*numpy.cos(self.angle)-self.b*k*numpy.sin(t)*numpy.sin(self.angle)
        y=self.y0+self.b*k*numpy.sin(t)*numpy.cos(self.angle)+self.a*k*numpy.cos(t)*numpy.sin(self.angle)
        return (x,y)
    def tonormalizedof(self,blob):
        mean=(blob.t.T*(self.mean-blob.mean).T).T
        covar=blob.t.T*self.covar*blob.t
        return c_blob(mean,covar)
    def tomeanof(self,blob):
        mean=self.mean-blob.mean#(blob.u*(self.mean-blob.mean).T).T
        covar=self.covar
        return c_blob(mean,covar)
    def tomeanrotof(self,blob):
        mean=(blob.u.T*(self.mean-blob.mean).T).T
        covar=blob.u.T*self.covar*blob.u
        return c_blob(mean,covar)
    def unittomeanrotof(self):
        mean=(self.u.T*(-self.mean).T).T
        return mean
    def __str__(self):
        para=['mean: %s'%str(self.mean),'covar: %s'%str(self.covar),'a:%s'%str(self.a),'b:%s'%str(self.b),'angledeg:%s'%str(self.angledeg%(360))]
        return '\n'.join(para)

def ell(a,b,theta):
    s=numpy.mat(numpy.array(((a,0),(0,b))))
    r=numpy.mat(numpy.array(((numpy.cos(theta),numpy.sin(theta)),(-numpy.sin(theta),numpy.cos(theta)))))
    covar=r*s*s.T*r.T
    return covar
def gmixdist(gmix):
    #    labels={int(k):v for k,v in gmix.labels.items()}
    #labels={k:v for k,v in labels.items() if k in [0,1,2]}
    blobs={}
    for ilabel,label in gmix.labels.items():
        blob=c_blob(mean=gmix.means_[ilabel],covar=gmix.covariances_[ilabel])
        blobs[label]=blob
    distij={}
    for i,j in itertools.combinations(blobs,2):
        bi=blobs[i] # 1
        bj=blobs[j] # 0
        distij[tuple(sorted([i,j]))]=bi.dist(bj)
    return distij
if __name__=="__main__":
    if 0:
        gmix=numpy.load("Q6_ef_gmix.npz")
        distij=gmixdist(gmix)
        #print(distij)
#        pyplot.show()
    if 1:
        blobs={}
        k=1
        covar1=ell(a=6*k,b=1*k,theta=numpy.pi/6)
        covar2=ell(a=14*k,b=9*k,theta=-numpy.pi*1/4)
        mean0=(1,6)
        mean1=(-3,-14)
        blobs['1']=c_blob(mean=mean0,covar=covar1)#numpy.mat(numpy.array([[4,2],[2,7]])))
        blobs['0']=c_blob(mean=mean1,covar=covar2)#numpy.mat(numpy.array([[14,-2],[-2,3.9]])))
    if 0:
        fig=pyplot.figure(figsize=(10,10))
        ax=fig.subplots(2,2)
        bnorm3={}
        for k,blob in blobs.items():
            ell0=blob.ellipse()
            #b1=blob.tomeanof(blobs['1'])
            #ell1=b1.ellipse()
            #b2=blob.tomeanrotof(blobs['1'])
            #ell2=b2.ellipse()
            #b3=blob.tonormalizedof(blobs['1'])
            #print('b3',k,b3)
            #ell3=b3.ellipse()
            #bnorm3[k]=b3
            ax[0,0].plot(blob.mean[0,0],blob.mean[0,1],color='b',marker='$%s$'%k,markersize=20)
            ax[0,0].add_patch(ell0)
        for i,j in distij:
            if '0' in (i,j) and '1' in (i,j):
                axi=ax[0,1]
            elif '2' in (i,j) and '1' in (i,j):
                axi=ax[1,0]
            elif '2' in (i,j) and '0' in (i,j):
                axi=ax[1,1]
            
            elli0=blobs[i].ellipse()
            ellj0=blobs[j].ellipse()
            elli=blobs[i].ellipse(distij[(i,j)])
            ellj=blobs[j].ellipse(distij[(i,j)])
            axi.plot(blobs[i].mean[0,0],blobs[i].mean[0,1],color='b',marker='$%s$'%i,markersize=20)
            axi.plot(blobs[j].mean[0,0],blobs[j].mean[0,1],color='b',marker='$%s$'%j,markersize=20)
            axi.add_patch(elli0)
            axi.add_patch(ellj0)
            axi.add_patch(elli)
            axi.add_patch(ellj)

        #    ax[1,0].add_patch(ell2)
        #    ax[1,1].add_patch(ell3)
        scalemax(ax[0,0])
        scalemax(ax[0,1])
        scalemax(ax[1,0])
        scalemax(ax[1,1])
        pyplot.show()
#    for k,b1 in bnorm3.items():
#        print('dist b1',bnorm3['1'].dist(b1))

    
    if 1:
        k=0.37
        k=0.1221726*2#*numpy.sqrt(2)
        k=0.58#2.205
        k=0.9217
        k=1.3788601648310348

        covar11=ell(a=6*k,b=1*k,theta=numpy.pi/6)
        covar22=ell(a=14*k,b=9*k,theta=-numpy.pi*1/4)
        b11=c_blob(mean=mean0,covar=covar11)#numpy.mat(numpy.array([[4,2],[2,7]])))
        b22=c_blob(mean=mean1,covar=covar22)#numpy.mat(numpy.array([[4,2],[2,7]])))
        ell11=b11.ellipse()
        ell22=b22.ellipse()
        fig=pyplot.figure(figsize=(10,10))
        ax=fig.subplots()
        ax.add_patch(ell11)
        ax.add_patch(ell22)
        x,y=b11.ellipsexy()
        ax.plot(x,y,'r')
        x,y=b22.ellipsexy()
        ax.plot(x,y,'k')
        scalemax(ax)
        pyplot.show()
