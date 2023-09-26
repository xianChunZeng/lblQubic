import sys
from matplotlib import pyplot
import numpy
from chipcalibration.gatesim import *
from scipy.optimize import curve_fit
import numpy.linalg


#  cnot =   xx(0,pxtgt)  zz(pzctl,0)  zz(0,-pztgt) cr zz(0,pztgt) @s
#  cr   =   zz(0,pztgt) zz(-pzctl,0) xx(0,-pxtgt) cnot zz(0,-pztgt)
def cr(pzctl,pztgt,pxtgt):
    return zz(0,pztgt)@xx(0,pxtgt)@zz(-pzctl,0)@cnot()@zz(0,-pztgt)

def fun2(xyscan,pzctl,pztgt,pxtgt,zzdeglist):#=[[246,147],[314,234],[336,37]]):
    k,s=sini()
    for iz,(zzdeg1,zzdeg2) in enumerate(zzdeglist):
        if iz!=0:
            s=xx()@s
        s=zzdeg(zzdeg1,zzdeg2)@s
    state=cr(pztgt=pztgt,pzctl=pzctl,pxtgt=pxtgt)@s

    xyval=numpy.zeros((len(xyscan),len(state)))
    for i,p in enumerate(xyscan):
        xyval[i,:]=tomoxy(state,p)[:,0]
    return xyval

def fitfun(xyscan,pzctl,pztgt,pxtgt,zzdeglist=[[246,147],[314,234],[336,37]]):
    y=fun2(xyscan[0:int(len(xyscan)/4)],pzctl=pzctl,pztgt=pztgt,pxtgt=pxtgt,zzdeglist=zzdeglist)
    return y.flatten('F')

if __name__=="__main__":
    meas=numpy.load('crfullxydata.npz',allow_pickle=True)
    xymeas=meas['v'].tolist()
    xyscan=meas['xyrot']#*180/numpy.pi
    xymeas4=numpy.zeros((len(xymeas),len(xyscan)))
    xymeas4[0,:]=xymeas[('0','0')]#[:,2]
    xymeas4[1,:]=xymeas[('0','1')]#[:,2]
    xymeas4[2,:]=xymeas[('1','0')]#[:,2]
    xymeas4[3,:]=xymeas[('1','1')]#[:,2]
    xdata=numpy.tile(xyscan.T,4)
    ydata=xymeas4.T.flatten('F')
    #popt,pcov=curve_fit(fitfun,xdata=xdata,ydata=ydata,p0=(0,0.1,0.2),bounds=(([-numpy.pi,-numpy.pi,-numpy.pi],[numpy.pi,numpy.pi,numpy.pi])))#([-360,-360,-360],[360,360,360])))
    popt,pcov=curve_fit(fitfun,xdata=xdata,ydata=ydata,p0=(0,0,0),bounds=(([-360,-360,-360],[360,360,360])))
            #(([-numpy.pi,-numpy.pi,-numpy.pi],[numpy.pi,numpy.pi,numpy.pi]))
           # )
    xyfit=fitfun(xdata,*popt)
    pyplot.figure('ydata')
    pyplot.plot(ydata,'*')
    pyplot.plot(xyfit)
    print('popt deg',popt*180/numpy.pi,'perr',numpy.sqrt(numpy.diag(pcov))*180/numpy.pi)
    print('popt rad',popt)
    lx=len(xdata)//4
    pyplot.figure('xy')
    color=['r','g','b','k']
    for i in range(4):
        yi=ydata.reshape((4,-1))[i,:]
        pyplot.plot(xdata[0:lx],yi,'*',color=color[i])
        fi=xyfit.reshape((4,-1))[i,:]
        pyplot.plot(xdata[0:lx],fi,color=color[i])

    zero=fitfun(xdata,*(00,00,-200))#(120.46,46.68,-134))
    pyplot.show()

