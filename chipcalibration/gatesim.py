import numpy
import re

def deg2rad(p):
    return p*numpy.pi/180.0

def i():
    return numpy.array([[1,0]
        ,[0,1]])
def rzdeg(pdeg):
    return rz(deg2rad(pdeg))
def rz(p):
    return numpy.array([[numpy.exp(-1j*p/2),0]
        ,[0,numpy.exp(1j*p/2)]])

def rxdeg(pdeg):
    return rx(deg2rad(pdeg))
def rx(p):
    c=numpy.cos(p/2)
    s=numpy.sin(p/2)
    return numpy.array([[c,-1j*s],[-1j*s,c]])

def rydeg(pdeg):
    return rydeg(deg2rad(pdeg))
def ry(p):
    return rzdeg(90)@rxdeg(p)@rzdeg(-90)

def yydeg(p1deg,p2deg):
    return yy(deg2rad(p1deg),deg2rad(p2deg))
def yy(p1,p2):
    return numpy.kron(ry(p1),ry(p2))

def tomoxx():
    return yydeg(-90,-90)
def tomoyy():
    return xxdeg(90,90)

def zzdeg(p1deg,p2deg):
    return zz(deg2rad(p1deg),deg2rad(p2deg))
def zz(p1,p2):
    return numpy.kron(rz(p1),rz(p2))

def xxdeg(p1deg=90,p2deg=90):
    return xx(deg2rad(p1deg),deg2rad(p2deg))
def xx(p1=numpy.pi/2,p2=numpy.pi/2):
    return numpy.kron(rx(p1),rx(p2))

def cnot():
    return numpy.array(
        [[1,0,0,0]
        ,[0,1,0,0]
        ,[0,0,0,1]
        ,[0,0,1,0]])

def sini(s00=1,s01=0,s10=0,s11=0):
    ret=([('0','0'),('0','1'),('1','0'),('1','1')],numpy.array([[s00,s01,s10,s11]]).T)
    return ret

def prob(s):
    return numpy.square(abs(s))

def tomoxydeg(s,pdeg):
    return tomoxy(s,deg2rad(pdeg))
def tomoxy(s,p):
    return prob(xx()@zz(p,p)@s)

#def op1(z1p1deg,z1p2deg,z2p1deg,z2p2deg,z3p1deg,z3p2deg):
#    return cnot()@z(z3p1deg,z3p2deg)@x()@z(z2p1deg,z2p2deg)@x()@z(z1p1deg,z1p2deg)
#
#def op2(pdeg,p1deg,p2deg):
#    return x()@z(pdeg,pdeg)@z(p1deg,p2deg)
#
#def op3(z1p1deg,z1p2deg,z2p1deg,z2p2deg,z3p1deg,z3p2deg):
#    return z(z3p1deg,z3p2deg)@x()@z(z2p1deg,z2p2deg)@x()@z(z1p1deg,z1p2deg)
#
#def op4(pdeg,p1deg,p2deg,p4deg):
#    return x()@z(pdeg,pdeg)@z(p1deg,p2deg)@x(0,p4deg)@cnot()@z(0,-p2deg)
#
#
#def fullxymeas(pdeg,p1deg,p2deg,z1p1deg,z1p2deg,z2p1deg,z2p2deg,z3p1deg,z3p2deg,s00,s01,s10,s11):
#    '''fullxymeas Qubit Map:(control,target)'''
#    return numpy.array([prob(op2(p,p1deg,p2deg)@op1(z1p1deg,z1p2deg,z2p1deg,z2p2deg,z3p1deg,z3p2deg)@sini(s00,s01,s10,s11)) for p in pdeg])
#
#def fullxymeas_flat(pdeg,p1deg,p2deg,z1p1deg,z1p2deg,z2p1deg,z2p2deg,z3p1deg,z3p2deg,s00,s01,s10,s11):
#    return fullxymeas(pdeg[0:int(len(pdeg)/4)],p1deg,p2deg,z1p1deg,z1p2deg,z2p1deg,z2p2deg,z3p1deg,z3p2deg,s00,s01,s10,s11).flatten('F')
#
#def fullxymeasallparas(pdeg,p1deg,p2deg,p4deg,z1p1deg,z1p2deg,z2p1deg,z2p2deg,z3p1deg,z3p2deg,s00,s01,s10,s11):
#    '''fullxymeas Qubit Map:(control,target)'''
#    return numpy.array([prob(op4(p,p1deg,p2deg,p4deg)@op3(z1p1deg,z1p2deg,z2p1deg,z2p2deg,z3p1deg,z3p2deg)@sini(s00,s01,s10,s11)) for p in pdeg])
#
#def fullxymeasallparas_flat(pdeg,p1deg,p2deg,p4deg,z1p1deg,z1p2deg,z2p1deg,z2p2deg,z3p1deg,z3p2deg,s00,s01,s10,s11):
#    return fullxymeasallparas(pdeg[0:int(len(pdeg)/4)],p1deg,p2deg,p4deg,z1p1deg,z1p2deg,z2p1deg,z2p2deg,z3p1deg,z3p2deg,s00,s01,s10,s11).flatten('F')
#
#def cnotfullxyfit(filename,z1p1deg=246,z1p2deg=147,z2p1deg=314,z2p2deg=234,z3p1deg=336,z3p2deg=37,s00=1,s01=0,s10=0,s11=0):
#    pattern_c='control(.*?)_'
#    pattern_t='target(.*?)_'
#    qubitid_c=re.search(pattern_c,filename).group(1)
#    qubitid_t=re.search(pattern_t,filename).group(1)
#    f=open(filename)
#    s=f.read().replace("'",'"')
#    f.close()
#    measphdeg=[]
#    Meas=[]
#    m00=[]
#    m01=[]
#    m10=[]
#    m11=[]
#    for l in s.split('\n')[:-1]:
#        m1=re.match("cnotphdeg (\d+) measphdeg (\d+) Sim (\{[\s\S]+\}) Meas (\{[\s\S]*\})",l)
#        m2=re.match("measphdeg( )(\d+) Sim (\{[\s\S]+\}) Meas (\{[\s\S]*\})",l)
#        if m1:
#            m=m1
#        if m2:
#            m=m2
#        if m:
#            g=m.groups()
#            meas=json.loads(g[3])
#            m00.append(meas['00'] if '00' in meas.keys() else 0)
#            m01.append(meas['01'] if '01' in meas.keys() else 0)
#            m10.append(meas['10'] if '10' in meas.keys() else 0)
#            m11.append(meas['11'] if '11' in meas.keys() else 0)
#            measphdeg.append(float(g[1]))
#        else:
#            print('not match',l)
#    xall=measphdeg
#    xfit=measphdeg[:-1]
#    z1p1deg=z1p1deg
#    z1p2deg=z1p2deg
#    z2p1deg=z2p1deg
#    z2p2deg=z2p2deg
#    z3p1deg=z3p1deg
#    z3p2deg=z3p2deg
#    s00=s00
#    s01=s01
#    s10=s10
#    s11=s11
#    xdata=numpy.tile(numpy.array(xfit),4)
#    '''TrueQ Meas Qubit Map (sorting qubit labels in ascending order)
#    Example: Q6Q5CNOT, Q6 for control, Q5 for target
#    TrueQ Meas Qubit Map: (Q5,Q6)
#    m00: Q5 |0>, Q6 |0>
#    m01: Q5 |0>, Q6 |1>
#    m10: Q5 |1>, Q6 |0>
#    m11: Q5 |1>, Q6 |1>
#    '''
#    if int(qubitid_c[1:])<int(qubitid_t[1:]):
#        mm00,mm01,mm10,mm11=m00,m01,m10,m11
#    else:
#        mm00,mm01,mm10,mm11=m00,m10,m01,m11
#    '''fullxymeas Qubit Map:(control,target)
#    mm00: control |0>, target |0>
#    mm01: control |0>, target |1>
#    mm10: control |1>, target |0>
#    mm11: control |1>, target |1>
#    '''
#    ydata=numpy.array([mm00[:-1],mm01[:-1],mm10[:-1],mm11[:-1]]).flatten()
#    popt,pcov=curve_fit(lambda pdeg,p1deg,p2deg,p4deg: fullxymeasallparas_flat(pdeg,p1deg,p2deg,p4deg,z1p1deg=z1p1deg,z1p2deg=z1p2deg,z2p1deg=z2p1deg,z2p2deg=z2p2deg,z3p1deg=z3p1deg,z3p2deg=z3p2deg,s00=s00,s01=s01,s10=s10,s11=s11),xdata,ydata,bounds=(([-360,-360,-360],[360,360,360])))#,p0=[0,0,0])
#    results=fullxymeasallparas_flat(xdata,*popt,z1p1deg=z1p1deg,z1p2deg=z1p2deg,z2p1deg=z2p1deg,z2p2deg=z2p2deg,z3p1deg=z3p1deg,z3p2deg=z3p2deg,s00=s00,s01=s01,s10=s10,s11=s11).reshape((4,-1))
#    print('popt',popt,'perr',numpy.sqrt(numpy.diag(pcov)))
#    print('results',results)
#
#    pyplot.figure()
#    pyplot.plot(xfit,results[0],label='Fit_c0t0',color='b',linestyle='--')
#    pyplot.plot(xfit,results[1],label='Fit_c0t1',color='orange',linestyle='--')
#    pyplot.plot(xfit,results[2],label='Fit_c1t0',color='g',linestyle='--')
#    pyplot.plot(xfit,results[3],label='Fit_c1t1',color='r',linestyle='--')
#    pyplot.plot(xall,mm00,label='Meas_c0t0',color='b',marker='*',linestyle='None')
#    pyplot.plot(xall,mm01,label='Meas_c0t1',color='orange',marker='*',linestyle='None')
#    pyplot.plot(xall,mm10,label='Meas_c1t0',color='g',marker='*',linestyle='None')
#    pyplot.plot(xall,mm11,label='Meas_c1t1',color='r',marker='*',linestyle='None')
#    pyplot.ylim([-0.1,1.1])
#    pyplot.axvline(0,linestyle='--',color='k')
#    pyplot.axvline(270,linestyle='--',color='k')
#    pyplot.axvline(360,linestyle='--',color='k')
#    pyplot.legend()
#    pyplot.xlabel('Rotation Angle (deg)')
#    pyplot.ylabel('Population')
#    pyplot.title('Full XY-Plane Measurement')
#    pyplot.savefig('cnotfullxyfit_'+filename[15:-4]+'.png')
#    pyplot.close()
#    return popt,numpy.sqrt(numpy.diag(pcov))
#
if __name__=="__main__":
    #    filename=sys.argv[1]
    #popt,perr=cnotfullxyfit(filename,z1p1deg=246,z1p2deg=147,z2p1deg=314,z2p2deg=234,z3p1deg=336,z3p2deg=37,s00=1,s01=0,s10=0,s11=0)
#    state=cnot()@zzdeg(246,147)@xx()@zzdeg(314,234)@xx()@zzdeg(336,37)@sini()
#    zz=prob(state)
#    xyscan=numpy.arange(0,numpy.pi*2,0.001)
#    xyval=numpy.zeros((len(xyscan),len(state)))
#    for i,p in enumerate(xyscan):
#        xyval[i,:]=tomoxy(state,p)[:,0]
#    from matplotlib import pyplot
#    pyplot.plot(xyscan,xyval)
#    pyplot.show()
    print(xxdeg(0,200))

