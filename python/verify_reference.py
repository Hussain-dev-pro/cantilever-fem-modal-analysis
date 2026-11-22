# Python reference for the FEM cantilever project (mirrors the MATLAB code)
import numpy as np, scipy.linalg as sl, json
from scipy.optimize import brentq
E=210e9; rho=7850; L=0.6; mt=1.2; zeta=0.02; g=9.81
n_op=955; f1x=n_op/60; EXC={'1x':(f1x,0.5),'2x':(2*f1x,1.5),'VP':(6*f1x,2.0)}   # Hz, base accel amplitude m/s^2
BAND=0.10; SEP=0.20; FAT=71; S_lim=FAT*(2e6/1e7)**(1/3)/2   # amplitude fatigue limit
def shs(B,t): return B*B-(B-2*t)**2,(B**4-(B-2*t)**4)/12,B/2
def beamFE(A,I,ne,mtip=mt,lumped=False):
    le=L/ne; nd=2*(ne+1); K=np.zeros((nd,nd)); M=np.zeros((nd,nd))
    ke=E*I/le**3*np.array([[12,6*le,-12,6*le],[6*le,4*le**2,-6*le,2*le**2],[-12,-6*le,12,-6*le],[6*le,2*le**2,-6*le,4*le**2]])
    if lumped: me=rho*A*le/2*np.diag([1,0,1,0])
    else: me=rho*A*le/420*np.array([[156,22*le,54,-13*le],[22*le,4*le**2,13*le,-3*le**2],[54,13*le,156,-22*le],[-13*le,-3*le**2,-22*le,4*le**2]])
    for e in range(ne):
        i=slice(2*e,2*e+4); K[i,i]+=ke; M[i,i]+=me
    M[-2,-2]+=mtip
    return K[2:,2:],M[2:,2:],ke
def modes(A,I,ne=40,mtip=mt,lumped=False):
    K,M,_=beamFE(A,I,ne,mtip,lumped)
    if lumped:   # Guyan-condense rotations (zero mass)
        w=np.arange(0,K.shape[0],2); t=np.arange(1,K.shape[0],2)
        Kc=K[np.ix_(w,w)]-K[np.ix_(w,t)]@np.linalg.solve(K[np.ix_(t,t)],K[np.ix_(t,w)])
        lam=sl.eigh(Kc,M[np.ix_(w,w)],eigvals_only=True); return np.sqrt(lam)/2/np.pi
    lam,Phi=sl.eigh(K,M); return np.sqrt(lam)/2/np.pi,Phi,K,M
def exact(A,I,mtip,n=4):
    a=mtip/(rho*A*L)
    F=lambda b: 1+np.cos(b)*np.cosh(b)+a*b*(np.cos(b)*np.sinh(b)-np.sin(b)*np.cosh(b))
    bs=[];x=np.linspace(0.05,16,20000);v=F(x)
    for i in range(len(x)-1):
        if v[i]*v[i+1]<0: bs.append(brentq(F,x[i],x[i+1]))
        if len(bs)==n: break
    return np.array(bs)**2*np.sqrt(E*I/(rho*A*L**4))/2/np.pi
def response(A,I,c,ne=40,scale=1.0):
    K,M,ke=beamFE(A,I,ne); lam,Phi=sl.eigh(K,M); w=np.sqrt(lam); r=np.zeros(len(K)); r[0::2]=1
    C=M@Phi@np.diag(2*zeta*w)@Phi.T@M; S=0;T=0;parts={}
    for k,(fx,a) in EXC.items():
        W=2*np.pi*fx*scale; u=np.linalg.solve(K-W**2*M+1j*W*C,-M@r*a)
        Fe=ke@np.r_[0,0,u[:2]]; s=abs(Fe[1])*c/I/1e6; tacc=abs(a-W**2*u[-2]); S+=s;T+=tacc; parts[k]=(s,tacc,a)
    return S,T,parts
def worst(A,I,c):
    sc=np.linspace(1-BAND,1+BAND,41); res=[response(A,I,c,scale=s) for s in sc]
    i=int(np.argmax([r[0] for r in res])); return res[i][0],max(r[1] for r in res),sc[i]
def robust_margin(A,I):
    m=1e9
    for es in (0.85,1.0,1.15):
        for ms in (0.8,1.0,1.2):
            global E
            E0=E; E=E0*es; f=modes(A,I,40,mt*ms)[0][:6]; E=E0
            m=min(m,margin(f))
    return m
def margin(f):
    m=1e9
    for fx,_ in EXC.values():
        lo,hi=fx*(1-BAND),fx*(1+BAND)
        for fi in f:
            if fi>300: continue
            d=0 if lo<=fi<=hi else (lo-fi)/lo if fi<lo else (fi-hi)/hi
            m=min(m,d)
    return m
R={}
# --- verification (original section) ---
A0,I0,c0=shs(.025,.0025)
R['S_lim']=S_lim
for tip,key in [(0,'notip'),(mt,'tip')]:
    fe=modes(A0,I0,40,tip)[0][:4]; ex=exact(A0,I0,tip); R['ver_'+key]=dict(fe=list(fe),ex=list(ex),err=list((fe-ex)/ex*100))
conv={'ne':[1,2,4,8,16,32,64]}
ex=exact(A0,I0,mt)
for lab,lum in [('cons',False),('lump',True)]:
    errs=[]
    for ne in conv['ne']:
        f=modes(A0,I0,ne,mt,lum); f=f[0] if not lum else f
        e=[(f[j]-ex[j])/ex[j]*100 if j<len(f) else np.nan for j in range(4)]; errs.append(e)
    conv[lab]=errs
R['conv']=conv
# static check
K,M,_=beamFE(A0,I0,40); Fv=np.zeros(len(K)); Fv[-2]=1.0; d=np.linalg.solve(K,Fv)[-2]; R['static']=dict(fe=d,exact=L**3/(3*E*I0))
# baseline
f,Phi,K,M=modes(A0,I0); r=np.zeros(len(K)); r[0::2]=1; Lm=Phi.T@M@r; meff=Lm**2/(r@M@r)
S,T,parts=response(A0,I0,c0); Sw,Tw,sw=worst(A0,I0,c0)
R['base']=dict(f=list(f[:4]),meff=list(meff[:4]*100),margin=margin(f[:6]),S=S,T=T,parts={k:list(v) for k,v in parts.items()},Sworst=Sw,Tworst=Tw,mass=rho*A0*L)
# --- catalogue redesign ---
cat=[]
for B in [20,25,30,35,40,50,60]:
    for t in [2,2.5,3,4,5]:
        if t>B/5: continue
        A,I,c=shs(B/1000,t/1000); f=modes(A,I)[0][:6]; mg=margin(f); Sw,Tw,_=worst(A,I,c)
        defl=(mt*g)*L**3/(3*E*I)+rho*A*g*L**4/(8*E*I)
        rm=robust_margin(A,I); ok=(mg>=SEP) and (rm>=0.10) and (Sw<=0.5*S_lim) and (t>=2.5) and (defl<=L/300)
        cat.append(dict(sec=f'SHS {B}x{t}',B=B,t=t,mass=rho*A*L,f1=f[0],f2=f[1],margin=mg*100,rmargin=rm*100,Sworst=Sw,Tworst=Tw,defl=defl*1e3,ok=ok))
best=min([c for c in cat if c['ok']],key=lambda c:(round(c['mass'],3),-c['margin']))
R['cat']=cat; R['best']=best
A1,I1,c1=shs(best['B']/1000,best['t']/1000); f,Phi,K,M=modes(A1,I1); Lm=Phi.T@M@r; meff=Lm**2/(r@M@r)
S,T,parts=response(A1,I1,c1); Sw,Tw,_=worst(A1,I1,c1)
R['final_rm']=robust_margin(A1,I1)
R['base_rm']=robust_margin(A0,I0)
R['final']=dict(f=list(f[:4]),meff=list(meff[:4]*100),margin=margin(f[:6]),S=S,T=T,parts={k:list(v) for k,v in parts.items()},Sworst=Sw,Tworst=Tw,mass=rho*A1*L)
# ruler rig
Er=200e9; br=.025; hr=.001; Lr=.25; mph=.19; Ar=br*hr; Ir=br*hr**3/12
def exr(mtip):
    a=mtip/(7850*Ar*Lr); F=lambda b:1+np.cos(b)*np.cosh(b)+a*b*(np.cos(b)*np.sinh(b)-np.sin(b)*np.cosh(b))
    x=np.linspace(.05,16,20000);v=F(x);bs=[brentq(F,x[i],x[i+1]) for i in range(len(x)-1) if v[i]*v[i+1]<0][:3]
    return list(np.array(bs)**2*np.sqrt(Er*Ir/(7850*Ar*Lr**4))/2/np.pi)
R['rig']=dict(notip=exr(0),phone=exr(mph))
json.dump(R,open('ref.json','w'),indent=1,default=float)
def p(x): return json.dumps(x,default=lambda o: round(float(o),3))
for k in ['S_lim','ver_notip','ver_tip','static','base','final','best','rig']: print(k,p(R[k]))
print('conv cons',p(R['conv']['cons'])); print('conv lump',p(R['conv']['lump']))
for c in cat: print(p(c))
