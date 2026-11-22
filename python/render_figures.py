import numpy as np, scipy.linalg as sl, json, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib import font_manager as fm
for f in ['Regular','Bold','SemiBold']: fm.fontManager.addfont(f'/usr/share/fonts/truetype/google-fonts/ofl/firasans/FiraSans-{f}.ttf')
plt.rcParams.update({'font.family':'Fira Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
exec(open('verify_reference.py').read().split('R={}')[0])
os.makedirs('figs',exist_ok=True); R=json.load(open('ref.json'))
NAVY='#1F3A5F';RED='#C0392B';GREEN='#1E8449'
A0,I0,c0=shs(.025,.0025); A1,I1,c1=shs(.035,.0025)
f0=modes(A0,I0)[0]; f1=modes(A1,I1)[0]
# 1 interference
fig,ax=plt.subplots(figsize=(7.2,3.6)); rpm=np.linspace(0,1300,10)
for k,(lab,mult) in enumerate([('1×',1),('2×',2),('Vane pass 6×',6)]): ax.plot(rpm,rpm/60*mult,lw=1.6,label=lab)
ax.axvspan(.9*955,1.1*955,color='#FDEBD0',label='Operating band ±10%')
ax.axhline(f0[0],color=RED,lw=2);ax.text(30,f0[0]+2,f'Original f₁ = {f0[0]:.1f} Hz',color=RED)
ax.axhline(f1[0],color=GREEN,lw=2);ax.text(30,f1[0]+2,f'Redesign f₁ = {f1[0]:.1f} Hz',color=GREEN)
ax.plot([955*2/2],[2*955/60],'o',ms=10,mfc='none',mec=RED,mew=2)
ax.set_xlim(0,1300);ax.set_ylim(0,120);ax.set_xlabel('Pump speed (rpm)');ax.set_ylabel('Frequency (Hz)');ax.legend(loc='upper left',fontsize=8.5);ax.grid(alpha=.3)
ax.set_title('Interference diagram: the original bracket\'s first mode sits on the 2× line at 955 rpm',fontsize=9.5,loc='left')
fig.tight_layout();fig.savefig('figs/interference.png',dpi=220);plt.close()
# 2 convergence
ne=np.array(R['conv']['ne']);C=np.array(R['conv']['cons'],float);Lm=np.array(R['conv']['lump'],float)
fig,ax=plt.subplots(figsize=(7.2,3.6));cols=['#1F77B4','#FF7F0E','#2CA02C','#D62728']
for j in range(4):
    ax.loglog(ne,np.abs(C[:,j]),'-o',color=cols[j],ms=4,label=f'Mode {j+1} consistent')
    ax.loglog(ne,np.abs(Lm[:,j]),'--s',color=cols[j],ms=3.5,alpha=.7,label=f'Mode {j+1} lumped')
x=np.array([4,64]);ax.loglog(x,2e0*(x/4)**-4,'k:',lw=1);ax.text(20,2e-3,'slope −4',fontsize=8.5);ax.loglog(x,1.0*(x/4)**-2,'k:',lw=1);ax.text(24,3.5e-2,'slope −2',fontsize=8.5)
ax.set_ylim(1e-9,1e2);ax.set_xlabel('Number of elements');ax.set_ylabel('|Error vs exact| (%)');ax.grid(alpha=.3,which='both');ax.legend(fontsize=7.5,ncol=2,loc='lower left')
ax.set_title('Mesh convergence against the exact solution (consistent ~h⁴, lumped ~h²)',fontsize=9.5,loc='left')
fig.tight_layout();fig.savefig('figs/convergence.png',dpi=220);plt.close()
# 3 mode shapes
f,Phi,K,M=modes(A0,I0); x=np.linspace(0,600,41); r=np.zeros(len(K)); r[0::2]=1; meff=(Phi.T@M@r)**2/(r@M@r)*100
fig,axs=plt.subplots(2,2,figsize=(7.2,4.4))
for j,ax in enumerate(axs.flat):
    w=np.r_[0,Phi[0::2,j]]; w=w/np.max(np.abs(w))*np.sign(w[-1] if abs(w[-1])>.05 else 1)
    ax.axhline(0,color='#999',lw=.8);ax.plot(x,w,color=NAVY,lw=2.2);ax.plot(600,w[-1],'s',ms=8,mfc='#E67E22',mec='k')
    ax.add_patch(plt.Rectangle((-30,-1.3),30,2.6,color='#AAB7B8'));ax.set_xlim(-35,640);ax.set_ylim(-1.3,1.3)
    ax.set_title(f'Mode {j+1}: {f[j]:.1f} Hz  (effective mass {meff[j]:.1f}%)',fontsize=9.5)
    if j>1: ax.set_xlabel('Position along bracket (mm)')
fig.tight_layout();fig.savefig('figs/mode_shapes.png',dpi=220);plt.close()
# 4 FRF compare
def frf(A,I,c,fh):
    K,M,ke=beamFE(A,I,40); lam,Phi=sl.eigh(K,M); w=np.sqrt(lam); r=np.zeros(len(K)); r[0::2]=1
    C=M@Phi@np.diag(2*zeta*w)@Phi.T@M; T=[];S=[]
    for fx in fh:
        W=2*np.pi*fx; u=np.linalg.solve(K-W**2*M+1j*W*C,-M@r); T.append(abs(1-W**2*u[-2])); S.append(abs((ke@np.r_[0,0,u[:2]])[1])*c/I/1e6)
    return np.array(T),np.array(S)
fh=np.linspace(2,150,1500); T0,S0=frf(A0,I0,c0,fh); T1,S1=frf(A1,I1,c1,fh)
fig,axs=plt.subplots(2,1,figsize=(7.2,4.6),sharex=True)
for ax,a,b,yl in [(axs[0],T0,T1,'Tip / base acceleration'),(axs[1],S0,S1,'Root stress per 1 m/s² (MPa)')]:
    ax.semilogy(fh,a,color=RED,lw=1.6,label='Original SHS 25×2.5');ax.semilogy(fh,b,color=GREEN,lw=1.6,label='Redesign SHS 35×2.5')
    for k,(fx,_) in EXC.items(): ax.axvspan(fx*.9,fx*1.1,color='#FDEBD0',alpha=.8)
    ax.set_ylabel(yl);ax.grid(alpha=.3,which='both')
axs[0].legend(fontsize=8.5);axs[1].set_xlabel('Excitation frequency (Hz)');axs[0].set_title('Response to skid vibration (shaded = 1×, 2× and vane-pass bands)',fontsize=9.5,loc='left')
fig.tight_layout();fig.savefig('figs/frf_compare.png',dpi=220);plt.close()
# 5 design space
cat=R['cat'];fig,ax=plt.subplots(figsize=(7.2,3.6))
for c in cat:
    ok=bool(c['ok']);ax.plot(c['mass'],c['rmargin'],'o' if ok else 'x',color=GREEN if ok else RED,ms=7 if ok else 6,mew=1.5)
for c in cat:
    if c['sec'] in ('SHS 25x2.5','SHS 20x2.5','SHS 35x2.5'): ax.annotate(c['sec'].replace('x','×'),(c['mass'],c['rmargin']),xytext=(8,-4 if c['sec']!='SHS 35x2.5' else 6),textcoords='offset points',fontsize=9,weight='bold')
b=R['best'];ax.plot(b['mass'],b['rmargin'],'o',ms=15,mfc='none',mec='k',mew=2)
ax.axhline(10,color='#444',ls='--',lw=1);ax.text(4.4,10.6,'robust margin 10%',fontsize=8.5)
ax.set_xlabel('Bracket mass (kg)');ax.set_ylabel('Robust separation margin (%)');ax.grid(alpha=.3)
ax.plot([],[],'o',color=GREEN,label='passes all checks');ax.plot([],[],'x',color=RED,label='fails a check');ax.legend(fontsize=8.5,loc='upper right')
ax.set_title('Design space: 35 square hollow sections checked; lightest passing = SHS 35×2.5',fontsize=9.5,loc='left')
fig.tight_layout();fig.savefig('figs/design_space.png',dpi=220);plt.close()
# 6 animations (twins)
def forced(A,I,c,W,a):
    K,M,ke=beamFE(A,I,40); lam,Phi=sl.eigh(K,M); w=np.sqrt(lam); r=np.zeros(len(K)); r[0::2]=1
    C=M@Phi@np.diag(2*zeta*w)@Phi.T@M; u=np.linalg.solve(K-W**2*M+1j*W*C,-M@r*a)
    return np.r_[0,u[0::2]],abs((ke@np.r_[0,0,u[:2]])[1])*c/I/1e6,abs(a-W**2*u[-2])
W=2*np.pi*EXC['2x'][0];a=EXC['2x'][1];w0,s0,t0=forced(A0,I0,c0,W,a);w1,s1,t1=forced(A1,I1,c1,W,a)
fig,axs=plt.subplots(1,2,figsize=(12,5.2),dpi=80);mag=60;H=[]
for ax,w,s,t,ttl,col in [(axs[0],w0,s0,t0,'ORIGINAL  SHS 25×2.5',RED),(axs[1],w1,s1,t1,'REDESIGN  SHS 35×2.5',GREEN)]:
    ax.add_patch(plt.Rectangle((-35,-70),35,140,color='#AAB7B8'));env=np.abs(w)*mag*1e3;ax.plot(x,env,':',color='#bbb');ax.plot(x,-env,':',color='#bbb')
    l,=ax.plot(x,0*x,color=NAVY,lw=5);m,=ax.plot(600,0,'s',ms=14,mfc='#E67E22',mec='k');H.append((l,m,w))
    ax.set_xlim(-45,650);ax.set_ylim(-70,70);ax.set_title(ttl,color=col,fontsize=15,weight='bold');ax.set_xlabel('Position (mm)');ax.set_ylabel(f'Deflection ×{mag} (mm)');ax.grid(alpha=.3)
    ax.text(320,-58,f'Tip acceleration {t:.1f} m/s² (base {a})\nRoot stress {s:.1f} MPa (fatigue limit {S_lim:.1f})',ha='center',color=col,weight='bold')
fig.text(.06,.02,f'Skid shakes at 2× running speed = {EXC["2x"][0]:.1f} Hz. Original f₁ = {f0[0]:.1f} Hz → resonance; redesign f₁ = {f1[0]:.1f} Hz. Same magnification, slowed down. (Python twin of FEM_Animation.m)',fontsize=10)
fig.tight_layout(rect=(0,.05,1,1))
def up(n):
    th=2*np.pi*n/20
    for l,m,w in H: y=np.real(w*np.exp(1j*th))*mag*1e3; l.set_ydata(y); m.set_ydata([y[-1]])
FuncAnimation(fig,up,frames=60).save('FEM_forced.mp4',writer=FFMpegWriter(fps=20,bitrate=1500));plt.close()
fig,axs=plt.subplots(2,2,figsize=(10,5.6),dpi=80);HL=[]
for j,ax in enumerate(axs.flat):
    wj=np.r_[0,Phi[0::2,j]];wj=wj/np.max(np.abs(wj));ax.add_patch(plt.Rectangle((-30,-1.4),30,2.8,color='#AAB7B8'))
    ax.plot(x,wj,':',color='#bbb');ax.plot(x,-wj,':',color='#bbb');l,=ax.plot(x,0*x,color=NAVY,lw=3.5);m,=ax.plot(600,0,'s',ms=10,mfc='#E67E22',mec='k');HL.append((l,m,wj))
    ax.set_xlim(-35,640);ax.set_ylim(-1.4,1.4);ax.set_title(f'Mode {j+1}: {f[j]:.1f} Hz',weight='bold')
fig.suptitle('Original bracket SHS 25×2.5 + 1.2 kg box: first four mode shapes (Python twin of FEM_Animation.m)',fontsize=11);fig.tight_layout()
def up2(n):
    s=np.sin(2*np.pi*n/20)
    for l,m,w in HL: l.set_ydata(w*s); m.set_ydata([w[-1]*s])
FuncAnimation(fig,up2,frames=40).save('FEM_modes.mp4',writer=FFMpegWriter(fps=20,bitrate=1500));plt.close()
print('ok')
