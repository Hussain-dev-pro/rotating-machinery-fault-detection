# Explanatory figures for the report (theory / design; real-data figures come from MATLAB screenshots)
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
NAVY='#1F3A5F'; BLUE='#2E86C1'; RED='#C0392B'; GRN='#1E8449'; ORG='#D68910'; GREY='#5D6D7E'
fs=12000; fr=1797/60; n=9; d=0.3126; D=1.537; r=d/D
BPFO=n/2*(1-r); BPFI=n/2*(1+r); BSF=D/(2*d)*(1-r**2); FTF=0.5*(1-r)

# 1 bearing geometry
fig,ax=plt.subplots(figsize=(9,4.2)); ax.set_aspect('equal'); ax.axis('off')
ax.add_patch(Circle((0,0),2.0,fc='#D6DBDF',ec=NAVY,lw=2)); ax.add_patch(Circle((0,0),1.55,fc='white',ec=NAVY,lw=1))
ax.add_patch(Circle((0,0),1.05,fc='#D6DBDF',ec=NAVY,lw=1)); ax.add_patch(Circle((0,0),0.65,fc='white',ec=NAVY,lw=2))
for k in range(9):
    a=2*np.pi*k/9; ax.add_patch(Circle((1.3*np.cos(a),1.3*np.sin(a)),0.24,fc='#85929E',ec=NAVY))
ax.plot([1.55*np.cos(-np.pi/2)],[1.55*np.sin(-np.pi/2)-0.02],'v',color=RED,ms=12)
ax.text(0,-2.35,'outer-race fault (6 o\'clock)',ha='center',color=RED,fontsize=9)
ax.annotate('',xy=(0.65,0),xytext=(0,0),arrowprops=dict(arrowstyle='->'));ax.text(0.1,0.1,'shaft',fontsize=8)
ax.annotate('d (ball)',xy=(1.3*np.cos(2*np.pi*2/9)+0.24,1.3*np.sin(2*np.pi*2/9)),xytext=(2.4,1.7),arrowprops=dict(arrowstyle='->'),fontsize=9)
ax.annotate('D (pitch circle)',xy=(-1.3,0),xytext=(-3.9,1.2),arrowprops=dict(arrowstyle='->'),fontsize=9)
ax.add_patch(Circle((0,0),1.3,fill=False,ls='--',ec=GREY))
txt=('SKF 6205-2RS (CWRU test bearing)\n n = 9 balls, d = 0.3126 in, D = 1.537 in\n\n'
     f'BPFO = n/2 (1 - d/D) fr = {BPFO:.4f} fr\nBPFI = n/2 (1 + d/D) fr = {BPFI:.4f} fr\n'
     f'BSF  = D/2d (1 - (d/D)^2) fr = {BSF:.4f} fr\nFTF  = 1/2 (1 - d/D) fr = {FTF:.4f} fr\n\n'
     f'At 1797 rpm (fr = {fr:.2f} Hz):\nBPFO {BPFO*fr:.1f} Hz, BPFI {BPFI*fr:.1f} Hz, BSF {BSF*fr:.1f} Hz')
ax.text(3.0,0,txt,va='center',fontsize=9.5,family='DejaVu Sans Mono')
ax.set_xlim(-4,8.3); ax.set_ylim(-2.6,2.3)
fig.savefig('figs/bearing_geometry.png',dpi=170,bbox_inches='tight'); plt.close()

# 2 why envelope: synthetic demonstration
rng=np.random.default_rng(3); t=np.arange(fs)/fs; rate=BPFO*fr
x=0.05*rng.standard_normal(fs)+0.03*np.sin(2*np.pi*fr*t)
for ti in np.arange(0,1,1/rate)+rng.normal(0,3e-4,len(np.arange(0,1,1/rate))):
    m=t>=ti; tt=t[m]-ti; x[m]+=0.3*np.exp(-800*tt)*np.sin(2*np.pi*3200*tt)
X=np.fft.fft(x); f=np.fft.fftfreq(fs,1/fs); mask=(f>2700)&(f<3700); z=np.fft.ifft(2*X*mask); e=np.abs(z)
w=np.hanning(fs); A=np.abs(np.fft.rfft(x*w))*2/w.sum(); E=np.abs(np.fft.rfft((e-e.mean())*w))*2/w.sum(); fp=np.fft.rfftfreq(fs,1/fs)
fig,axs=plt.subplots(2,2,figsize=(10,5.6))
axs[0,0].plot(t[:1200]*1e3,x[:1200],'k',lw=.6); axs[0,0].set_title('1  Raw signal: impacts hidden in noise',fontsize=10); axs[0,0].set_xlabel('time (ms)')
axs[0,1].plot(fp/1e3,A,color=GREY,lw=.6); axs[0,1].axvspan(2.7,3.7,color=BLUE,alpha=.15); axs[0,1].set_title('2  Raw FFT: energy at the 3.2 kHz resonance, no BPFO line',fontsize=10); axs[0,1].set_xlabel('frequency (kHz)')
axs[1,0].plot(t[:1200]*1e3,np.real(z)[:1200],color=BLUE,lw=.6); axs[1,0].plot(t[:1200]*1e3,e[:1200],color=RED,lw=1.2)
axs[1,0].set_title('3  Band-pass 2.7–3.7 kHz + envelope (red)',fontsize=10); axs[1,0].set_xlabel('time (ms)')
axs[1,1].plot(fp,E,color=BLUE,lw=.9); axs[1,1].set_xlim(0,500)
for h in (1,2,3,4): axs[1,1].axvline(h*rate,color=RED,ls=':',lw=1)
axs[1,1].set_title(f'4  Envelope spectrum: clear lines at BPFO = {rate:.0f} Hz and harmonics',fontsize=10); axs[1,1].set_xlabel('frequency (Hz)')
fig.tight_layout(); fig.savefig('figs/why_envelope.png',dpi=160); plt.close()

# 3 pipeline
fig,ax=plt.subplots(figsize=(10,2.6)); ax.axis('off'); ax.set_xlim(0,10); ax.set_ylim(0,2.4)
steps=[('Accelerometer\n12 kHz, 1 s','#D6EAF8'),('Spectral\nkurtosis band\n(70 candidates)','#D6EAF8'),('Hilbert\nenvelope\n(FFT method)','#D6EAF8'),
       ('Envelope\nspectrum\n(Hann)','#D6EAF8'),('Local-prominence\nindicators\nRo, Ri, Rb (dB)','#D5F5E3'),('Two-stage rule\nRMS alarm +\nsignature','#FDEBD0')]
for i,(s,c) in enumerate(steps):
    x0=0.1+i*1.65; ax.add_patch(FancyBboxPatch((x0,0.5),1.4,1.4,boxstyle='round,pad=0.05',fc=c,ec=NAVY))
    ax.text(x0+0.7,1.2,s,ha='center',va='center',fontsize=8.6)
    if i<5: ax.annotate('',xy=(x0+1.62,1.2),xytext=(x0+1.45,1.2),arrowprops=dict(arrowstyle='->',lw=1.4))
ax.text(5,0.15,'Healthy  /  Inner race  /  Ball  /  Outer race',ha='center',fontsize=10,color=NAVY,weight='bold')
fig.savefig('figs/pipeline.png',dpi=170,bbox_inches='tight'); plt.close()

# 4 ablation
lab=['v1  global floor\n(first version)','v2a local prominence\nno RMS stage','v2b local prominence\nRAW FFT + RMS','v2  local prominence\nENVELOPE + RMS']
tr=[60.8,76.2,79.6,90.8]; te=[53.3,63.7,74.6,87.5]
fig,ax=plt.subplots(figsize=(9,3.8)); xx=np.arange(4)
ax.bar(xx-0.18,tr,0.36,color='#AED6F1',label='train (0–1 hp)'); b=ax.bar(xx+0.18,te,0.36,color=NAVY,label='TEST (unseen 2–3 hp)')
for i,v in enumerate(te): ax.text(i+0.18,v+1.5,f'{v:.1f}%',ha='center',fontsize=9,weight='bold')
ax.axhline(25,color=GREY,ls=':',lw=1); ax.text(3.55,26.5,'random guess 25 %',fontsize=8,color=GREY,ha='right')
ax.set_xticks(xx); ax.set_xticklabels(lab,fontsize=8.5); ax.set_ylabel('balanced accuracy (%)'); ax.set_ylim(0,100); ax.legend(loc='upper left',fontsize=8.5,frameon=False)
fig.tight_layout(); fig.savefig('figs/ablation.png',dpi=170); plt.close()
print('ok', BPFO,BPFI,BSF,FTF)
