"""Reproduce figure_2_regional_spectra. Uses included data; see README.txt."""
from plot_common import *

def main():
    # Spectra and structural ratios: fixed region colors throughout.
    fig,axes=plt.subplots(2,2,figsize=(6,5.55))
    fig.subplots_adjust(left=.11,right=.99,top=.96,bottom=.10,wspace=.42,hspace=.53)
    ax,bx,cx,dx=axes.flat
    for name,col,ls,txt in [('PP-maps',PP,'-','Neat PP'),('PE-map',PE,'--','Neat PE')]:
     d=D[name];sel=(d['x']>=740)&(d['x']<=1490);v=d['corrected']/d['corrected'][sel].max();ax.plot(d['x'][sel],v[sel],color=col,ls=ls,lw=1.0,label=txt)
    for center,lab in [(809,'809'),(974,'974'),(1128,'1128'),(1295,'1295'),(1416,'1416')]:
     ax.axvline(center,color='.87',lw=.5);ax.text(center,1.11,lab,ha='center',fontsize=7)
    ax.set_xlim(740,1490);ax.set_xticks([800,1000,1200,1400]);ax.set_ylim(-.04,1.48)
    ax.set_xlabel(r'Raman shift (cm$^{-1}$)');ax.set_ylabel('Normalized intensity')
    ax.set_title('a  Neat-polymer references',loc='left',pad=7)
    ax.legend(frameon=False,fontsize=8,loc='upper right',handlelength=1.6)
    for i in range(6):
     q=REG[REG['sample']==f'E{i}'].set_index('region');lo=q.loc['lower_PE_marker','PP809_841_shared_baseline'];hi=q.loc['higher_PE_marker','PP809_841_shared_baseline']
     bx.plot([i,i],[lo,hi],color='.6',lw=.8)
     bx.scatter(i,lo,s=25,facecolor='white',edgecolor=PP,label='Lower quartile' if i==0 else None,zorder=4)
     bx.scatter(i,hi,s=25,color=PE,marker='s',label='Higher quartile' if i==0 else None,zorder=4)
    bx.set_xticks(range(6),[f'E{i}' for i in range(6)]);bx.set_ylabel('Pooled $A_{809}/A_{841}$');bx.set_ylim(.92,1.64)
    bx.set_title('b  Marker-defined regions',loc='left',pad=7)
    bx.legend(frameon=False,fontsize=8,loc='upper right',labelspacing=.2);bx.grid(axis='y',color='.92')
    spec=E4_SPECTRA.copy();x=spec.wavenumber_cm_1.to_numpy()
    ppbase={'integration':[790,860],'left':[765,785],'right':[865,885]}
    for col,color,ls,lab in [('PP_side',PP,'-','PP-rich side'),('interface',INTERFACE,'--','Interface'),('PE_side',PE,':','PE-rich side')]:
     v=spec[col].to_numpy();corr=local_corrected(v[None,:],x,ppbase)[0];area=v@band_weights(x,dict(ppbase,integration=[833,849]))[0];sel=(x>=785)&(x<=868)
     cx.plot(x[sel],corr[sel]/area,color=color,ls=ls,lw=1.25,label=lab)
    cx.set_xlabel(r'Raman shift (cm$^{-1}$)');cx.set_ylabel('Intensity / integrated $A_{841}$')
    cx.set_title('c  E4 regional PP bands',loc='left',pad=7)
    cx.set_ylim(-.008,.215)
    cx.legend(frameon=False,fontsize=7.7,loc='upper right',labelspacing=.2,handlelength=1.5)
    sel=(x>=1260)&(x<=1460)
    dx.plot(x[sel],spec.interface[sel],color=INTERFACE,lw=1.25,label='Observed interface')
    dx.plot(x[sel],spec.interface_mixture_prediction[sel],color=GRAY,ls='--',lw=1.15,label='Mixture prediction')
    dx.set_xlabel(r'Raman shift (cm$^{-1}$)');dx.set_ylabel('Mean intensity (counts)')
    dx.set_title('d  E4 held-out PE markers',loc='left',pad=7)
    dx.legend(frameon=False,fontsize=7.6,loc='upper right',labelspacing=.2,handlelength=1.5);dx.set_ylim(5,67)
    save(fig,'Figure_2_Regional_Spectra')
    

if __name__ == '__main__':
    main()
