"""Reproduce supplement_Raman_workflow. Uses included data; see README.txt."""
from plot_common import *

def main():
    median_width=float(E4_TRANSECTS.width_10_90_um.median())
    # Completed and planned branches remain explicit, with color secondary to line style.
    fig,ax=plt.subplots(figsize=(6,4.25));ax.set_xlim(0,10);ax.set_ylim(0,7);ax.axis('off')
    fig.subplots_adjust(left=.01,right=.99,bottom=.02,top=.99)
    def box(x,y,w,h,title,body,planned=False):
     p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.08,rounding_size=.04',facecolor='#FFF9F1' if planned else '#F0F7FA',edgecolor=PE if planned else PP,linestyle='--' if planned else '-',linewidth=1.0);ax.add_patch(p)
     ax.text(x+w/2,y+h-.17,title,ha='center',va='top',fontweight='bold',fontsize=9.0,linespacing=1.08)
     ax.text(x+w/2,y+h-.60,body,ha='center',va='top',fontsize=8.6,linespacing=1.2)
    def arrow(a,b,dash=False):
     ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=9,color=PE if dash else PP,linestyle='--' if dash else '-',linewidth=.9))
    ax.text(2.35,6.75,'COMPLETED',ha='center',fontsize=10,fontweight='bold')
    ax.text(7.63,6.75,'PLANNED',ha='center',fontsize=10,fontweight='bold')
    box(.15,4.73,4.4,1.55,'Multiscale Raman maps','E0–E5: 2, 1 and 0.1 µm steps\nNeat PP and PE: 2 µm controls')
    box(.15,2.55,4.4,1.55,'Chemical and boundary analysis','Native maps and matched fields\nE4 transition width; regional spectra')
    arrow((2.35,4.70),(2.35,4.21))
    box(5.42,4.73,4.4,1.55,'Prepare boundary targets','Register Raman-selected regions\nVerify sections and fiducials',True)
    box(5.42,2.55,4.4,1.55,'STEM–EELS measurements','Establish dose limits with controls\nAcquire spectra across regions',True)
    box(5.42,.38,4.4,1.55,'Cross-method comparison','Measure spatial response\nTest width and structural hypotheses',True)
    arrow((4.65,3.30),(5.30,5.41),True);arrow((7.62,4.70),(7.62,4.21),True);arrow((7.62,2.52),(7.62,2.03),True)
    ax.text(.15,1.70,'Established by Raman',fontweight='bold',fontsize=9.2)
    ax.text(.15,1.30,f'Microscale chemical heterogeneity\nE4 apparent width: {median_width:.2f} µm median',fontsize=8.8,va='top',linespacing=1.2)
    ax.text(.15,.55,'Intrinsic interphase thickness and\nlocal crystallinity are undetermined.',fontsize=8.5,va='top',linespacing=1.2)
    save(fig,'Supplement_Raman_Workflow')
    

if __name__ == '__main__':
    main()
