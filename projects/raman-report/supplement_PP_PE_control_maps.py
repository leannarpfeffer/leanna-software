"""Reproduce supplement_PP_PE_control_maps. Uses included data; see README.txt."""
from plot_common import *

def main():
    fig,(ax,bx)=plt.subplots(1,2,figsize=(6,2.85))
    fig.subplots_adjust(left=.10,right=.99,top=.84,bottom=.27,wspace=.27)
    for a,name,title in [(ax,'PP-maps','a  Neat PP'),(bx,'PE-map','b  Neat PE')]:
     d=D[name];im=mapplot(a,d)
     a.set_title(f"{title}\nQ pooled = {S.loc[name,'bounded_pool']:.3f}",pad=5,fontsize=9.2,linespacing=1.2)
    cbax=fig.add_axes([.31,.065,.43,.027]);cb=fig.colorbar(im,cax=cbax,orientation='horizontal',ticks=np.linspace(0,1,6))
    cb.ax.tick_params(labelsize=8,length=2,pad=2);cb.set_label('PE-associated index Q (uncalibrated)',fontsize=9,labelpad=2)
    save(fig,'Supplement_PP_PE_Control_Maps')

if __name__ == '__main__':
    main()
