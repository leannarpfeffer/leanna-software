"""Reproduce figure_3A_maps_E0_E2. Uses included data; see README.txt."""
from plot_common import *

def main():
    nominal=[.5,1,5,10,20,50]
    coarse=['E0-map2','E1-map','E2-map1','E3-map','E4-map','E5-map']
    for group,idx in [('A',range(3))]:
     fig,axes=plt.subplots(3,3,figsize=(6,6.10))
     fig.subplots_adjust(left=.09,right=.99,top=.91,bottom=.13,wspace=.40,hspace=.54)
     for j,t in enumerate(['2 µm step\n50× objective','1 µm step\n100× objective','0.1 µm step\n100× objective']):
      fig.text([.207,.522,.837][j],.992,t,ha='center',va='top',fontsize=9.0,fontweight='bold',linespacing=1.13)
     for row,i in enumerate(idx):
      for j,name in enumerate([coarse[i],f'E{i}-2-1.0',f'E{i}-small']):
       d=D[name];ax=axes[row,j];im=mapplot(ax,d)
       ax.set_title(f"E{i}: {d['width']:.0f} × {d['height']:.0f} µm",pad=4,fontsize=9.1)
       if j==1:
        f=D[f'E{i}-small'];xy=(f['xmin']-d['xmin'],f['ymin']-d['ymin'])
        rr=Rectangle(xy,f['width'],f['height'],fill=False,edgecolor='white',linestyle='--',linewidth=.9)
        rr.set_path_effects([patheffects.Stroke(linewidth=1.5,foreground='black'),patheffects.Normal()]);ax.add_patch(rr)
     cbax=fig.add_axes([.25,.060,.52,.017]);cb=fig.colorbar(im,cax=cbax,orientation='horizontal',ticks=np.linspace(0,1,6))
     cb.ax.tick_params(labelsize=8,length=2,pad=2)
     cb.set_label(r'PE-associated index $Q = A_{1295}^{+}/(A_{1295}^{+}+A_{809}^{+})$',labelpad=2,fontsize=9)
     save(fig,f'Figure_3{group}_Native_Raman_Maps')
    

if __name__ == '__main__':
    main()
