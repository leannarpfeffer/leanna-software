"""Reproduce figure_3D_E4_interface_width. Uses included data; see README.txt."""
from plot_common import *

def main():
    # E4 interface: keep seven profiles distinct and matching markers in widths.
    tr=E4_TRANSECTS.copy();pr=E4_PROFILES.copy()
    median_width=float(tr.width_10_90_um.median())
    colors=plt.get_cmap('viridis')(np.linspace(.12,.86,len(tr)))
    color_by_transect=dict(zip(tr.transect,colors))
    fig=plt.figure(figsize=(6,5.90))
    gs=fig.add_gridspec(2,2,left=.09,right=.99,top=.94,bottom=.10,wspace=.40,hspace=.55,height_ratios=[1,1])
    ax=fig.add_subplot(gs[0,0]);bx=fig.add_subplot(gs[0,1]);cx=fig.add_subplot(gs[1,:])
    im=mapplot(ax,D['E4-small']);ax.set_title('a  Boundary-normal transects',loc='left',pad=6)
    for row in tr.itertuples():
     t=np.array([-2.5,2.5]);xx=row.seed_x_um+t*row.normal_x;yy=row.seed_y_um+t*row.normal_y
     line=ax.plot(xx,yy,color='white',lw=.75)[0];line.set_path_effects([patheffects.Stroke(linewidth=1.3,foreground='black'),patheffects.Normal()])
     ax.text(xx[0]-.05,yy[0],str(row.transect),color='white',fontsize=7.8,ha='right',va='center',path_effects=[patheffects.withStroke(linewidth=1.3,foreground='black')])
    for k,g in pr.groupby('transect'):
     bx.plot(g.distance_from_fitted_midpoint_um,g.normalized_local_contrast,color=color_by_transect[k],lw=.8,alpha=.9)
    t=np.linspace(-2.5,2.5,200);sigma=np.median(tr.sigma_um)
    bx.plot(t,.5*erfc(t/(np.sqrt(2)*sigma)),color='black',lw=1.7,label='Median fitted width')
    bx.axhline(.1,color='.7',lw=.6,ls=':');bx.axhline(.9,color='.7',lw=.6,ls=':')
    bx.annotate('',xy=(median_width/2,.5),xytext=(-median_width/2,.5),arrowprops={'arrowstyle':'<->','color':'black','lw':.9})
    bx.text(0,.56,f'{median_width:.2f} µm',ha='center',fontsize=8.5)
    bx.set_xlim(-2.5,2.5);bx.set_xticks([-2,0,2]);bx.set_ylim(-.08,1.13)
    bx.set_xlabel('Distance from fitted midpoint (µm)');bx.set_ylabel('Normalized local contrast')
    profile_count='Seven' if len(tr)==7 else str(len(tr))
    bx.set_title(f'b  {profile_count} boundary profiles',loc='left',pad=6);bx.legend(frameon=False,fontsize=7.6,loc='upper right')
    cx.axhspan(tr.width_10_90_um.min(),tr.width_10_90_um.max(),color='#EDF5F8')
    cx.axhline(median_width,color=PP,ls='--',lw=.9,label=f'Median {median_width:.2f} µm')
    cx.scatter(tr.transect-.06,tr.width_10_90_um,color=colors,s=27,label='Five-line strip',zorder=3)
    cx.scatter(tr.transect+.06,tr.width_single_line_um,facecolor='white',edgecolor=colors,s=27,label='Single line',zorder=3)
    width_min=min(tr.width_10_90_um.min(),tr.width_single_line_um.min())
    width_max=max(tr.width_10_90_um.max(),tr.width_single_line_um.max())
    cx.set_ylim(min(1.1,width_min-.05),max(1.8,width_max+.05))
    cx.set_xticks(tr.transect);cx.set_xlim(tr.transect.min()-.5,tr.transect.max()+.5)
    cx.set_xlabel('Transect along the same E4 boundary');cx.set_ylabel('10–90% width (µm)')
    cx.set_title('c  Spatial variation in apparent width',loc='left',pad=6)
    cx.legend(frameon=False,fontsize=8.1,loc='lower center',ncol=3,bbox_to_anchor=(.5,-.45),handletextpad=.5,columnspacing=1)
    save(fig,'Figure_3D_E4_Interface_Width')
    

if __name__ == '__main__':
    main()
