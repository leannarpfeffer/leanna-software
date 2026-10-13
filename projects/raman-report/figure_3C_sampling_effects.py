"""Reproduce figure_3C_sampling_effects. Uses included data; see README.txt."""
from plot_common import *

def main():
    # Sampling results: stable data with visual distinctions in hue and marker.
    fig,(ax,bx)=plt.subplots(1,2,figsize=(6,3.80))
    fig.subplots_adjust(left=.09,right=.99,bottom=.26,top=.90,wspace=.40)
    y=np.arange(6)
    for col,delta,marker,face,title in [('r1_coarse_full',-.16,'s',PP,'Full 1 µm'),('r1_coarse_roi',0,'o',GRAY,'Matched 1 µm'),('r1_fine',.16,'D',PE,'Native 0.1 µm')]:
     ax.scatter(MATCH.loc[[f'E{i}' for i in range(6)],col],y+delta,marker=marker,s=24,edgecolor=face,facecolor=face,label=title,zorder=4)
    ax.set_xscale('log');ax.set_xlim(.07,1.7);ax.set_xticks([.1,.2,.5,1]);ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_yticks(y,[f'E{i}' for i in range(6)]);ax.invert_yaxis();ax.set_xlabel('Pooled $A_{1295}/A_{1330}$\n(dimensionless)',labelpad=4)
    ax.set_title('a  Full and matched fields',loc='left',pad=8);clean(ax)
    ax.legend(loc='upper left',bbox_to_anchor=(-.06,-.24),frameon=False,fontsize=8.3,ncol=1,handletextpad=.5,labelspacing=.25)
    for i in range(6):
     fine=MATCH.loc[f'E{i}','r1_fine']
     for st,dy,c,marker in [(1,-.13,PP,'o'),(2,.13,PE,'s')]:
      q=OFF[(OFF['sample']==f'E{i}')&(OFF.step_um==st)].r1.quantile([.05,.5,.95]).to_numpy();rel=100*(q/fine-1)
      bx.plot(rel[[0,2]],[i+dy]*2,color=c,lw=1.7)
      bx.plot(rel[1],i+dy,marker,color=c,markersize=3,label=f'{st} µm grid' if i==0 else None)
    bx.axvline(0,color='.6',ls='--',lw=.8);bx.set_xlim(-40,45);bx.set_xticks([-40,-20,0,20,40]);bx.set_yticks(y,[f'E{i}' for i in range(6)]);bx.invert_yaxis()
    bx.set_xlabel('Change from native pooled index (%)',labelpad=4);bx.set_title('b  Coarser-grid placement',loc='left',pad=8);clean(bx)
    bx.legend(frameon=False,loc='upper left',bbox_to_anchor=(-.05,-.24),fontsize=8.3,labelspacing=.25)
    save(fig,'Figure_3C_Sampling_Effects')
    

if __name__ == '__main__':
    main()
