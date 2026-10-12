from pathlib import Path
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Circle, Patch, Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/AME_local_revision_20261011'
FIG = OUT / 'Figures'
FIG.mkdir(parents=True, exist_ok=True)
DATA = ROOT / 'outputs/visual_upgrade_biorender_data_driven_v1_polished/derived_plot_data'
SUPP = ROOT / 'outputs/bmc_ai_second_stage_performance_package/additional_file_2_extracted'
SD = ROOT / 'source_data'
ORDER = ['A131','A129','A067','A097','A168','A001','A059','A098','A065','A173']
BLUE, TEAL, RED, GRAY, GOLD = '#18588A','#368777','#B44949','#DCE1E6','#C99A32'
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.titlesize':11,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','axes.labelcolor':'#222222','text.color':'#222222','legend.frameon':False})

def title(ax, letter, label):
    ax.set_title(label,loc='left',pad=14)
    ax.text(-.13,1.04,letter,transform=ax.transAxes,fontsize=15,fontweight='bold',va='bottom')

def save(fig,n):
    fig.savefig(FIG/f'figure {n}-revised.png',dpi=400,bbox_inches='tight',facecolor='white')
    fig.savefig(FIG/f'figure {n}-revised.pdf',bbox_inches='tight',facecolor='white')
    plt.close(fig)

def f1():
    fig=plt.figure(figsize=(9.5,6.9))
    gs=fig.add_gridspec(2,2,width_ratios=[1,1.5],height_ratios=[1,1.6],wspace=.75,hspace=.8)
    a=fig.add_subplot(gs[:,0]); a.axis('off'); title(a,'A','Audit framework')
    labels=['10,061\nBefore cross-database deduplication','7,262\nUnique verification records','178\nFull-text audit corpus','57\nReproducible formulas','27\nFurther evaluation','10\nPrimary candidates','1\nFixed-cutoff candidate']
    for i,lab in enumerate(labels):
        y=.96-i*.145
        a.text(.5,y,lab,ha='center',va='center',fontsize=9,bbox=dict(boxstyle='round,pad=.5',facecolor='white',edgecolor=BLUE if i<2 else TEAL,linewidth=1.2),transform=a.transAxes)
        if i<6:
            a.annotate('',xy=(.5,y-.096),xytext=(.5,y-.048),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color='#777777',linestyle='--' if i==1 else '-'))
    a.text(.5,-.04,'Dashed link: verification records\nand the completed audit corpus',ha='center',va='top',fontsize=8,transform=a.transAxes)
    b=fig.add_subplot(gs[0,1]);title(b,'B','Record and signature counts')
    b.barh([1,0],[2799,7084],color=[GRAY,GOLD],height=.5)
    b.set_yticks([1,0],['Cross-database\nduplicates','Difference before\nfull-text audit']);b.set_xlim(0,8200)
    b.set_xlabel('Records');b.set_xticks([0,2000,4000,6000,8000])
    for y,v in zip([1,0],[2799,7084]):b.text(v+110,y,f'{v:,}',va='center')
    c=fig.add_subplot(gs[1,1]); vals=[121,30,17,9]; yy=np.arange(4)[::-1]
    c.barh(yy,vals,color=[RED,GRAY,GRAY,RED],height=.5)
    c.set_yticks(yy,['Non-reproducible','Not advanced','Not primary','Not fixed-cutoff\ncandidate']);c.set_xlim(0,150);c.set_xlabel('Articles/signatures')
    for y,v in zip(yy,vals):c.text(v+2,y,str(v),va='center')
    fig.subplots_adjust(left=.09,right=.97,top=.94,bottom=.15)
    fig.text(.51,.03,'7,084 = 6,028 + 498 + 276 + 282\nComponent definitions and record linkage are given in the legend.',fontsize=8,ha='left')
    save(fig,1)

def f2():
    fig=plt.figure(figsize=(10,9.2));gs=fig.add_gridspec(2,2,width_ratios=[1,1.3],height_ratios=[1,1.15],wspace=.72,hspace=.65)
    a=fig.add_subplot(gs[0,0]);title(a,'A','Reporting completeness')
    d=pd.read_csv(DATA/'figure_2_scores_data.csv'); groups=list(d.group.unique()); vals=[d[d.group==g].score_0_10.dropna() for g in groups]
    vp=a.violinplot(vals,showmedians=True);rng=np.random.default_rng(3)
    for i,v in enumerate(vals):
        vp['bodies'][i].set_facecolor([TEAL,GRAY][i]);vp['bodies'][i].set_alpha(.6)
        a.scatter(rng.normal(i+1,.035,len(v)),v,s=8,color='#666666',alpha=.4)
    a.set_xticks([1,2],['Reproducible\nn=57','L0\nn=121']);a.set_ylabel('Reporting score (0-10)');a.set_ylim(-.3,10.5)
    c=fig.add_subplot(gs[1,0]);title(c,'C','Deployment elements')
    counts=np.array([12,0,0,47]); yy=np.arange(4)[::-1];pct=counts/57*100
    c.hlines(yy,0,100,color=GRAY);c.scatter(pct,yy,c=TEAL,s=50,zorder=3)
    c.set_yticks(yy,['Numeric cutoff','Code/calculator','Intended clinical use','Preprocessing details'])
    for y,n,p in zip(yy,counts,pct):c.text(103,y,f'{n}/57 ({p:.1f}%)',va='center',fontsize=8)
    c.set_xlim(-5,151);c.set_xticks([0,25,50,75,100]);c.set_xlabel('Reported availability (%)');c.spines['left'].set_visible(False);c.tick_params(axis='y',length=0)
    b=fig.add_subplot(gs[:,1]);title(b,'B','Characteristics of 57 signatures')
    char=pd.read_csv(SD/'Fig2_signature_characteristics.csv')
    panels=['Evidence level','Endpoint','Source data type','Omics/platform','Cutoff reporting','Preprocessing details reported','Code/calculator availability','Clinical intended use']
    short={'L3 score portable':'L3','L1 formula reproducible':'L1','L2 externally scoreable':'L2','L5 clinical deployability candidate':'L5 candidate','TCGA+GEO+Private/internal':'TCGA/GEO/private','TCGA+Private/internal':'TCGA/private','GEO+Private/internal':'GEO/private','Private/internal':'Private','Rule reported, numeric value unclear':'Rule; value unclear','Original code/calculator unclear or not reported':'Not reported/unclear','Preprocessing details reported':'Reported','Unclear/not reported':'Not reported/unclear','Numeric cutoff reported':'Numeric','Median/rule only':'Median/rule'}
    colors=[BLUE,TEAL,GOLD,'#B8868D','#8274A1']
    for i,p in enumerate(panels):
        y=(7-i)*2.8; sub=char[char.Panel==p].sort_values('Count',ascending=False);left=0; keys=[]
        b.text(0,y+.7,['Evidence level','Endpoint','Source data','Omics/platform','Cutoff reporting','Preprocessing','Code/calculator','Intended clinical use'][i],fontsize=9,fontweight='bold')
        for j,r in enumerate(sub.itertuples()):
            b.barh(y,r.Count,left=left,height=.52,color=colors[j],edgecolor='white')
            keys.append(f'{j+1}: {short.get(r.Category,r.Category)} (n={r.Count})')
            if r.Count>=5:b.text(left+r.Count/2,y,str(j+1),ha='center',va='center',fontsize=8,color='white')
            left+=r.Count
        b.text(0,y-.6,textwrap.fill('; '.join(keys),65),va='top',fontsize=7.2,linespacing=1.3)
    b.set_yticks([]);b.set_xlim(0,60);b.set_ylim(-2,21.5);b.set_xticks([0,15,30,45,60]);b.set_xlabel('Signature count');b.spines['left'].set_visible(False)
    fig.subplots_adjust(left=.16,right=.97,top=.94,bottom=.08);save(fig,2)

def f3():
    d=pd.read_csv(SD/'Fig3_methodological_applicability_assessment.csv');fig,axes=plt.subplots(1,3,figsize=(11.5,5.8));concerns=['Low concern','Unclear','High concern'];cols=[TEAL,GRAY,RED]
    for ax,group,letter,heading in zip(axes[:2],['Methodological','Applicability/deployability'],['A','B'],['Methodological concerns','Deployability concerns']):
        sub=d[d.Domain_group==group]; domains=list(sub.Domain.unique());title(ax,letter,heading)
        for i,dom in enumerate(domains):
            left=0
            for label,col in zip(concerns,cols):
                v=sub[(sub.Domain==dom)&(sub.Concern==label)].Percent.sum();ax.barh(len(domains)-1-i,v,left=left,color=col,height=.58,edgecolor='white')
                if v>=18:ax.text(left+v/2,len(domains)-1-i,f'{v:.0f}',ha='center',va='center',fontsize=8,color='white' if col!=GRAY else '#222222')
                left+=v
        ax.set_yticks(range(len(domains)),[textwrap.fill(t,21) for t in domains[::-1]],fontsize=8);ax.set_xlim(0,100);ax.set_xlabel('Signatures (%)')
    ax=axes[2];title(ax,'C','Largest evidence gaps')
    t=d[d.Concern.isin(['High concern','Unclear'])].groupby('Domain').Percent.sum().nlargest(5)
    for i,dom in enumerate(t.index):
        high=d[(d.Domain==dom)&(d.Concern=='High concern')].Percent.sum();unclear=d[(d.Domain==dom)&(d.Concern=='Unclear')].Percent.sum()
        ax.barh(4-i,high,color=RED,height=.58);ax.barh(4-i,unclear,left=high,color=GRAY,height=.58)
    ax.set_yticks(range(5),[textwrap.fill(x,20) for x in t.index[::-1]],fontsize=8);ax.set_xlim(0,100);ax.set_xlabel('High concern or unclear (%)')
    fig.legend([Patch(color=x) for x in cols],concerns,loc='lower center',ncol=3)
    fig.subplots_adjust(left=.16,right=.98,wspace=1.15,bottom=.16,top=.9);save(fig,3)

def f4():
    import generate_biorender_style_figures_v1_polished as old
    old.OUT=OUT/'figure_work';old.OUT.mkdir(exist_ok=True);old.DERIVED=old.OUT/'derived';old.DERIVED.mkdir(exist_ok=True)
    def export(fig,stem):
        for t in list(fig.texts):t.remove()
        for ax in fig.axes:
            for t in ax.texts:
                if t.get_text() in ['a','b','c','d']:t.set_text(t.get_text().upper())
            if ax.get_title(loc='left')=='A131 strict cutoff':ax.set_title('A131: fixed cutoff -0.366',loc='left')
        fig.axes[3].set_xticks([.5,1,2,4],['0.5','1','2','4'])
        fig.axes[3].minorticks_off()
        save(fig,4)
    old.save_fig=export;old.plot_figure_4()

def f5():
    c=pd.read_csv(DATA/'figure_5_clinical_increment_data.csv').set_index('Article_ID').reindex(ORDER)
    l=pd.read_csv(DATA/'figure_5_lrt_data.csv').set_index('Article_ID').reindex(ORDER)
    s=pd.read_csv(SUPP/'Scorecard_and_sensitivity/LUAD_deployability_scorecard_primary_candidates.csv').set_index('Article_ID').reindex(ORDER)
    dims=['formula_completeness_score','gene_coverage_score','external_scoreability_score','score_direction_portability_score','cutoff_transportability_score','high_risk_proportion_drift_score','endpoint_concordance_score','clinical_covariate_increment_score','assay_preprocessing_readiness_score','prospective_validation_readiness']
    fig,axs=plt.subplots(2,2,figsize=(9.8,7.8));a,b,cax,d=axs.flat;y=np.arange(10)[::-1]
    title(a,'A','Increment beyond clinical variables');a.hlines(y,c['min'],c['max'],color=GRAY,lw=3);a.scatter(c['median'],y,color=BLUE,zorder=3);a.set_yticks(y,ORDER);a.set_xlabel('Delta C-index (median and cohort range)')
    title(b,'B','Likelihood-ratio tests');b.bar(range(10),l.lrt_p05_fraction,color=BLUE);b.set_xticks(range(10),ORDER,rotation=45,ha='right');b.set_ylim(0,1.12);b.set_ylabel('Fraction with LRT P<0.05')
    title(cax,'C','Deployability scorecard');m=s[dims].to_numpy();cax.imshow(m,cmap=ListedColormap(['#EAB3AF',GRAY,BLUE]),vmin=0,vmax=2,aspect='auto');cax.set_yticks(range(10),ORDER)
    cax.set_xticks(range(10),['Formula','Coverage','Scoring','Direction','Cutoff','Drift','Endpoint','Clinical','Assay','Prospective'],rotation=55,ha='right',fontsize=8)
    for i,total in enumerate(s.total_deployability_score_0_20):cax.text(10,i,f'{int(total)}/20',va='center',fontsize=8)
    cax.set_xlim(-.5,11.5);cax.legend([Patch(color=x) for x in ['#EAB3AF',GRAY,BLUE]],['0','1','2'],title='Domain score',loc='upper center',bbox_to_anchor=(.5,-.37),ncol=3,fontsize=8,title_fontsize=8)
    title(d,'D','Independent evidence domains');d.axis('off');d.set_xlim(0,1);d.set_ylim(0,1)

    # Four compact visual encodings replace prose-heavy definitions. The 2x2
    # arrangement deliberately avoids arrows so the domains are not read as a ladder.
    domains=[
        (.25,.70,BLUE,'Prognostic\nsignal','association'),
        (.75,.70,TEAL,'Score\nportability','ranking'),
        (.25,.27,GOLD,'Cutoff\ntransportability','threshold'),
        (.75,.27,RED,'Clinical\ndeployability','clinical'),
    ]
    for x,y,col,label,glyph in domains:
        d.add_patch(Circle((x,y),.105,facecolor='white',edgecolor=col,lw=2.0))
        if glyph=='association':
            d.plot([x-.061,x+.061],[y,y],color='#6B7280',lw=1.7,solid_capstyle='round')
            d.plot([x-.061,x-.061],[y-.020,y+.020],color='#6B7280',lw=1.2)
            d.plot([x+.061,x+.061],[y-.020,y+.020],color='#6B7280',lw=1.2)
            d.scatter([x+.022],[y],s=34,color=col,zorder=4,edgecolor='white',linewidth=.5)
        elif glyph=='ranking':
            for off,shift in [(.040,.025),(0,.012),(-.040,-.004)]:
                d.plot([x-.055,x+.055],[y+off-shift,y+off+shift],color=col,lw=1.5,alpha=.9)
                d.scatter([x-.055,x+.055],[y+off-shift,y+off+shift],s=9,color=col,zorder=4)
        elif glyph=='threshold':
            xs=np.linspace(x-.072,x+.072,60)
            curve=y-.037+.060*np.exp(-((xs-x)/.041)**2)
            d.plot(xs,curve,color='#6B7280',lw=1.25)
            d.plot([x,x],[y-.055,y+.066],color=col,lw=1.7,ls='--')
            d.scatter([x-.047,x+.050],[y-.044,y-.044],s=16,color=[BLUE,RED],zorder=4)
        else:
            d.add_patch(Rectangle((x-.047,y-.057),.073,.105,facecolor='none',edgecolor=col,lw=1.5))
            d.plot([x-.028,x-.011,x+.047],[y-.004,y-.025,y+.033],color=col,lw=2.0,solid_capstyle='round')
            d.plot([x+.038,x+.067],[y-.048,y-.048],color='#6B7280',lw=1.3)
        d.text(x,y-.145,label,ha='center',va='top',fontsize=8.3,fontweight='bold',color=col,linespacing=1.05)
    fig.subplots_adjust(left=.09,right=.96,top=.92,bottom=.2,hspace=.65,wspace=.45);save(fig,5)

def f6():
    summ=pd.read_csv(DATA/'figure_6_strategy_summary.csv');heat=pd.read_csv(DATA/'figure_6_cohort_heatmap_data.csv');g=pd.read_csv(DATA/'figure_6_gene_stability_data.csv')
    fig,axs=plt.subplots(2,2,figsize=(10,8.3));a,b,c,d=axs.flat;colors=['#7964A8',BLUE,'#448F99',GOLD,TEAL]
    title(a,'A','Strategy-level performance');yy=np.arange(5)[::-1]
    for i,r in enumerate(summ.itertuples()):a.errorbar(r.median_C_index,yy[i],xerr=[[r.median_C_index-r.min_C_index],[r.max_C_index-r.median_C_index]],fmt='o',color=colors[i],capsize=3)
    a.set_yticks(yy,summ.strategy,fontsize=8);a.axvline(.5,color='#999999',ls='--');a.set_xlim(.42,.76);a.set_xlabel('C-index (median and cohort range)')
    title(b,'B','Cohort-level performance');h=heat.pivot(index='strategy',columns='GSE',values='C_index').reindex(summ.strategy);p=heat.pivot(index='strategy',columns='GSE',values='Cox_p').reindex(index=h.index,columns=h.columns)
    im=b.imshow(h,cmap='YlGnBu',vmin=.43,vmax=.74,aspect='auto');b.set_yticks(range(5),h.index,fontsize=7);b.set_xticks(range(9),h.columns,rotation=50,ha='right',fontsize=7)
    for i in range(5):
        for j in range(9):b.text(j,i,f'{h.iloc[i,j]:.2f}'+('*' if p.iloc[i,j]<.05 else ''),ha='center',va='center',fontsize=6.5,color='white' if h.iloc[i,j]>.64 else '#222222')
    fig.colorbar(im,ax=b,fraction=.035,pad=.03,label='C-index')
    title(c,'C','Gene-level stability');x=np.log(g.external_median_HR_per_SD);y=-np.log10(g.tcga_p.clip(lower=1e-300));sc=c.scatter(x,y,c=g.external_p05_n,cmap='viridis',s=30,alpha=.8,edgecolor='white',linewidth=.3)
    selected=g.head(5).index.tolist();offsets=[(-66,30),(-75,10),(12,20),(15,-25),(-72,-32)]
    for idx,off in zip(selected,offsets):c.annotate(g.loc[idx,'gene'],(x[idx],y[idx]),xytext=off,textcoords='offset points',fontsize=8,fontstyle='italic',arrowprops=dict(arrowstyle='-',lw=.5,color='#777777'))
    c.axvline(0,color='#999999',lw=.7);c.set_xlabel('External median log(HR per SD)');c.set_ylabel('-log10(TCGA P)');c.set_xlim(-.45,.5);c.set_ylim(-.3,8.5)
    fig.colorbar(sc,ax=c,fraction=.035,pad=.03,label='External cohorts with P<0.05')
    title(d,'D','Core-score performance and consistency');d.axis('off');t=summ.iloc[3:5]
    rows=[]
    for label,key,fmt in [('Median C-index','median_C_index','c'),('Minimum C-index','min_C_index','c'),('HR >1','HR_gt1_n','n'),('Cox P<0.05','p05_n','n'),('C-index >=0.55','C_ge_055_n','n')]:
        rows.append([label]+[f'{v:.3f}' if fmt=='c' else f'{int(v)}/9' for v in t[key]])
    table=d.table(cellText=rows,colLabels=['Metric','5 genes','40 genes'],colWidths=[.48,.25,.27],cellLoc='center',bbox=[0,.2,1,.76]);table.auto_set_font_size(False);table.set_fontsize(9)
    for (r,k),cell in table.get_celld().items():
        cell.set_edgecolor('#DDDDDD');cell.set_linewidth(.5)
        if r==0:cell.set_facecolor('#EEF2F4');cell.set_text_props(weight='bold')
    d.text(0,.06,'Both core scores are post hoc explorations\nusing external outcome information.',fontsize=8.5,transform=d.transAxes)
    fig.subplots_adjust(left=.19,right=.96,top=.93,bottom=.09,hspace=.64,wspace=.72);save(fig,6)

if __name__=='__main__':
    for func in [f1,f2,f3,f4,f5,f6]:func()
    print('Six revised figures exported from existing data.')
