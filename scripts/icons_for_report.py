import numpy as np
import matplotlib.pyplot as plt
import os
plt.rcParams["figure.figsize"]=(3,3)
plt.rcParams["axes.facecolor"]="white"
plt.rcParams["font.family"]="DejaVu Sans"
blue="#023d6b"
orange="#ff7f00"

outdir="/home/walterh/Desktop/graph_dumb"
os.makedirs(outdir,exist_ok=True)

def save_icon(fig,name): 
    fig.savefig(os.path.join(outdir,name),dpi=600,bbox_inches="tight",transparent=True) 
    plt.close(fig)

def timestamp_handling_icon():
    t=np.array([0,1,2,3,4,9,10,11,12,13],dtype=float)
    y=np.array([1.0,1.2,0.9,1.1,1.05,1.2,1.1,1.15,1.0,1.05])
    fig,ax=plt.subplots()
    ax.plot(t[:5],y[:5],lw=2,color=blue)
    ax.plot(t[5:],y[5:],lw=2,color=blue)
    ax.plot([4.2,8.8],[0.4,0.4],ls=(0,(3,3)),color=orange,lw=2)
    ax.text(6.5,0.5,">2 h gap",ha="center",va="bottom",fontsize=9,color=orange)
    ax.set_xlim(-0.5,13.5); ax.set_ylim(0.3,1.5); ax.axis("off")
    save_icon(fig,"icon_timestamp_handling.png")

def sentinel_values_icon():
    x=np.linspace(0,10,40)
    y=1+0.1*np.sin(x)
    y[8]=-9.999
    fig,ax=plt.subplots()
    ax.plot(x,y,marker="o",ms=3,lw=1,color=blue)
    ax.scatter([x[8]],[y[8]],s=60,edgecolor=orange,facecolor="none",linewidths=2)
    ax.text(x[8],y[8]-0.6,"sentinel",ha="center",va="top",fontsize=9,color=orange)
    y2=y.copy(); y2[y2<-5]=np.nan
    ax.plot(x,y2,lw=3,color=blue,alpha=0.6)
    ax.set_xlim(-0.5,10.5); ax.set_ylim(0.1,1.5); ax.axis("off")
    save_icon(fig,"icon_sentinel_values.png")

def binary_switches_icon():
    x=np.arange(0,12)
    y=(np.sin(x*0.7)>0).astype(float)
    fig,ax=plt.subplots()
    ax.step(x,y,where="post",lw=3,color=blue)
    ax.set_xlim(-0.2,11.2); ax.set_ylim(-0.2,1.2); ax.axis("off")
    ax.text(10.6,1.0,"1",fontsize=9,va="center")
    ax.text(10.6,0.0,"0",fontsize=9,va="center")
    save_icon(fig,"icon_binary_switches.png")

def flat_values_slopes_icon():
    x=np.linspace(0,10,100)
    y=np.piecewise(x,[x<3,(x>=3)&(x<7),x>=7],[lambda x:1.0,lambda x:1.0,lambda x:1.0+0.25*(x-7)])
    fig,ax=plt.subplots()
    ax.plot(x,y,lw=3,color=blue)
    ax.plot([0.2,2.8],[0.8,0.8],lw=3,color=orange)
    ax.text(1.5,0.78,"flat",ha="center",va="top",fontsize=9,color=orange)
    ax.plot([7.0,9.8],[1.0,1.7],lw=2,ls="--",color=orange)
    ax.text(8.4,1.78,"slope",ha="center",va="bottom",fontsize=9,color=orange)
    ax.set_xlim(-0.2,10.2); ax.set_ylim(0.6,1.9); ax.axis("off")
    save_icon(fig,"icon_flat_values_slopes.png")

def decimals_quantization_icon():
    x=np.linspace(0,10,120)
    raw=1+0.3*np.sin(1.1*x)+0.05*np.random.RandomState(3).randn(x.size)
    q=np.round(raw/0.05)*0.05
    fig,ax=plt.subplots()
    ax.plot(x,raw,lw=1,alpha=0.4,color=blue)
    ax.step(x,q,where="mid",lw=3,color=orange)
    ax.text(9.2,1.6,"Δ≈0.05",fontsize=9,ha="right",va="center",color=orange)
    ax.set_xlim(-0.2,10.2); ax.set_ylim(0.5,1.8); ax.axis("off")
    save_icon(fig,"icon_decimals_quantization.png")

timestamp_handling_icon()
sentinel_values_icon()
binary_switches_icon()
flat_values_slopes_icon()
decimals_quantization_icon()
print("Icons saved to:",outdir)












import numpy as np,matplotlib.pyplot as plt,os
plt.rcParams["figure.figsize"]=(3,3);plt.rcParams["axes.facecolor"]="white";plt.rcParams["font.family"]="DejaVu Sans"
blue="#023d6b";orange="#ff7f00";outdir="/home/walterh/Desktop/graph_dumb";os.makedirs(outdir,exist_ok=True)
def save(fig,name): fig.savefig(os.path.join(outdir,name),dpi=600,bbox_inches="tight",transparent=True);plt.close(fig)
def icon_seasonal_coverage():
    months=np.arange(12);cov=0.6+0.35*np.sin((months-1)/12*2*np.pi)+0.05*np.random.RandomState(3).randn(12);cov=np.clip(cov,0,1);nan=1-cov
    fig,ax=plt.subplots();ax.bar(months-0.18,cov,0.36,color=blue);ax.bar(months+0.18,nan,0.36,color=orange);ax.set_xlim(-0.8,11.8);ax.set_ylim(0,1.05);ax.set_xticks([0,2,5,8,11],["J","M","J","S","D"]);ax.set_yticks([]);ax.text(1,-0.12,"DJF",ha="center",va="top");ax.text(4,-0.12,"MAM",ha="center",va="top");ax.text(7,-0.12,"JJA",ha="center",va="top");ax.text(10,-0.12,"SON",ha="center",va="top");ax.legend(["coverage","NaN"],loc="upper right",frameon=False,fontsize=8);ax.axis("off");save(fig,"icon_seasonal_coverage.png")
def icon_plots_per_variable():
    rng=np.random.RandomState(4);x=np.linspace(0,10,120);raw=1+0.3*np.sin(1.1*x)+0.1*rng.randn(x.size);clean=raw.copy();clean[20:28]=np.nan;clean=clean
    fig=plt.figure(figsize=(3,3));gs=fig.add_gridspec(2,2);ax1=fig.add_subplot(gs[0,0]);ax2=fig.add_subplot(gs[0,1]);ax3=fig.add_subplot(gs[1,0]);ax4=fig.add_subplot(gs[1,1])
    ax1.plot(x,raw,lw=1,color=blue,alpha=0.5);ax1.plot(x,clean,lw=2,color=blue);ax1.set_title("raw vs cleaned",fontsize=8);ax1.axis("off")
    vals=np.round(clean[~np.isnan(clean)]/0.01)*0.01;ax2.hist(vals,bins=15,color=blue);ax2.set_title("decimals hist",fontsize=8);ax2.axis("off")
    d=np.abs(np.diff(np.nan_to_num(clean,nan=np.nanmean(clean))));ax3.plot(x[1:],d,lw=2,color=blue);ax3.set_title("|Δ|",fontsize=8);ax3.axis("off")
    H=np.outer(np.linspace(0.2,1,6),np.linspace(0.2,1,10));ax4.imshow(H,aspect="auto",origin="lower");ax4.set_title("seasonal heatmap",fontsize=8);ax4.axis("off");plt.tight_layout(pad=0.6);save(fig,"icon_plots_per_variable.png")
def icon_sentemqc():
    rng=np.random.RandomState(1);t=np.linspace(0,10,80);y=1+0.1*t+0.12*rng.randn(t.size);mean=1+0.1*t;band=0.25
    fig,ax=plt.subplots();ax.fill_between(t,mean-band,mean+band,color=blue,alpha=0.12);ax.plot(t,mean,lw=2,color=blue,alpha=0.7);ax.scatter(t,y,s=12,color=blue);mask=(y>mean+band)|(y<mean-band);ax.scatter(t[mask],y[mask],s=30,facecolors="none",edgecolors=orange,linewidths=2);ax.text(9.8,mean[-1]+band,"QC band",ha="right",va="bottom",fontsize=9);ax.axis("off");save(fig,"icon_sentemqc.png")
def icon_wrtds_proxy():
    rng=np.random.RandomState(2)
    t=np.linspace(0,10,90)
    true=1+0.4*np.sin(0.6*t)
    y=true+0.15*rng.randn(t.size)
    fit=1+0.38*np.sin(0.6*t+0.1)
    res=y-fit
    thr=1.5*np.median(np.abs(res-np.median(res)))
    flag=np.abs(res)>thr
    fig,ax=plt.subplots()
    ax.plot(t,fit,lw=2,color=blue)
    ax.scatter(t,y,s=10,color=blue,alpha=0.7)
    for ti,yi,fi,fl in zip(t,y,fit,flag):
        if fl:
            ax.vlines(ti,fi,yi,colors=orange,lw=2)
    ax.text(9.9,fit[-1],"local regression",ha="right",va="center",fontsize=9)
    ax.text(0.1,min(y.min(),fit.min())-0.05,"|residual|>robust z",ha="left",va="top",fontsize=8,color=orange)
    ax.axis("off")
    save(fig,"icon_wrtds_proxy.png")





icon_seasonal_coverage();icon_plots_per_variable();icon_sentemqc();icon_wrtds_proxy();print("Saved to",outdir)













