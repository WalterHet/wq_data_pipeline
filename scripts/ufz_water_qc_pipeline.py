# ufz_water_qc_pipeline.py
# (all SaQC calls audited and made version-safe; fallbacks provided if a method is missing.
#  Windows now use lowercase 'h' instead of deprecated 'H' everywhere.)

import os,math,re,warnings,stat
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator
try:
    import scipy.signal as spsig
    import scipy.ndimage as spimg
    from scipy.stats import chisquare
    SCIPY_OK=True
except Exception:
    SCIPY_OK=False
try:
    import saqc as _saqc
    SAQC_OK=True
except Exception:
    SAQC_OK=False
try:
    import sentemqc as sm
    SENTEM_OK=True
except Exception:
    SENTEM_OK=False

plt.rcParams['figure.dpi']=140
plt.rcParams['font.size']=10
plt.rcParams['axes.titlesize']=12
plt.rcParams['axes.labelsize']=11
plt.rcParams['legend.fontsize']=9
plt.rcParams['xtick.labelsize']=9
plt.rcParams['ytick.labelsize']=9

COLOR_RAW='#9e9e9e'
COLOR_BLUE='#023d6b'
COLOR_ORANGE='#ff7f00'
COLOR_GREY='#777777'
COLOR_GREEN='#0a7f2f'
COLOR_RED='#b22222'
COLOR_PURPLE='#6a51a3'
COLOR_BLACK='#000000'
CMAP_HEAT=ListedColormap(['#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#084594'])

pd.options.display.width=180
pd.options.display.max_columns=200
pd.options.display.max_rows=200

""" SETUP PATHS AND STRUCTURE """
def build_structure(base_root,station):
    root=Path(base_root)
    figs=root/'reports'/'figures'/sanitize(station)
    tabs=root/'reports'/'tables'/sanitize(station)
    logs=root/'reports'/'logs'/sanitize(station)
    proot=root/'processed'/sanitize(station)
    for p in [figs,tabs,logs,proot]:p.mkdir(parents=True,exist_ok=True)
    return {'figs':figs,'tabs':tabs,'logs':logs,'processed':proot,'root':root}

""" BASIC UTILS """
def sanitize(s):
    return re.sub(r'[^A-Za-z0-9._-]+','_',str(s)).strip('_')

def ensure_time_index(df,ts='timestamp',tz=None):
    if isinstance(df.index,pd.DatetimeIndex):idx=df.index
    else:idx=pd.to_datetime(df[ts],errors='coerce')
    if tz is not None and idx.tz is None:idx=idx.tz_localize(tz,ambiguous='NaT',nonexistent='NaT')
    df=df.copy()
    df.index=idx
    df=df[~df.index.duplicated(keep='first')].sort_index()
    return df

def infer_step(df_or_index):
    idx=df_or_index.index.view('int64') if hasattr(df_or_index,'index') else pd.Index(df_or_index).view('int64')
    if idx.size<3:return pd.Timedelta('15min')
    d=np.diff(idx);d=d[d>0]
    if d.size==0:return pd.Timedelta('15min')
    return pd.to_timedelta(int(np.median(d)),'ns')

def as_season(ts):
    m=ts.month
    if m in (12,1,2):return 'DJF'
    if m in (3,4,5):return 'MAM'
    if m in (6,7,8):return 'JJA'
    return 'SON'

def season_order():
    return ['DJF','MAM','JJA','SON']

def rle_bool(x):
    n=len(x)
    if n==0:return np.array([],dtype=int),np.array([],dtype=int),np.array([],dtype=bool)
    xb=np.asarray(x,dtype=bool)
    dif=np.diff(np.concatenate(([True],xb[1:]!=xb[:-1],[True])))
    idx=np.flatnonzero(dif)
    run_starts=idx[:-1]
    run_lengths=np.diff(idx)
    run_values=xb[run_starts]
    return run_starts,run_lengths,run_values

def rolling_lin_slope(y,x=None,win=25,minp=5):
    v=np.asarray(y,dtype=float)
    if x is None:
        t=np.arange(v.size,dtype=float)
    else:
        t=(pd.to_datetime(x).view('int64').astype(float)-pd.to_datetime(x).view('int64').astype(float).min())/1e9
    out=np.full(v.size,np.nan)
    k=int(win)
    if k<3:k=3
    h=k//2
    for i in range(v.size):
        a=max(0,i-h);b=min(v.size,i+h+1)
        seg=v[a:b];ts=t[a:b]
        m=np.isfinite(seg)
        if m.sum()>=minp:
            X=np.vstack([ts[m],np.ones(m.sum())]).T
            beta=np.linalg.lstsq(X,seg[m],rcond=None)[0]
            out[i]=beta[0]
    return out

def dec_frac(x):
    v=np.asarray(x,dtype=float)
    return np.modf(v)[0]%1.0

def quant_step_estimate(x,max_bins=50):
    v=pd.Series(np.asarray(x,dtype=float))
    d=v.diff().dropna().abs()
    d=d[(d>0)&np.isfinite(d)]
    if d.empty:return np.nan
    q=np.quantile(d,[0.1,0.25,0.5,0.75,0.9])
    cand=[q[0],q[1],q[2]]+[q[2]/i for i in range(2,10)]
    cand=[c for c in cand if c>0]
    best=np.nan;besth=0
    for c in cand:
        if c<=0:continue
        h=np.mean(np.isclose((d/c)-np.round(d/c),0,atol=1e-3))
        if h>besth:besth=h;best=c
    return best if besth>0.7 else np.nan

def _to_points(window_str, step):
    """convert '2h' style window to #points given sampling step (Timedelta)"""
    w=pd.to_timedelta(window_str)
    return max(1,int(round(w / step)))

""" DATA INGESTION """
def read_data(csv_path,timestamp_col='timestamp',station_field=None,station_value=None,usecols=None,parse_dates=True):
    df=pd.read_csv(csv_path,parse_dates=[timestamp_col] if parse_dates else None,usecols=usecols)
    df=ensure_time_index(df,ts=timestamp_col,tz=None)
    if station_field is not None and station_field in df.columns and station_value is not None:
        df=df[df[station_field]==station_value]
    return df

""" SENTINEL HANDLING """
def detect_and_mask_sentinels(s,values=(-9,-99,-999,-9999),min_count=5,min_frac=0.001):
    x=pd.Series(s).astype(float)
    counts={v:int(np.sum(x==v)) for v in values}
    n=len(x)
    use=set([v for v,c in counts.items() if (c>=min_count) or (n>0 and c/n>=min_frac)])
    m=np.zeros(len(x),dtype=bool)
    if len(use)>0:
        for v in use:m|=(x==v)
        x=x.mask(m,np.nan)
    return x,m,use

""" DUPLICATES """
def resolve_duplicates(df,how='median'):
    idx=df.index
    dup=idx.duplicated(keep=False)
    if not dup.any():return df,0
    aggf={'median':np.nanmedian,'mean':np.nanmean}.get(how,np.nanmedian)
    grouped=df.groupby(level=0).agg(aggf)
    return grouped,dup.sum()

""" GAP CLASSIFICATION """
def classify_gaps(idx,hours=2.0):
    t=pd.to_datetime(idx).view('int64')
    if t.size<2:return pd.Series(False,index=idx),pd.Series(pd.NaT,index=idx)
    d=np.diff(t)
    thr=int(hours*3600*1e9)
    g=np.concatenate(([False],d>thr))
    gap_time=pd.Series(pd.NaT,index=idx)
    if g.any():
        where=np.flatnonzero(g)
        for w in where:gap_time.iloc[w]=pd.to_datetime(t[w]-t[w-1])
    return pd.Series(g,index=idx),gap_time

""" BINARY SWITCHES """
def detect_binary_switches(x,zero_tol=1e-12,min_dur=1):
    v=np.asarray(x,dtype=float)
    z=np.isfinite(v)&(np.abs(v)<=zero_tol)
    s,l,val=rle_bool(z)
    events=[]
    for i in range(len(s)):
        if val[i] and l[i]>=min_dur:
            a=s[i];b=s[i]+l[i]-1
            pre=max(0,a-1);post=min(len(v)-1,b+1)
            preval=v[pre] if np.isfinite(v[pre]) else np.nan
            postval=v[post] if np.isfinite(v[post]) else np.nan
            if np.isfinite(preval) and np.isfinite(postval) and (preval!=0 or postval!=0):
                events.append({'start_idx':a,'end_idx':b,'pre':pre,'post':post,'pre_val':preval,'post_val':postval})
    return events

""" FLAT VALUES AND FLAT SLOPES """
def detect_flat_runs(x,dt_index,min_hours=2.0,abs_tol=0.0):
    v=pd.Series(x,index=dt_index).astype(float)
    eq=np.isclose(v.values,np.roll(v.values,1),atol=abs_tol,rtol=0.0)
    eq[0]=False
    s,l,val=rle_bool(eq)
    runs=[]
    for i in range(len(s)):
        if val[i]:
            a=s[i];b=s[i]+l[i]
            t1=v.index[a];t2=v.index[b] if b<len(v) else v.index[-1]
            dur=(t2-t1).total_seconds()/3600.0
            if dur>=min_hours:
                runs.append({'start':t1,'end':t2,'hours':dur,'value':float(v.iloc[a])})
    return runs

def detect_flat_slopes(x,dt_index,min_hours=2.0,win=25,abs_slope=0.0):
    s=rolling_lin_slope(x,x=dt_index,win=win,minp=max(5,win//3))
    v=pd.Series(s,index=dt_index)
    eq=np.isfinite(v.values)&(np.abs(v.values)<=abs_slope)
    s2,l2,val2=rle_bool(eq)
    runs=[]
    for i in range(len(s2)):
        if val2[i]:
            a=s2[i];b=s2[i]+l2[i]
            t1=v.index[a];t2=v.index[b] if b<len(v) else v.index[-1]
            dur=(t2-t1).total_seconds()/3600.0
            if dur>=min_hours:
                runs.append({'start':t1,'end':t2,'hours':dur,'slope':0.0})
    return runs

""" DECIMAL AND QUANTIZATION TESTS """
def decimal_uniformity(x,nbins=10):
    f=dec_frac(x)
    if not np.isfinite(f).any():return {'chisq':np.nan,'p':np.nan,'hist':np.zeros(nbins),'edges':np.linspace(0,1,nbins+1)}
    h,edges=np.histogram(f,bins=np.linspace(0,1,nbins+1))
    if SCIPY_OK:
        cs,p=chisquare(h)
    else:
        e=np.full_like(h,h.mean(),dtype=float)
        cs=np.sum((h-e)**2/(e+1e-9));p=np.nan
    return {'chisq':float(cs),'p':float(p),'hist':h,'edges':edges}

def find_long_unique_double_decimals(x,threshold_frac=0.6):
    f=np.round(dec_frac(x)*100).astype(int)
    f=f[np.isfinite(f)]
    if f.size==0:return {'dom':None,'frac':0.0}
    vals,cts=np.unique(f,return_counts=True)
    j=np.argmax(cts)
    dom=int(vals[j]);frac=float(cts[j]/f.size)
    if frac>=threshold_frac:return {'dom':dom,'frac':frac}
    return {'dom':dom,'frac':frac}

def infer_quantization(x):
    st=quant_step_estimate(x)
    return {'step':st,'is_quantized':bool(np.isfinite(st))}

""" SEASONAL COVERAGE AND SUMMARIES """
def expected_points_per_season(index,step):
    df=pd.DataFrame({'season':[as_season(t) for t in index]},index=index)
    res={}
    for s in season_order():
        m=(df['season']==s)
        if m.any():
            dur=(df.index[m][-1]-df.index[m][0]).total_seconds()+step.total_seconds()
            res[s]=max(1,int(round(dur/step.total_seconds())))
        else:
            res[s]=0
    return res

def seasonal_stats(x,index,step,flat_runs,flat_slope_runs):
    df=pd.DataFrame({'val':x},index=index)
    df['is_nan']=~np.isfinite(df['val'])
    df['season']=[as_season(t) for t in index]
    out=[]
    exp=expected_points_per_season(index,step)
    for s in season_order():
        sub=df[df['season']==s]
        n=len(sub);nnan=int(sub['is_nan'].sum())
        cov=0.0 if exp[s]==0 else 100.0*(n-nnan)/max(1,exp[s])
        fval=sum(r['hours']>=2.0 for r in flat_runs if as_season(r['start'])==s)
        fslp=sum(r['hours']>=2.0 for r in flat_slope_runs if as_season(r['start'])==s)
        out.append({'season':s,'coverage_pct':cov,'nan_pct':0.0 if n==0 else 100.0*nnan/max(1,n),'flat_values_events':fval,'flat_slopes_events':fslp})
    return pd.DataFrame(out).set_index('season').loc[season_order()]

""" PLOTTING BASICS """
def set_xaxis_dates(ax):
    loc=mdates.AutoDateLocator()
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))

def plot_series_basic(idx,y,title,ylabel,savepath=None,log=False,lines=None,markers=None,legend=True,annot=None):
    fig=plt.figure(figsize=(12,4))
    ax=plt.gca()
    if markers is None:markers=[]
    if lines is None:lines=[]
    ax.plot(idx,y,color=COLOR_RAW,lw=0.7,label='Raw')
    for ln in lines:
        ax.plot(idx,ln.get('y'),lw=ln.get('lw',0.8),color=ln.get('color',COLOR_BLUE),label=ln.get('label',''))
    for mk in markers:
        ax.plot(idx,mk.get('x'),mk.get('y'),mk.get('style','.'),ms=mk.get('ms',1.5),color=mk.get('color',COLOR_ORANGE),label=mk.get('label',''))
    if log:
        ax.set_yscale('log')
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Time')
    set_xaxis_dates(ax)
    if legend:ax.legend(ncol=5,frameon=False,loc='best')
    if annot is not None:ax.text(0.01,0.02,annot,transform=ax.transAxes,ha='left',va='bottom')
    plt.tight_layout()
    if savepath is not None:plt.savefig(savepath,bbox_inches='tight')
    plt.close(fig)

def plot_decimal_hist(frac_info,savepath=None):
    h=frac_info['hist'];e=frac_info['edges']
    fig=plt.figure(figsize=(7,3))
    ax=plt.gca()
    centers=0.5*(e[:-1]+e[1:])
    ax.bar(centers,h,width=(e[1]-e[0])*0.95,color=COLOR_BLUE,edgecolor='none')
    ax.set_xlabel('Fractional part [0..1)')
    ax.set_ylabel('Count')
    ax.set_title('Decimal distribution')
    plt.tight_layout()
    if savepath is not None:plt.savefig(savepath,bbox_inches='tight')
    plt.close(fig)

def plot_quant_diffs(x,savepath=None):
    s=pd.Series(x).astype(float).diff().abs()
    fig=plt.figure(figsize=(7,3))
    ax=plt.gca()
    ax.plot(s.index,s.values,lw=0.6,color=COLOR_BLUE)
    set_xaxis_dates(ax)
    ax.set_title('Absolute successive differences')
    ax.set_xlabel('Time');ax.set_ylabel('|Δ|')
    plt.tight_layout()
    if savepath is not None:plt.savefig(savepath,bbox_inches='tight')
    plt.close(fig)

def plot_flat_events(idx,x,flat_runs,flat_slope_runs,savepath=None):
    fig=plt.figure(figsize=(12,4))
    ax=plt.gca()
    ax.plot(idx,x,color=COLOR_RAW,lw=0.7,label='Raw')
    for r in flat_runs:ax.axvspan(r['start'],r['end'],color=COLOR_ORANGE,alpha=0.2)
    for r in flat_slope_runs:ax.axvspan(r['start'],r['end'],color=COLOR_BLUE,alpha=0.18)
    ax.set_title('Flat values (orange) and flat slopes (blue)')
    set_xaxis_dates(ax)
    ax.set_xlabel('Time');ax.set_ylabel('Value')
    plt.tight_layout()
    if savepath is not None:plt.savefig(savepath,bbox_inches='tight')
    plt.close(fig)

def plot_heatmap_coverage(df_cov,years,savepath=None,title='Coverage per season [%]'):
    arr=np.array([df_cov.loc[s,'coverage_pct'] if s in df_cov.index else 0.0 for s in season_order()])[:,None]
    fig=plt.figure(figsize=(4,3))
    ax=plt.gca()
    im=ax.imshow(arr,aspect='auto',cmap=CMAP_HEAT,vmin=0,vmax=100)
    ax.set_yticks(range(len(season_order())));ax.set_yticklabels(season_order())
    ax.set_xticks([0]);ax.set_xticklabels([years])
    ax.set_title(title)
    cb=plt.colorbar(im,ax=ax);cb.ax.set_ylabel('%',rotation=90)
    plt.tight_layout()
    if savepath is not None:plt.savefig(savepath,bbox_inches='tight')
    plt.close(fig)

def plot_gaussian_compare(idx,x,sigma_hours,step,savepath=None):
    if SCIPY_OK:
        y=pd.Series(x,index=idx).astype(float)
        n=len(y)
        if n==0:return
        dt=max(1,int(round(pd.Timedelta(step).total_seconds()/3600.0)))
        sig=max(1,int(round(sigma_hours/dt)))
        ys=spimg.gaussian_filter1d(y.fillna(method='pad').fillna(method='bfill').values,sigma=sig,mode='nearest')
        res=y.values-ys
        fig=plt.figure(figsize=(12,6))
        ax=plt.subplot(2,1,1)
        ax.plot(idx,y.values,color=COLOR_RAW,lw=0.6,label='Raw')
        ax.plot(idx,ys,color=COLOR_BLUE,lw=0.9,label='Gaussian')
        ax.set_title(f'Gaussian filter σ≈{sigma_hours}h');set_xaxis_dates(ax);ax.legend(frameon=False,loc='best')
        ax=plt.subplot(2,1,2)
        ax.plot(idx,res,color=COLOR_ORANGE,lw=0.7,label='Residual');set_xaxis_dates(ax);ax.legend(frameon=False,loc='best')
        plt.tight_layout()
        if savepath is not None:plt.savefig(savepath,bbox_inches='tight')
        plt.close(fig)

""" SAQC ADAPTER (VERSION-SAFE) """
def _mask_from_runs(runs,index):
    if len(index)==0:return pd.Series([],dtype=bool,index=index)
    m=np.zeros(len(index),dtype=bool)
    for r in runs:
        m |= (index>=r['start']) & (index<=r['end'])
    return pd.Series(m,index=index)

def _fallback_flag_plateau(series, min_length='2h'):
    step=infer_step(series)
    hours=pd.to_timedelta(min_length).total_seconds()/3600.0
    runs=detect_flat_runs(series.values, series.index, min_hours=hours, abs_tol=0.0)
    return _mask_from_runs(runs, series.index)

def _fallback_flag_constants(series, window='2h'):
    # treat as short plateaus >= window
    return _fallback_flag_plateau(series, min_length=window)

def _fallback_flag_by_variance(series, window='6h', thresh=1e-12):
    step=infer_step(series)
    k=_to_points(window, step)
    v=pd.Series(series,copy=True).rolling(k,min_periods=max(3,k//3),center=True).var()
    return (v<=thresh).fillna(False)

def _fallback_flag_zscore(series, window='24h', thresh=4.0):
    step=infer_step(series)
    k=_to_points(window, step)
    s=pd.Series(series,copy=True)
    med=s.rolling(k,min_periods=max(5,k//3),center=True).median()
    mad=(np.abs(s-med)).rolling(k,min_periods=max(5,k//3),center=True).median()
    z=0.6745*(s-med)/(mad.replace(0,np.nan))
    return (np.abs(z)>=thresh).fillna(False)

def _fallback_flag_jumps(series, thresh=5.0, window='1h'):
    # simple first-difference threshold
    s=pd.Series(series,copy=True)
    d=s.diff().abs()
    return (d>=thresh).fillna(False)

def _fallback_flag_isolated(series, gap_window='3h', group_window='30min'):
    idx=series.index
    step=infer_step(series)
    gw=pd.to_timedelta(gap_window)
    gw_pts=_to_points(gap_window, step)
    grp_pts=_to_points(group_window, step)
    # find gaps > gw, then runs of finite between big gaps with short length
    finite=np.isfinite(series.values)
    starts,lengths,vals=rle_bool(finite)
    m=np.zeros(len(idx),dtype=bool)
    # compute time diffs to detect big gaps
    t=idx.view('int64')
    dif=np.diff(t)
    big_gap=np.concatenate(([False],dif>gw.value))
    for i in range(len(starts)):
        if not vals[i]:continue
        a=starts[i];b=a+lengths[i]-1
        # surrounding gaps:
        left_big = (a>0) and big_gap[a]
        right_big = (b+1<len(big_gap)) and big_gap[b+1]
        short = lengths[i] <= grp_pts
        if short and left_big and right_big:
            m[a:b+1]=True
    return pd.Series(m,index=idx)

def _fallback_flag_unilof(series, n=20, thresh=2.0, probability=False):
    # fallback to robust z-score
    return _fallback_flag_zscore(series, window='24h', thresh=4.0)

def apply_saqc_suite(df,col,freq=None,range_min=None,range_max=None,plot_path=None,scheme='simple'):
    data=df[[col]].copy()
    if freq is not None:
        data=data.resample(freq).median()
    idx=data.index
    # collect fallback masks (OR’ed later)
    fallback_masks=[]
    saqc_mask=np.zeros(len(idx),dtype=bool)
    series_out=data[col].copy()

    if SAQC_OK:
        qc=_saqc.SaQC(data=data,scheme=scheme)
        try:
            if (range_min is not None) or (range_max is not None):
                qc=qc.flagRange(col,min=range_min,max=range_max)
        except Exception:
            s=series_out
            fb=((s<range_min) | (s>range_max)) if (range_min is not None and range_max is not None) else ((s<range_min) if range_min is not None else (s>range_max))
            fallback_masks.append(fb.fillna(False))

        try:
            qc=qc.flagMissing(col)
        except Exception:
            fallback_masks.append(series_out.isna())

        # flagConstants
        try:
            if hasattr(qc,'flagConstants'):
                qc=qc.flagConstants(col,thresh=0.0,window='2h')
            else:
                raise AttributeError
        except Exception:
            fallback_masks.append(_fallback_flag_constants(series_out,'2h'))

        # flagPlateau
        try:
            if hasattr(qc,'flagPlateau'):
                qc=qc.flagPlateau(col,min_length='2h')
            else:
                raise AttributeError
        except Exception:
            fallback_masks.append(_fallback_flag_plateau(series_out,'2h'))

        # flagByVariance
        try:
            if hasattr(qc,'flagByVariance'):
                qc=qc.flagByVariance(col,window='6h',thresh=1e-12)
            else:
                raise AttributeError
        except Exception:
            fallback_masks.append(_fallback_flag_by_variance(series_out,'6h',1e-12))

        # flagZScore
        try:
            if hasattr(qc,'flagZScore'):
                qc=qc.flagZScore(col,method='rolling',window='24h',thresh=4.0)
            else:
                raise AttributeError
        except Exception:
            fallback_masks.append(_fallback_flag_zscore(series_out,'24h',4.0))

        # flagJumps
        try:
            if hasattr(qc,'flagJumps'):
                qc=qc.flagJumps(col,thresh=5.0,window='1h')
            else:
                raise AttributeError
        except Exception:
            fallback_masks.append(_fallback_flag_jumps(series_out,thresh=5.0,window='1h'))

        # flagIsolated
        try:
            if hasattr(qc,'flagIsolated'):
                qc=qc.flagIsolated(col,gap_window='3h',group_window='30min')
            else:
                raise AttributeError
        except Exception:
            fallback_masks.append(_fallback_flag_isolated(series_out,'3h','30min'))

        # flagUniLOF
        try:
            if hasattr(qc,'flagUniLOF'):
                qc=qc.flagUniLOF(col,n=20,thresh=2.0,probability=False)
            else:
                raise AttributeError
        except Exception:
            fallback_masks.append(_fallback_flag_unilof(series_out,n=20,thresh=2.0,probability=False))

        # Extract SaQC flags (object or numeric), then to boolean mask
        try:
            pf=qc.flags[col]
            if getattr(pf,'dtype',None)==object:
                saqc_mask = (pf.astype(str)!='UNFLAGGED').values
            else:
                saqc_mask = (pf.fillna(0).values!=0)
        except Exception:
            saqc_mask=np.zeros(len(idx),dtype=bool)

        # Plot if requested (SaQC’s own plot)
        if plot_path is not None:
            try: qc.plot(col,path=str(plot_path))
            except Exception: pass

        series_out=qc.data[col].copy()

    else:
        # No SaQC available: only fallbacks
        if (range_min is not None) or (range_max is not None):
            s=series_out
            fb=((s<range_min) | (s>range_max)) if (range_min is not None and range_max is not None) else ((s<range_min) if range_min is not None else (s>range_max))
            fallback_masks.append(fb.fillna(False))
        fallback_masks.append(series_out.isna())
        fallback_masks.append(_fallback_flag_constants(series_out,'2h'))
        fallback_masks.append(_fallback_flag_plateau(series_out,'2h'))
        fallback_masks.append(_fallback_flag_by_variance(series_out,'6h',1e-12))
        fallback_masks.append(_fallback_flag_zscore(series_out,'24h',4.0))
        fallback_masks.append(_fallback_flag_jumps(series_out,thresh=5.0,window='1h'))
        fallback_masks.append(_fallback_flag_isolated(series_out,'3h','30min'))
        fallback_masks.append(_fallback_flag_unilof(series_out,n=20,thresh=2.0,probability=False))

    # Combine SaQC mask with fallbacks
    if len(fallback_masks)>0:
        fb_stack=np.column_stack([m.values if hasattr(m,'values') else np.asarray(m) for m in fallback_masks])
        fb_mask=fb_stack.any(axis=1)
    else:
        fb_mask=np.zeros(len(idx),dtype=bool)

    combined=saqc_mask | fb_mask
    flags_out=pd.Series(np.where(combined,255,0),index=idx)

    return series_out,flags_out

""" WRTDS-PROXY AND BUSTER """
def wrtds_proxy(idx,y,discharge=None,ht_days=90,hq_frac=0.1,hs_season=0.5,minp=30):
    t=pd.to_datetime(idx)
    tnum=(t.view('int64')-t.view('int64').min())/86400e9
    s=np.sin(2*np.pi*t.dayofyear/365.25);c=np.cos(2*np.pi*t.dayofyear/365.25)
    if discharge is None:
        X=np.column_stack([tnum,s,c])
    else:
        q=pd.Series(discharge,index=t).astype(float).values
        X=np.column_stack([tnum,s,c,q])
    y=np.asarray(y,dtype=float)
    m=np.isfinite(y)&np.all(np.isfinite(X),axis=1)
    yv=y[m];Xv=X[m];tv=tnum[m]
    if yv.size<max(minp,10):return np.full_like(y,np.nan),np.full_like(y,np.nan)
    yhat=np.full_like(y,np.nan);res=np.full_like(y,np.nan)
    for i in range(y.size):
        if not m[i]:continue
        dt=np.abs(tv-tnum[i]);wt=np.exp(-(dt/ht_days)**2)
        if X.shape[1]==4:
            dq=np.abs(Xv[:,3]-X[i,3]);wq=np.exp(-(dq/max(1e-9,np.nanmedian(np.abs(Xv[:,3]-np.nanmedian(Xv[:,3])))*3))**2)
            w=wt*wq
        else:
            w=wt
        ws=w/np.nanmax(w) if np.nanmax(w)>0 else w
        idxw=ws>1e-3
        if idxw.sum()<minp:continue
        A=np.column_stack([np.ones(idxw.sum()),Xv[idxw]])
        W=np.diag(ws[idxw])
        beta=np.linalg.lstsq(W@A,W@yv[idxw],rcond=None)[0]
        yhat[i]=(np.array([1.0]+list(X[i]))@beta)
        res[i]=y[i]-yhat[i]
    return yhat,res

def wrtds_buster(idx,y,residuals,z_thresh=4.0,plot_path=None):
    r=pd.Series(residuals,index=idx)
    z=(r-r.median())/(1.4826*(np.abs(r-r.median())).median()+1e-9)
    spikes=np.abs(z)>=z_thresh
    fig=plt.figure(figsize=(12,5))
    ax=plt.gca()
    ax.plot(idx,y,color=COLOR_RAW,lw=0.6,label='Raw')
    ax.plot(idx,r,color=COLOR_BLUE,lw=0.8,label='Residual')
    ax.plot(idx,np.where(spikes,r,np.nan),'.',ms=2.0,color=COLOR_ORANGE,label='Anomaly')
    set_xaxis_dates(ax)
    ax.set_title('WRTDS-proxy residuals and anomalies')
    ax.legend(frameon=False,loc='best')
    plt.tight_layout()
    if plot_path is not None:plt.savefig(plot_path,bbox_inches='tight')
    plt.close(fig)
    return spikes

""" EVENT TABLES AND MERGE """
def events_from_binary(binary_events,idx):
    rec=[]
    for e in binary_events:
        try:
            rec.append({'type':'binary_switch','start':idx[int(e['start_idx'])],'end':idx[int(e['end_idx'])],'pre':idx[int(e['pre'])],'post':idx[int(e['post'])],'pre_val':e['pre_val'],'post_val':e['post_val']})
        except Exception:
            continue
    return pd.DataFrame(rec) if len(rec)>0 else pd.DataFrame(columns=['type','start','end','pre','post','pre_val','post_val'])

def events_from_runs(runs,typ):
    rec=[]
    for r in runs:rec.append({'type':typ,'start':r['start'],'end':r['end'],'hours':r['hours']})
    return pd.DataFrame(rec) if len(rec)>0 else pd.DataFrame(columns=['type','start','end','hours'])

def merge_flags(*arrays):
    arr=[np.asarray(a) for a in arrays if a is not None]
    if len(arr)==0:return None
    m=None
    for a in arr:
        if m is None:m=np.zeros_like(a,dtype=bool)
        m|=(a.astype(bool))
    return m

""" PIPELINE PER VARIABLE """
def process_variable(df,col,station,outs,range_min=None,range_max=None,zero_tol=1e-12,gap_hours=2.0,flat_hours=2.0,flat_slope_win=25,flat_slope_abs=0.0,dec_bins=10,gauss_sigma_hours=6.0,apply_saqc=True,apply_sentem=True,wrtds_q_col=None):
    s0=df[col].astype(float)
    s1,m_sent,used=detect_and_mask_sentinels(s0)
    df_tmp=pd.DataFrame({col:s1},index=df.index)
    df_tmp,dup_count=resolve_duplicates(df_tmp,how='median')
    idx=df_tmp.index
    step=infer_step(df_tmp)
    gaps,gap_dt=classify_gaps(idx,hours=gap_hours)
    x=df_tmp[col].copy()
    x[gaps]=np.nan

    bin_events=detect_binary_switches(x.values,zero_tol=zero_tol,min_dur=1)
    flat_runs=detect_flat_runs(x.values,idx,min_hours=flat_hours,abs_tol=0.0)
    flat_slope_runs=detect_flat_slopes(x.values,idx,min_hours=flat_hours,win=flat_slope_win,abs_slope=flat_slope_abs)

    dec_info=decimal_uniformity(x.values,nbins=dec_bins)
    quant_info=infer_quantization(x.values)
    dom_dec=find_long_unique_double_decimals(x.values,threshold_frac=0.6)
    seas=seasonal_stats(x.values,idx,step,flat_runs,flat_slope_runs)

    p_series=None;p_flags=None
    if apply_saqc:
        p_series,p_flags=apply_saqc_suite(pd.DataFrame({col:x},index=idx),col,freq=None,range_min=range_min,range_max=range_max,plot_path=outs['figs']/f"{sanitize(col)}__saqc_plot.png",scheme='simple')

    sm_df=None
    if apply_sentem and SENTEM_OK:
        try:
            smap={e['col']:(e['code'],e['is_nitrate'],e['label'],e['unit']) for e in sm.VARIABLE_MAP}
            if col in smap:
                code,is_nitrate,lab,unit=smap[col]
                cfg=sm.SENTEMQC_CONFIG.get(code)
                df_sm=sm.apply_sentemqc_to_series(df[col].astype(float),code,cfg,is_nitrate)
                sm_df=df_sm
        except Exception:
            sm_df=None

    figs=outs['figs']/sanitize(col)
    tabs=outs['tabs']
    figs.mkdir(parents=True,exist_ok=True)

    plot_series_basic(idx,s0.values,f'{col} Raw',col,figs/f'{sanitize(col)}__00_raw.png',log=False,legend=False)
    ann=f'Sentinels used: {sorted(list(used))}  dup={dup_count}  step≈{str(step)}'
    plot_series_basic(idx,x.values,f'{col} After sentinels+gaps',col,figs/f'{sanitize(col)}__01_clean.png',log=False,legend=True,annot=ann)
    plot_decimal_hist(dec_info,figs/f'{sanitize(col)}__02_decimal_hist.png')
    plot_quant_diffs(pd.Series(x,index=idx),figs/f'{sanitize(col)}__03_quant_diffs.png')
    plot_flat_events(idx,x.values,flat_runs,flat_slope_runs,figs/f'{sanitize(col)}__04_flat_events.png')
    plot_gaussian_compare(idx,x.values,gauss_sigma_hours,step,figs/f'{sanitize(col)}__05_gaussian_compare.png')
    cov_years=f"{idx.min().year}-{idx.max().year}"
    plot_heatmap_coverage(seas,cov_years,figs/f'{sanitize(col)}__06_seasonal_coverage.png',title=f'{station} — {col} coverage [%]')

    wr_spikes=None
    if wrtds_q_col is not None and wrtds_q_col in df.columns:
        yhat,res=wrtds_proxy(idx,x.values,discharge=df[wrtds_q_col].values,ht_days=90,hq_frac=0.1,hs_season=0.5,minp=30)
        wr_spikes=wrtds_buster(idx,x.values,res,4.0,figs/f'{sanitize(col)}__07_wrtds_buster.png')

    saqc_flag_mask=None
    if p_flags is not None:
        if hasattr(p_flags,'values'):
            saqc_flag_mask=(p_flags.values.astype(str)!='UNFLAGGED') if p_flags.dtype==object else (p_flags.values!=0)

    final_mask=merge_flags(~np.isfinite(x.values),saqc_flag_mask,wr_spikes.values if hasattr(wr_spikes,'values') else wr_spikes)
    accepted=np.where(final_mask,np.nan,x.values)
    plot_series_basic(idx,accepted,f'{col} Accepted mask',col,figs/f'{sanitize(col)}__08_accepted.png',log=False,legend=False)

    ev_bin=events_from_binary(bin_events,idx)
    ev_flat=events_from_runs(flat_runs,'flat_values')
    ev_slp=events_from_runs(flat_slope_runs,'flat_slopes')
    ev=pd.concat([ev_bin,ev_flat,ev_slp],ignore_index=True) if len(ev_bin)+len(ev_flat)+len(ev_slp)>0 else pd.DataFrame(columns=['type','start','end'])
    evp=tabs/f'{sanitize(col)}__events.csv'
    ev.to_csv(evp,index=False)
    seas.to_csv(tabs/f'{sanitize(col)}__seasonal_summary.csv')

    meta={'station':station,'col':col,'step':str(step),'sentinel_used':sorted(list(used)),'duplicates':int(dup_count),'wrtds_ok':bool(wrtds_q_col is not None and wr_spikes is not None)}
    outdf=pd.DataFrame({f'{sanitize(col)}__raw':s0.reindex(idx).values,f'{sanitize(col)}__clean':x.values,f'{sanitize(col)}__accepted':accepted},index=idx)
    if p_flags is not None:outdf[f'{sanitize(col)}__saqc_flag']=np.array(saqc_flag_mask,dtype=bool)
    if sm_df is not None:
        outdf[f'{sanitize(col)}__sm_masked']=sm_df['value_masked'].reindex(idx).values
        outdf[f'{sanitize(col)}__sm_flagged']=sm_df['is_flagged'].reindex(idx).values
        outdf[f'{sanitize(col)}__sm_flagreason']=sm_df['flag_reason'].reindex(idx).astype(str).values

    return outdf,meta,seas,ev

""" DRIVER """
def run_pipeline(csv_path,station_string,columns,base_out,station_field=None,timestamp_col='timestamp',range_map=None,wrtds_q_col=None,apply_saqc=True,apply_sentem=True,gap_hours=2.0):
    df=read_data(csv_path,timestamp_col=timestamp_col,station_field=station_field,station_value=station_string,parse_dates=True)
    outs=build_structure(base_out,station_string)
    all_out=[];all_meta=[];all_seas=[];all_ev=[]
    for col in columns:
        if col not in df.columns:continue
        rmin=None;rmax=None
        if range_map and col in range_map:
            rmin,rmax=range_map[col]
        o,m,ss,ev=process_variable(df,col,station_string,outs,range_min=rmin,range_max=rmax,wrtds_q_col=wrtds_q_col,apply_saqc=apply_saqc,apply_sentem=apply_sentem,gap_hours=gap_hours)
        all_out.append(o);all_meta.append(m);ss=ss.copy();ss['col']=col;all_seas.append(ss);ev=ev.copy();ev['col']=col;all_ev.append(ev)
    if len(all_out)==0:return
    out_big=pd.concat(all_out,axis=1)
    out_path=outs['processed']/f'{sanitize(station_string)}__qc_timeseries.parquet'
    out_big.to_parquet(out_path,index=True)
    meta=pd.DataFrame(all_meta)
    meta.to_csv(outs['tabs']/f'{sanitize(station_string)}__meta.csv',index=False)
    if len(all_seas)>0:
        seas_big=pd.concat(all_seas,axis=0).reset_index().rename(columns={'index':'season'})
        seas_big.to_csv(outs['tabs']/f'{sanitize(station_string)}__seasonal_all.csv',index=False)
    if len(all_ev)>0:
        ev_big=pd.concat(all_ev,axis=0).reset_index(drop=True)
        ev_big.to_csv(outs['tabs']/f'{sanitize(station_string)}__events_all.csv',index=False)
    return {'timeseries':out_path,'meta':outs['tabs']/f'{sanitize(station_string)}__meta.csv'}

""" MAIN """
if __name__=='__main__':
    base_out='./water_qc_output'
    csv_path='../Erkenruh_Einruhr.csv'
    timestamp_col='timestamp'
    station_field=None
    station_string='Erkenruh_Einruhr'
    columns=[
        'SurfaceWaterConcentration_O2 [mg*L-1]',
        'SurfaceWaterpH [pH]',
        'SurfaceWaterTurbidity [NTU]',
        'SurfaceWaterConcentration_NO3_Trios [mg*L-1]',
        'SurfaceWaterConcentration_NO3_YSI [mg*L-1]'
    ]
    range_map={
        'SurfaceWaterConcentration_O2 [mg*L-1]':(0.0,40.0),
        'SurfaceWaterpH [pH]':(0.0,13.0),
        'SurfaceWaterTurbidity [NTU]':(0.0,4000.0),
        'SurfaceWaterConcentration_NO3_Trios [mg*L-1]':(0.0,35.0),
        'SurfaceWaterConcentration_NO3_YSI [mg*L-1]':(0.0,35.0)
    }
    wrtds_q_col=None
    apply_saqc=True
    apply_sentem=True
    gap_hours=2.0
    out=run_pipeline(csv_path,station_string,columns,base_out,station_field=station_field,timestamp_col=timestamp_col,range_map=range_map,wrtds_q_col=wrtds_q_col,apply_saqc=apply_saqc,apply_sentem=apply_sentem,gap_hours=gap_hours)
    print(out)
