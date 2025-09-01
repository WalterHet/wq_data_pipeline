import re,os,math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
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
COLOR_QCBAND_GREY='#777777'
MINFRAC=0.3
GAP_FACTOR=10
YFLOOR_PERCENTILE=0.5

"""""" UTILS """"""
def sanitize(s):
    return re.sub(r'[^A-Za-z0-9._-]+','_',str(s)).strip('_')
def rolling_mean_std_basic(x,window,center=True,median=False,minfrac=MINFRAC):
    minp=max(3,int(np.ceil(window*minfrac)))
    if median:ma=x.rolling(window=window,center=center,min_periods=minp).median()
    else:ma=x.rolling(window=window,center=center,min_periods=minp).mean()
    sd=x.rolling(window=window,center=center,min_periods=minp).std(ddof=0)
    return ma,sd
def rolling_mean_std_tri5(x,mincount=3):
    v=x.values.astype(float);n=len(v);ma=np.full(n,np.nan);sd=np.full(n,np.nan)
    for i in range(n):
        a=max(0,i-2);b=min(n-1,i+2);seg=v[a:b+1];msk=np.isfinite(seg)
        if msk.sum()<mincount:continue
        w=np.array([1,2,3,2,1],dtype=float)[max(0,2-(i-a)):min(5,3+(b-i))]
        w=w[:msk.size][msk].astype(float);w=w/np.sum(w);seg=seg[msk]
        m=np.sum(w*seg);s=np.sqrt(np.sum(w*(seg-m)**2));ma[i]=m;sd[i]=s
    return pd.Series(ma,index=x.index),pd.Series(sd,index=x.index)
def apply_sensor_error_flags(code,series):
    f=pd.Series(100,index=series.index,dtype=int)
    if code in (157787,):f.loc[series<0]=200
    if code in (2477034,2477787):f.loc[series<0]=200
    if code in (400,400100011):
        f.loc[(series==30)|(series==25)|(series==0)|(series>=40)]=200
        f.loc[~np.isfinite(series)]=200
    if code==410:
        f.loc[(series==0)|(series>=13)|(series<0)]=200
        f.loc[~np.isfinite(series)]=200
    return f
def combine_flags(df_flags):
    g=pd.Series(100,index=df_flags.index,dtype=int)
    for c in df_flags.columns:g.loc[df_flags[c]==200]=200
    return g
def reasons_from_flags(df_flags):
    r=pd.Series("",index=df_flags.index,dtype=object)
    for c in df_flags.columns:
        m=(df_flags[c]==200)
        if m.any():r.loc[m]=r.loc[m].astype(str)+c+" "
    return r.str.strip()
def infer_base_step_ns(idx):
    if len(idx)<3:return np.int64(15*60*1e9)
    d=np.diff(idx.view('int64'));d=d[d>0]
    if d.size==0:return np.int64(15*60*1e9)
    return np.int64(np.median(d))
def break_on_gaps(idx,y,gap_ns):
    y=np.array(y, dtype=float);n=len(y)
    if n<=1:return y
    d=np.diff(idx.view('int64'));out=y.copy()
    for i in range(1,n):
        if d[i-1]>gap_ns:out[i]=np.nan
    return out

"""""" CONFIG """"""
SENTEMQC_CONFIG={157787:{"w1":960,"sf1":3.5,"c1":True,"ta1":10.0,"bs1":1.5,"w2":960,"sf2":3.5,"c2":True,"ta2":10.0,"bs2":1.5,"w3":48,"sf3":1.7,"c3":True,"ta3":2.5,"bs3":2.5,"w4":5,"sf4":1.3,"c4":True,"ta4":0.5,"bs4":0.3,"w5":5,"sf5":1.3,"c5":True,"ta5":0.5,"bs5":0.3,"uncertainty_pct":0.05,"tri5":True},2477034:{"w1":960,"sf1":3.0,"c1":True,"ta1":1.5,"bs1":0.05,"w2":960,"sf2":2.6,"c2":True,"ta2":1.2,"bs2":0.05,"w3":48,"sf3":2.5,"c3":True,"ta3":0.35,"bs3":0.2,"w4":5,"sf4":0.3,"c4":True,"ta4":0.05,"bs4":0.025,"w5":5,"sf5":0.9,"c5":True,"ta5":0.05,"bs5":0.025,"uncertainty_pct":0.03,"tri5":True},2477787:{"w1":960,"sf1":2.5,"c1":True,"ta1":0.4,"bs1":0.05,"w2":960,"sf2":1.5,"c2":True,"ta2":0.3,"bs2":0.05,"w3":48,"sf3":1.0,"c3":True,"ta3":0.35,"bs3":0.01,"w4":5,"sf4":0.2,"c4":True,"ta4":0.03,"bs4":0.01,"w5":5,"sf5":0.9,"c5":True,"ta5":0.05,"bs5":0.025,"uncertainty_pct":0.03,"tri5":True},400:{"w1":960,"sf1":2.5,"c1":True,"ta1":0.4,"bs1":0.2,"w2":960,"sf2":2.5,"c2":True,"ta2":0.4,"bs2":0.2,"w3":48,"sf3":1.7,"c3":True,"ta3":0.35,"bs3":0.35,"w4":5,"sf4":1.3,"c4":True,"ta4":0.05,"bs4":0.05,"w5":5,"sf5":1.3,"c5":True,"ta5":0.05,"bs5":0.05,"uncertainty_pct":0.02,"tri5":True},410:{"w1":960,"sf1":2.0,"c1":True,"ta1":0.3,"bs1":0.2,"w2":960,"sf2":2.0,"c2":True,"ta2":0.3,"bs2":0.2,"w3":12,"sf3":1.6,"c3":True,"ta3":0.05,"bs3":0.05,"w4":5,"sf4":1.3,"c4":True,"ta4":0.05,"bs4":0.05,"w5":5,"sf5":1.3,"c5":True,"ta5":0.05,"bs5":0.05,"uncertainty_con":0.1,"tri5":True}}
CALIBRATION_OFFSETS={"SurfaceWaterConcentration_O2 [mg*L-1]":0.0,"SurfaceWaterpH [pH]":0.0,"SurfaceWaterTurbidity [NTU]":0.0,"SurfaceWaterConcentration_NO3_Trios [mg*L-1]":0.0,"SurfaceWaterConcentration_NO3_YSI [mg*L-1]":0.0}
VARIABLE_MAP=[{"col":"SurfaceWaterConcentration_O2 [mg*L-1]","code":400,"unit":"mg/L","label":"Dissolved Oxygen","is_nitrate":False},{"col":"SurfaceWaterpH [pH]","code":410,"unit":"","label":"pH","is_nitrate":False},{"col":"SurfaceWaterTurbidity [NTU]","code":157787,"unit":"NTU","label":"Turbidity","is_nitrate":False},{"col":"SurfaceWaterConcentration_NO3_Trios [mg*L-1]","code":2477034,"unit":"mg/L","label":"NO3 Trios","is_nitrate":True},{"col":"SurfaceWaterConcentration_NO3_YSI [mg*L-1]","code":2477787,"unit":"mg/L","label":"NO3 YSI","is_nitrate":True}]

"""""" CORE QC """"""
def apply_sentemqc_to_series(ts,code,config,is_nitrate):
    ts_in=pd.Series(ts.astype(float))
    mask_pre=(~np.isfinite(ts_in))|(ts_in<=0)
    s=ts_in.mask(mask_pre,np.nan)
    df=pd.DataFrame({"OBS_in":ts_in,"OBS_raw":s})
    if is_nitrate:df["OBS"]=df["OBS_raw"]*(14.0/62.0)
    else:df["OBS"]=df["OBS_raw"]
    df["flag_sensor"]=apply_sensor_error_flags(code,df["OBS"])
    flags=["flag_sensor"]
    for run in [1,2,3,4,5]:
        if run==5 and config.get("tri5",False):ma,sd=rolling_mean_std_tri5(df["OBS"])
        else:ma,sd=rolling_mean_std_basic(df["OBS"],window=config["w"+str(run)],center=config["c"+str(run)],median=(run==3),minfrac=MINFRAC)
        top=ma+config["ta"+str(run)]+config["sf"+str(run)]*sd
        bot=ma-config["bs"+str(run)]-config["sf"+str(run)]*sd
        fu=pd.Series(100,index=df.index,dtype=int)
        if "uncertainty_pct" in config:
            topunc=(1.0+config["uncertainty_pct"])*df["OBS"];botunc=(1.0-config["uncertainty_pct"])*df["OBS"]
            fu.loc[((topunc>top)&(botunc>top))|((botunc<bot)&(topunc<bot))]=200
        elif "uncertainty_con" in config:
            topunc=df["OBS"]+config["uncertainty_con"];botunc=df["OBS"]-config["uncertainty_con"]
            fu.loc[((topunc>top)&(botunc>top))|((botunc<bot)&(topunc<bot))]=200
        else:
            fu.loc[(df["OBS"]>top)|(df["OBS"]<bot)]=200
        df["flag_run"+str(run)]=fu;flags.append("flag_run"+str(run))
        if run==5:
            valid=df["OBS"].notna().astype(int);minp=max(3,int(np.ceil(config["w5"]*MINFRAC)))
            support=valid.rolling(window=config["w5"],center=config["c5"],min_periods=1).sum()
            top=top.where(support>=minp);bot=bot.where(support>=minp)
            df["qcband_top"]=top;df["qcband_bottom"]=bot
    df["flag_global"]=combine_flags(df[flags]);df["is_flagged"]=df["flag_global"]==200;df["flag_reason"]=reasons_from_flags(df[flags])
    if is_nitrate:
        df["qcband_top"]=df["qcband_top"]*(62.0/14.0);df["qcband_bottom"]=df["qcband_bottom"]*(62.0/14.0)
    df["value_masked"]=np.where(df["is_flagged"],np.nan,df["OBS_raw"])
    df["mask_pre"]=mask_pre.values
    return df[["OBS_in","OBS_raw","value_masked","flag_global","is_flagged","flag_reason","qcband_top","qcband_bottom","mask_pre"]]

"""""" PLOTTER """"""
def plot_combined_log_gap(idx,raw,is_flagged,top,bot,title,ylabel,mask_pct,avail_pct,flag_pct,savepath=None):
    raw=np.asarray(raw,dtype=float);is_flagged=np.asarray(is_flagged,dtype=bool)
    acc=np.where((~is_flagged)&np.isfinite(raw)&(raw>0),raw,np.nan);flg=np.where(is_flagged&np.isfinite(raw)&(raw>0),raw,np.nan)
    t=np.asarray(top,dtype=float);b=np.asarray(bot,dtype=float);t=np.where(t>0,t,np.nan);b=np.where(b>0,b,np.nan)
    step_ns=infer_base_step_ns(idx.values);gap_ns=np.int64(GAP_FACTOR*step_ns)
    raw_line=break_on_gaps(idx.values,raw,gap_ns);t_line=break_on_gaps(idx.values,t,gap_ns);b_line=break_on_gaps(idx.values,b,gap_ns)
    fig=plt.figure(figsize=(12,5))
    plt.plot(idx,raw_line,color=COLOR_RAW,lw=0.6,label="Raw")
    if np.isfinite(t_line).any() and np.isfinite(b_line).any():
        plt.plot(idx,t_line,color=COLOR_QCBAND_GREY,lw=0.8,label="QC band top")
        plt.plot(idx,b_line,color=COLOR_QCBAND_GREY,lw=0.8,label="QC band bottom")
    if np.isfinite(acc).any():plt.plot(idx,acc,'.',ms=0.6,color=COLOR_ORANGE,label="Accepted")
    if np.isfinite(flg).any():plt.plot(idx,flg,'x',ms=0.6,color=COLOR_BLUE,label="Flagged")
    pos_vals=np.concatenate([v[np.isfinite(v)&(v>0)] for v in [raw_line,acc,flg,t_line,b_line] if v is not None]) if any([raw_line is not None]) else np.array([])
    if pos_vals.size>0:
        base=np.nanpercentile(pos_vals,max(1e-9,min(100.0,max(0.0,YFLOOR_PERCENTILE))))
        ymin=max(1e-8,base*0.5)
        plt.ylim(bottom=ymin)
    plt.yscale("log");plt.title(title);plt.xlabel("Time");plt.ylabel(ylabel)
    ax=plt.gca();loc=mdates.AutoDateLocator();ax.xaxis.set_major_locator(loc);ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
    plt.legend(ncol=5,frameon=False,loc="upper right")
    txt=f"Availability: {avail_pct:.1f}%   Masked≤0: {mask_pct:.1f}%   Flagged: {flag_pct:.1f}%"
    plt.text(0.01,0.02,txt,transform=ax.transAxes,ha="left",va="bottom")
    plt.tight_layout()
    if savepath is not None:plt.savefig(savepath,bbox_inches='tight')
    plt.close(fig)

"""""" RUNNERS """"""
def run_sentemqc_on_dataframe(df,variable_map=VARIABLE_MAP,calibration_offsets=CALIBRATION_OFFSETS,config=SENTEMQC_CONFIG,start_time=None,end_time=None,plot_dir=None):
    d=df.copy()
    if start_time is not None and end_time is not None:
        d=d.loc[start_time:end_time]
    for k,v in calibration_offsets.items():
        if k in d.columns:d[k]=d[k].astype(float)+float(v)
    results={}
    for entry in variable_map:
        c=entry["col"]
        if c not in d.columns:continue
        code=entry["code"];cfg=config[code];is_nitrate=entry["is_nitrate"]
        series=d[c].astype(float)
        res=apply_sentemqc_to_series(series,code,cfg,is_nitrate)
        base=sanitize(c);results[base]=res
        d[base+"__flag_global"]=res["flag_global"].astype(int)
        d[base+"__is_flagged"]=res["is_flagged"].astype(bool)
        d[base+"__flag_reason"]=res["flag_reason"].astype(str)
        d[base+"__qcband_top"]=res["qcband_top"].astype(float)
        d[base+"__qcband_bottom"]=res["qcband_bottom"].astype(float)
        d[base+"__masked"]=res["value_masked"].astype(float)
        if plot_dir is not None:
            lab=entry["label"];unit=entry["unit"]
            N=len(res);avail_pct=100.0*np.isfinite(res["OBS_in"].values).sum()/max(1,N)
            mask_pct=100.0*(res["mask_pre"].values.astype(bool).sum())/max(1,N)
            flag_pct=100.0*((res["is_flagged"].values & np.isfinite(res["OBS_raw"].values)).sum())/max(1,N)
            Path(plot_dir).mkdir(parents=True,exist_ok=True)
            plot_combined_log_gap(res.index,res["OBS_raw"].values,res["is_flagged"].values,res["qcband_top"].values,res["qcband_bottom"].values,f"{lab} — QC results (orange=accepted, blue=flagged, log-scale)",f"{lab} [{unit}]" if unit else lab,mask_pct,avail_pct,flag_pct,savepath=Path(plot_dir)/f"{sanitize(c)}__sentemqc.png")
    return d,results
def run_sentemqc_on_file(base_csv,start_time=None,end_time=None,plot_dir=None):
    df=pd.read_csv(base_csv,parse_dates=["timestamp"]).set_index("timestamp").sort_index()
    return run_sentemqc_on_dataframe(df,start_time=start_time,end_time=end_time,plot_dir=plot_dir)

"""""" MAIN """"""
if __name__=='__main__':
    base_csv="./Erkenruh_Einruhr.csv"
    plot_dir="./water_qc_output/reports/figures/Erkenruh_Einruhr/sentemqc"
    _df,_res=run_sentemqc_on_file(base_csv,plot_dir=plot_dir)
    out_csv="./water_qc_output/processed/Erkenruh_Einruhr/sentemqc_export.csv"
    Path(Path(out_csv).parent).mkdir(parents=True,exist_ok=True)
    _df.to_csv(out_csv)
    print({"export":out_csv,"plots":str(plot_dir)})

