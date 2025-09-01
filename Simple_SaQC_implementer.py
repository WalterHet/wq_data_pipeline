#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 15:59:43 2025

@author: walterh
"""



import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "/home/walterh/thesis/2025_07_16_data_extention/data/Erkenruh_Einuhr/combined/Erkenruh_Einruhr.csv",
    parse_dates=["timestamp"]
)

df = df.set_index("timestamp")

# --------- CONFIG: set your time range here ---------
start_time = "2015-05-01"
end_time   = "2016-07-31"
df = df.loc[start_time:end_time]
# ---------------------------------------------------

for value in df.columns:
    if pd.api.types.is_numeric_dtype(df[value]):
        plt.figure(figsize=(10, 4))
        plt.plot(df.index, df[value], marker="", linestyle="-")
        plt.title(value)
        plt.xlabel("Time")
        plt.ylabel(value)
        plt.tight_layout()
        plt.show()
        
        
# classify all gaps aboce 2 hours as nan data
# detect possible nan value eg -9, -9999, -999 only if they are frequent
# detect binary switches: from value to 0, then from o to value

# Cretae the following as a heatmap of percentage of coverage in percent per seasona. fo one station just a string 
# report nan values per season in percent
# report flat vlaues per season above 2 hours in percent
# report flat slopes per season above 2 hours


#Check if the decimals are equally distirbuted to avoid eg 0.1, 0.2
#Check for douled dates
#Check for long unique doulbe decimals


#Check for possible wrtds, devlop a wrtd buster
#Check for filter or gaussian filter

        
#Please a
        
        
import pandas as pd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import saqc
df=pd.read_csv("/home/walterh/thesis/2025_07_16_data_extention/data/Erkenruh_Einuhr/combined/Erkenruh_Einruhr.csv",parse_dates=["timestamp"]).set_index("timestamp")
start_time="2010-05-01"
end_time="2020-07-31"
df=df.loc[start_time:end_time]
num_cols=[c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
step=int(df.index.to_series().diff().dropna().dt.total_seconds().round().mode().iloc[0])
df_reg=df.resample(f"{step}S").asfreq()
qc_reg=saqc.SaQC(df_reg)
for col in num_cols:
    qc_reg=qc_reg.flagConstants(col,thresh=0.0,window="3H",target=f"{col}__const")
    qc_reg=qc_reg.flagByVariance(col,window="24H",thresh=1e-8,target=f"{col}__lvar")
def plateau_mask(series,eps,min_len):
    b=series.diff().abs().fillna(0)<=eps
    g=(b!=b.shift()).cumsum()
    run=b.groupby(g).transform("sum")>=min_len
    return b&run
def rel_eps(s):
    m=float(s.abs().median()) if np.isfinite(s.abs().median()) else 1.0
    return 1e-6*max(1.0,m)
min_plateau_seconds=6*3600
min_len=max(1,int(round(min_plateau_seconds/step)))
for col in num_cols:
    s=df[col]
    s_reg=df_reg[col]
    flags_const=qc_reg.flags.get(f"{col}__const",pd.Series(index=df_reg.index,dtype=float)).reindex(df.index)
    flags_lvar=qc_reg.flags.get(f"{col}__lvar",pd.Series(index=df_reg.index,dtype=float)).reindex(df.index)
    m_const=flags_const.eq(255.0)
    m_lvar=flags_lvar.eq(255.0)
    eps=rel_eps(s_reg)
    m_plat_reg=plateau_mask(s_reg,eps,min_len)
    m_plat=m_plat_reg.reindex(df.index).fillna(False)
    y_const=s.mask(m_const)
    y_lvar=s.mask(m_lvar)
    y_plat=s.mask(m_plat)
    y_any=s.mask(m_const|m_lvar|m_plat)
    plt.figure(figsize=(12,4))
    plt.plot(s.index,s,linewidth=1)
    plt.plot(y_const.index,y_const,linewidth=1)
    plt.title(f"{col} (constant masked)")
    plt.xlabel("Time")
    plt.ylabel(col)
    plt.tight_layout()
    plt.show()
    plt.figure(figsize=(12,4))
    plt.plot(s.index,s,linewidth=1)
    plt.plot(y_lvar.index,y_lvar,linewidth=1)
    plt.title(f"{col} (low-variance masked)")
    plt.xlabel("Time")
    plt.ylabel(col)
    plt.tight_layout()
    plt.show()
    plt.figure(figsize=(12,4))
    plt.plot(s.index,s,linewidth=1)
    plt.plot(y_plat.index,y_plat,linewidth=1)
    plt.title(f"{col} (plateau masked)")
    plt.xlabel("Time")
    plt.ylabel(col)
    plt.tight_layout()
    plt.show()
    plt.figure(figsize=(12,4))
    plt.plot(s.index,s,linewidth=1)
    plt.plot(y_any.index,y_any,linewidth=1)
    plt.title(f"{col} (any flat/plateau masked)")
    plt.xlabel("Time")
    plt.ylabel(col)
    plt.tight_layout()
    plt.show()

        

