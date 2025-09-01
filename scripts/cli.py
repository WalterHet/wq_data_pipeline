import argparse,sys,os
from pathlib import Path
from scripts.ufz_water_qc_pipeline import run_pipeline
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--csv",required=True)
    p.add_argument("--station",required=True)
    p.add_argument("--out",default="water_qc_output")
    p.add_argument("--vars",nargs="+",required=True)
    p.add_argument("--wrtds-q-col",default=None)
    a=p.parse_args()
    csv_path=Path(a.csv).resolve()
    base_out=Path(a.out).resolve()
    columns=a.vars
    range_map={}
    run_pipeline(csv_path=csv_path,station_string=a.station,columns=columns,range_map=range_map,base_out=base_out,wrtds_q_col=a.wrtds_q_col)
if __name__=="__main__":
    sys.exit(main())

