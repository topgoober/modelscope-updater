"""Exact accumulation windows and convective maxima from complete hourly histories."""
import datetime as dt,json,gzip,re,pathlib
import numpy as np
from sources import source_url,indexed,CACHE
from ingest import decode,save_frame,update_catalog,serial
ROOT=pathlib.Path(__file__).parent

def temporal(model,run,hour,member):
 source=source_url(model,run,hour,member)
 if model=='hrrr':source=source.replace('wrfprs','wrfsfc')
 if model=='rrfs':source=source.replace('.prslev.','.2dfld.')
 p=CACHE/f'history-{model}-{run:%Y%m%d%H}-{hour}-{member}.grib2'
 if not p.exists():
  if model in ['ecmwf','eps']:selector=lambda r:r.get('param')=='tp'
  else:selector=lambda line: bool(re.search(r':(APCP|MXUPHL|MSLMA|PRMSL|MSLET):',line))
  indexed(source,model in ['ecmwf','eps'],member if model=='eps' else None,selector,output=p)
 fields,meta,_,_=decode(p)
 import os
 if os.environ.get('MODELSCOPE_EPHEMERAL_HISTORY')=='1':p.unlink(missing_ok=True)
 # Save the 0–3 km 1-hour max internally as well as public products.
 import eccodes as ec
 # decode() exposes exact hourly maxima; no max is inferred from instantaneous UH.
 return fields,meta

def augment(frame,full_maxima=True):
 model,member,hour=frame['model'],frame['member'],frame['hour'];run=dt.datetime.fromisoformat(frame['run'].replace('Z',''))
 frame['fields']={k:np.array(v,dtype=float) for k,v in frame['fields'].items()}
 needed={hour}
 for window in [1,6,12,24]:
  if hour>=window and (model not in ['sref','ecmwf','eps'] or (hour-window)%3==0):needed.add(hour-window)
 if full_maxima and model in ['hrrr','hiresw','rrfs']:needed.update(range(1,hour+1))
 history={};errors=[]
 for h in sorted(needed):
  try:history[h]=temporal(model,run,h,member)
  except Exception as e:errors.append(f'F{h:03}: {e}')
 fields=frame['fields'];meta=frame['fieldMeta'];size=frame['grid']['width']*frame['grid']['height']
 for h,(f,m) in history.items():
  if h==hour:
   for key in ['mslp','qpf1','qpf6','qpf12','qpf24','qpf_total','uh25_1']:
    if key in f:fields[key]=f[key];meta[key]=m.get(key,dict(kind='native'))
 for window in [1,6,12,24]:
  if hour<window:continue
  current=history.get(hour,({},{}))[0].get('qpf_total',fields.get('qpf_total'))
  prior=np.zeros(size) if hour==window else history.get(hour-window,({},{}))[0].get('qpf_total')
  if current is not None and prior is not None:
   result=current-prior;result[result<-.1]=np.nan;result=np.maximum(result,0);fields[f'qpf{window}']=result;meta[f'qpf{window}']=dict(kind='derived',method='Difference of same-run accumulated precipitation',start=hour-window,end=hour)
 if model in ['hrrr','hiresw','rrfs']:
  for key,source_key,steps in [('uh25_3','uh25_1',list(range(max(1,hour-2),hour+1))),('uh25_run','uh25_1',list(range(1,hour+1))),('uh03_run','_uh03_1',list(range(1,hour+1)))]:
   if key=='uh25_3' and hour<3:continue
   if steps and all(h in history and source_key in history[h][0] for h in steps):
    fields[key]=np.maximum.reduce([history[h][0][source_key] for h in steps]);meta[key]=dict(kind='derived',method='Maximum over complete hourly maximum fields',start=hour-3 if key=='uh25_3' else 0,end=hour)
 frame['historyErrors']=errors
 return frame

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--models',nargs='+');p.add_argument('--hour',type=int);p.add_argument('--no-run-max',action='store_true');args=p.parse_args()
 c=json.loads((ROOT/'dist/catalog.json').read_text())
 for e in c['frames']:
  if args.models and e['model'] not in args.models:continue
  if args.hour is not None and e['hour'] not in [args.hour,args.hour+3 if e['model']=='sref' else args.hour]:continue
  print('History',e['model'],e['member'],e['hour'],flush=True)
  f=json.load(gzip.open(ROOT/'dist'/e['file']));f=augment(f,not args.no_run_max);update_catalog([save_frame(f)]);print('Updated',e['model'],'history gaps:',len(f['historyErrors']),flush=True)
