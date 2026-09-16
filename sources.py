"""Public upstream adapters. Every returned frame retains its exact run and member."""
import datetime as dt, json, pathlib, re, urllib.parse, urllib.request, time, os
from concurrent.futures import ThreadPoolExecutor
CACHE=pathlib.Path(os.environ.get('MODELSCOPE_CACHE_DIR','/tmp/modelscope-grib-cache')); CACHE.mkdir(parents=True,exist_ok=True)
LEVELS=[1000,975,950,925,900,875,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100]
VARS='HGT TMP DPT RH SPFH UGRD VGRD PRES PRMSL MSLET MSLMA GUST CAPE CIN HLCY REFD REFC APCP TCDC LCDC MCDC HCDC PWAT VIS CSNOW CRAIN CFRZR CICEP LFTX MASSDEN COLMD MXUPHL UPHL LTNG SBT113 SBT114 SBT123 SBT124'.split()

def read(url,byte_range=None):
 headers={'User-Agent':'Modelscope/2.0'}
 if byte_range: headers['Range']=f'bytes={byte_range[0]}-{byte_range[1]}'
 for attempt in range(3):
  try:
   with urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=120) as r:
    if byte_range and r.status!=206:raise ValueError('Upstream did not honor byte-range request')
    body=r.read(450_000_001)
    if len(body)>450_000_000:raise ValueError('Response exceeds download limit')
    if byte_range and len(body)!=byte_range[1]-byte_range[0]+1:raise ValueError('Incomplete range')
    return body
  except Exception:
   if attempt==2:raise
   time.sleep(2+attempt*2)

def noaa_url(model,run,hour,member='ctl'):
 date,cycle=run.strftime('%Y%m%d'),run.strftime('%H')
 if model=='gfs':script,file,directory='gfs_0p25',f'gfs.t{cycle}z.pgrb2.0p25.f{hour:03}',f'/gfs.{date}/{cycle}/atmos'
 elif model=='hrrr':script,file,directory='hrrr_2d',f'hrrr.t{cycle}z.wrfprsf{hour:02}.grib2',f'/hrrr.{date}/conus'
 elif model=='nam':script,file,directory='nam',f'nam.t{cycle}z.awphys{hour:02}.tm00.grib2',f'/nam.{date}'
 elif model=='rap':script,file,directory='rap',f'rap.t{cycle}z.awp130pgrbf{hour:02}.grib2',f'/rap.{date}'
 elif model=='hiresw':script,file,directory='hiresw_conus',f'hiresw.t{cycle}z.arw_5km.f{hour:02}.conus.grib2',f'/hiresw.{date}'
 else:raise ValueError('Unsupported filter model')
 params={'file':file,'dir':directory,'subregion':'','leftlon':-125,'rightlon':-65,'toplat':50,'bottomlat':24}
 params.update({f'var_{v}':'on' for v in VARS})
 # Include parcel layers and special fields; select all levels to avoid silently missing diagnostics.
 params['all_lev']='on'
 return f'https://nomads.ncep.noaa.gov/cgi-bin/filter_{script}.pl?'+urllib.parse.urlencode(params)

def selected_noaa(line):
 parts=line.split(':')
 if len(parts)<6:return False
 var,level=parts[3:5]
 if var not in VARS:return False
 match=re.fullmatch(r'(\d+) mb',level)
 if match:return int(match[1]) in LEVELS and var in ['HGT','TMP','DPT','RH','SPFH','UGRD','VGRD']
 return not any(s in level for s in ['below ground','hybrid level','sigma'])

def indexed(url,ecmwf=False,member=None,selector=None,output=None):
 rows=[]
 if ecmwf:
  for line in read(url.removesuffix('.grib2')+'.index').splitlines():
   r=json.loads(line)
   if member is not None and str(r.get('number'))!=str(member):continue
   if r.get('levtype')=='pl':ok=r.get('param') in ['t','r','q','u','v','gh','z'] and int(r.get('levelist',0)) in LEVELS
   else:ok=r.get('param') in ['2t','2d','10u','10v','10fg','msl','sp','tp','tcc','tcwv','mucape','ptype','z']
   if selector is not None:ok=selector(r)
   if ok:rows.append((r['_offset'],r['_offset']+r['_length']-1))
 else:
  lines=read(url+'.idx').decode().splitlines()
  for i,line in enumerate(lines):
   if (selector(line) if selector else selected_noaa(line)):
    end=int(lines[i+1].split(':')[1])-1 if i+1<len(lines) else None
    if end is None:
     # Last record length is encoded in its GRIB2 header.
     start=int(line.split(':')[1]);head=read(url,(start,start+15));end=start+int.from_bytes(head[8:16],'big')-1
    rows.append((int(line.split(':')[1]),end))
 if not rows:raise ValueError('No requested records in upstream inventory')
 ranges=[]
 for start,end in sorted(rows):
  if ranges and start==ranges[-1][1]+1:ranges[-1]=(ranges[-1][0],end)
  else:ranges.append((start,end))
 if output is not None:
  # Stream one range at a time to bound memory on small ingestion workers.
  with open(output,'wb') as out:
   for byte_range in ranges:out.write(read(url,byte_range))
  return output
 with ThreadPoolExecutor(max_workers=4) as pool: parts=list(pool.map(lambda r:read(url,r),ranges))
 return b''.join(parts)

def source_url(model,run,hour,member='det'):
 date,cycle=run.strftime('%Y%m%d'),run.strftime('%H')
 if model in ['ecmwf','eps']:
  stream,typ=('oper','fc') if model=='ecmwf' else ('enfo','ef')
  source=f'https://data.ecmwf.int/forecasts/{date}/{cycle}z/ifs/0p25/{stream}/{date}{cycle}0000-{hour}h-{stream}-{typ}.grib2'
 elif model=='hrrr':
  source=f'https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod/hrrr.{date}/conus/hrrr.t{cycle}z.wrfprsf{hour:02}.grib2'
 elif model=='sref':
  core,pert=member.split('-');source=f'https://nomads.ncep.noaa.gov/pub/data/nccf/com/sref/prod/sref.{date}/{cycle}/pgrb/sref_{core}.t{cycle}z.pgrb212.{pert}.f{hour:02}.grib2'
 elif model=='rrfs':
  source=f'https://nomads.ncep.noaa.gov/pub/data/nccf/com/rrfs/para/rrfs.{date}/{cycle}/rrfs.t{cycle}z.prslev.3km.f{hour:03}.conus.grib2'
 elif model=='gfs':source=f'https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.{date}/{cycle}/atmos/gfs.t{cycle}z.pgrb2.0p25.f{hour:03}'
 elif model=='nam':source=f'https://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod/nam.{date}/nam.t{cycle}z.awphys{hour:02}.tm00.grib2'
 elif model=='rap':source=f'https://nomads.ncep.noaa.gov/pub/data/nccf/com/rap/prod/rap.{date}/rap.t{cycle}z.awp130pgrbf{hour:02}.grib2'
 elif model=='hiresw':source=f'https://nomads.ncep.noaa.gov/pub/data/nccf/com/hiresw/prod/hiresw.{date}/hiresw.t{cycle}z.arw_5km.f{hour:02}.conus.grib2'
 else:raise ValueError('Unknown model')
 return source

def fetch(model,run,hour,member='det'):
 path=CACHE/f'{model}-{run:%Y%m%d%H}-{hour:03}-{member}.grib2'
 source=source_url(model,run,hour,member)
 if not path.exists():
  temp=path.with_suffix('.tmp')
  indexed(source,model in ['ecmwf','eps'],member if model=='eps' else None,output=temp)
  with temp.open('rb') as f:
   if f.read(4)!=b'GRIB':raise ValueError('Upstream did not return GRIB data')
  temp.replace(path)
 return path,source
