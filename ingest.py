"""Decode real GRIB records into lazy-loadable, provenance-aware forecast frames."""
import datetime as dt, gzip,json,pathlib,warnings
import numpy as np
import eccodes as ec
ec.codes_grib_multi_support_on()
from scipy.spatial import cKDTree
from catalog import MODELS,PRODUCTS
from sources import LEVELS,fetch
from diagnostics import diagnostics,dewpoint,mixing
ROOT=pathlib.Path(__file__).parent
from region_grid import display_grids
GRID,PGRID,PROFILE_STRIDE=display_grids()

def sphere(lon,lat):
 lon,lat=np.deg2rad(lon),np.deg2rad(lat)
 return np.column_stack((np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)))
def coordinates(g):
 x,y=np.meshgrid(g['west']+np.arange(g['width'])*g['step'],g['north']-np.arange(g['height'])*g['step']);return x.ravel(),y.ravel()
COARSE=(np.arange(PGRID['height'])[:,None]*PROFILE_STRIDE*GRID['width']+np.arange(PGRID['width'])[None,:]*PROFILE_STRIDE).ravel()
xx,yy=coordinates(GRID);UPSAMPLE=(np.minimum(np.rint((PGRID['north']-yy)/PGRID['step']).astype(int),PGRID['height']-1)*PGRID['width']+np.minimum(np.rint((xx-PGRID['west'])/PGRID['step']).astype(int),PGRID['width']-1))

def get(h,key,default=None):
 try:return ec.codes_get(h,key)
 except Exception:return default

def decode(path):
 flds={};meta={};profiles={};samplers={};wind_rotation={};inventory=[];accum=[];maxima=[]
 with open(path,'rb') as stream:
  while (h:=ec.codes_grib_new_from_file(stream)) is not None:
   try:
    name=get(h,'shortName');typ=get(h,'typeOfLevel');lev=get(h,'level');units=get(h,'units');start=get(h,'startStep',0);end=get(h,'endStep',0)
    key=None;scale=1.;offset=0.;profile_key=None
    if typ=='isobaricInhPa' and lev in LEVELS:
     k={'t':'t','r':'rh','q':'q','u':'u','v':'v','gh':'z','z':'z','dpt':'td'}.get(name)
     if k:profile_key=(int(lev),k);offset=-273.15 if k in ['t','td'] else 0;scale=1/9.80665 if name=='z' else 1
    if profile_key is None:
     surface={('2t',2):'tmp',('2d',2):'dpt',('2r',2):'rh',('10u',10):'u10',('10v',10):'v10'}
     key=surface.get((name,lev))
     if key in ['tmp','dpt']:offset=-273.15
     if name in ['sp','pres'] and typ=='surface':key='ps';scale=.01
     if name in ['orog','gh','z'] and typ=='surface':key='terrain';scale=1/9.80665 if name=='z' else 1
     if name in ['prmsl','msl','mslet','mslma']:key='mslp';scale=.01
     if name in ['gust','10fg']:key='gust';scale=1.943844
     if name in ['cape','cin'] and typ=='surface':key=name
     if name=='mucape':key='mucape_native'
     if name=='cape' and typ=='heightAboveGroundLayer' and set([get(h,'topLevel'),get(h,'bottomLevel')])=={0,3000}:key='cape03_native'
     if name in ['cape','cin'] and typ=='pressureFromGroundLayer':key=f'{name}_layer{int(lev/100)}'
     if name=='refc':key='refc'
     if name=='refd' and typ=='heightAboveGround' and lev==1000:key='ref1'
     if name=='vis':key='visibility';scale=.001
     if (name=='pwat' and typ in ['atmosphere','atmosphereSingleLayer','entireAtmosphere']) or name=='tcwv':key='pwat'
     if name in ['tcc','lcc','mcc','hcc']:key={'tcc':'cloud','lcc':'cloud_low','mcc':'cloud_mid','hcc':'cloud_high'}[name];scale=100 if units in ['(0 - 1)','Proportion','Numeric'] else 1
     if name in ['crain','csnow','cicep','cfrzr']:key=name
     if name=='ptype':key='ptype_ecmwf'
     if name=='lftx':key='li_native'
     if name=='hlcy':key=f'srh{int(lev/1000)}_native'
     if name=='SBT114':key='ir';offset=-273.15
     if name=='ltng':key='lightning_native'
     if get(h,'parameterCategory')==20 and get(h,'parameterNumber')==0 and typ=='heightAboveGround':key='smoke_sfc';scale=1e9
     if get(h,'parameterCategory')==20 and get(h,'parameterNumber')==1 and typ=='atmosphereSingleLayer':key='smoke_col';scale=1e6
     if name=='tp':key='_qpf';scale=1000 if units=='m' else 1
     if get(h,'parameterCategory')==7 and get(h,'parameterNumber')==199 and typ=='heightAboveGroundLayer':key='_uh'
    if not key and not profile_key:continue
    gridkey=get(h,'md5GridSection',str(get(h,'numberOfPoints')))
    if gridkey not in samplers:
     lon=ec.codes_get_array(h,'longitudes');lat=ec.codes_get_array(h,'latitudes');tree=cKDTree(sphere(lon,lat));dist,ix=tree.query(sphere(xx,yy));samplers[gridkey]=(ix,dist>.008)
     if get(h,'uvRelativeToGrid')==1 and get(h,'gridType')=='lambert':
      la1=np.deg2rad(get(h,'Latin1InDegrees'));la2=np.deg2rad(get(h,'Latin2InDegrees'));n=np.sin(la1) if abs(la1-la2)<1e-6 else np.log(np.cos(la1)/np.cos(la2))/np.log(np.tan(np.pi/4+la2/2)/np.tan(np.pi/4+la1/2));angle=n*np.deg2rad(((lon[ix]-get(h,'LoVInDegrees')+180)%360)-180);wind_rotation[gridkey]=angle
     del lon,lat,tree
    ix,mask=samplers[gridkey];values=ec.codes_get_values(h)
    if get(h,'bitmapPresent')==1:
     bitmap=ec.codes_get_array(h,'bitmap');values[bitmap==0]=np.nan
    a=values[ix]*scale+offset;a[mask|(~np.isfinite(a))|(np.abs(a)>1e10)]=np.nan
    if key=='cin' or (key and key.startswith('cin_')):a=-np.abs(a)
    definition=dict(name=get(h,'name'),level=f'{lev} {typ}',units=units,start=start,end=end,kind='native')
    if profile_key:
     lev,k=profile_key;profiles.setdefault(lev,{})[k]=a
     if k in ['u','v']:profiles[lev]['rotation']=wind_rotation.get(gridkey)
    elif key=='_qpf':accum.append(dict(start=start,end=end,values=a))
    elif key=='_uh':maxima.append(dict(start=start,end=end,top=get(h,'topLevel'),bottom=get(h,'bottomLevel'),values=a))
    else:
     flds[key]=a;meta[key]=definition
     if key in ['u10','v10']:flds['_rotation']=wind_rotation.get(gridkey)
    inventory.append(dict(key=key or f'{profile_key[1]}{profile_key[0]}',**definition))
   finally:ec.codes_release(h)
 if isinstance(flds.get('_rotation'),np.ndarray) and 'u10'in flds and 'v10'in flds:
  a=flds['_rotation'];u,v=flds['u10'],flds['v10'];flds['u10'],flds['v10']=u*np.cos(a)+v*np.sin(a),-u*np.sin(a)+v*np.cos(a)
 flds.pop('_rotation',None)
 for lev,pr in profiles.items():
  a=pr.pop('rotation',None)
  if a is not None and 'u'in pr and 'v'in pr:u,v=pr['u'],pr['v'];pr['u'],pr['v']=u*np.cos(a)+v*np.sin(a),-u*np.sin(a)+v*np.cos(a)
  if 'td' not in pr and 't'in pr:
   if 'rh'in pr:
    from diagnostics import es
    vapor=es(pr['t']+273.15)*np.clip(pr['rh']/100,.00001,1);r=.62195691*vapor/(lev-vapor);pr['td']=dewpoint(lev,r)-273.15
   elif 'q'in pr:pr['td']=dewpoint(lev,pr['q']/(1-pr['q']))-273.15
  if 'ps'in flds:
   below=(lev>flds['ps'])
   for k,a in pr.items():a[below]=np.nan
  for k in ['t','z','u','v']:
   if k in pr:flds[f'{k}{lev}']=pr[k];meta[f'{k}{lev}']=dict(kind='native',level=f'{lev} hPa',name=k)
 for item in accum:
  duration=item['end']-item['start']
  if duration in [1,6,12,24]:flds[f'qpf{duration}']=item['values'];meta[f'qpf{duration}']=dict(kind='native',start=item['start'],end=item['end'])
  if item['start']==0:flds['qpf_total']=item['values'];meta['qpf_total']=dict(kind='native',start=0,end=item['end'])
 for item in maxima:
  duration=item['end']-item['start'];bounds={item['top'],item['bottom']}
  k='uh25_1' if duration==1 and bounds=={2000,5000} else 'uh25_3' if duration==3 and bounds=={2000,5000} else 'uh25_run' if item['start']==0 and bounds=={2000,5000} else 'uh03_run' if item['start']==0 and bounds=={0,3000} else None
  if duration==1 and bounds=={0,3000}:flds['_uh03_1']=item['values']
  if k:flds[k]=item['values'];meta[k]=dict(kind='native',start=item['start'],end=item['end'])
 if all(k in flds for k in ['crain','csnow','cicep','cfrzr']):
  ptype=np.zeros(len(xx));valid=np.ones(len(xx),bool)
  for k,code in [('crain',1),('csnow',2),('cicep',3),('cfrzr',4)]:valid&=np.isfinite(flds[k]);ptype=np.where(flds[k]>.5,code,ptype)
  ptype[~valid]=np.nan;flds['ptype']=ptype
 if 'ptype_ecmwf'in flds:
  raw=flds['ptype_ecmwf'];out=np.full(raw.shape,np.nan)
  for codes,value in [([0],0),([1,3],1),([5,6,7,8],2),([4],3),([2],4)]:out[np.isin(raw,codes)]=value
  flds['ptype']=out
 return flds,meta,profiles,inventory

def derive(flds,meta,profiles):
 import metpy.calc as mp
 from metpy.units import units
 if all(k in flds for k in ['tmp','dpt']):
  flds['rh']=np.clip(mp.relative_humidity_from_dewpoint((flds['tmp']+273.15)*units.K,(flds['dpt']+273.15)*units.K).m*100,0,100)
 if all(k in flds for k in ['tmp','rh','u10','v10']):
  flds['apparent']=np.asarray(mp.apparent_temperature((flds['tmp']+273.15)*units.K,flds['rh']/100,np.hypot(flds['u10'],flds['v10'])*units('m/s'),mask_undefined=False).to('degC').m)
 for lev in [700,850]:
  if all(f'{k}{lev}'in flds for k in ['t','u','v']):
   t=flds[f't{lev}'].reshape(GRID['height'],GRID['width']);dy=-GRID['step']*111195;dx=GRID['step']*111195*np.cos(np.deg2rad(yy.reshape(t.shape)));gy,gx=np.gradient(t);flds[f'adv{lev}']=-(flds[f'u{lev}']*(gx/dx).ravel()+flds[f'v{lev}']*(gy/dy).ravel())*3600
 if 'terrain' not in flds and all(k in flds for k in ['ps','tmp','dpt']):
  terrain=np.full(len(xx),np.nan)
  for pp in sorted(profiles,reverse=True):
   pr=profiles[pp]
   if all(k in pr for k in ['z','t','td']):
    tmean=(flds['tmp']+pr['t'])/2+273.15
    rmean=(mixing(flds['ps'],flds['dpt']+273.15)+mixing(pp,pr['td']+273.15))/2
    tv=tmean*(1+rmean/.62195691)/(1+rmean)
    estimate=pr['z']-287.05/9.80665*tv*np.log(flds['ps']/pp)
    terrain=np.where(np.isnan(terrain)&(pp<=flds['ps'])&np.isfinite(estimate),estimate,terrain)
  flds['terrain']=terrain;meta['terrain']=dict(kind='derived',method='Hydrostatic surface-height estimate from the lowest above-ground pressure level and surface pressure.')
 available=sorted([p for p,pr in profiles.items() if all(k in pr for k in ['t','td','z','u','v'])],reverse=True)
 if len(available)<5 or not all(k in flds for k in ['tmp','dpt','ps','terrain','u10','v10']):return None
 pr={k:np.array([profiles[p][k][COARSE] for p in available]) for k in ['t','td','z','u','v']};pr['z']=pr['z']-flds['terrain'][COARSE]
 for k in pr:pr[k][pr['z']<0]=np.nan
 sfc={k:flds[k][COARSE] for k in ['tmp','dpt','ps','terrain','u10','v10']}
 with warnings.catch_warnings():
  warnings.simplefilter('ignore');d,parcel=diagnostics(np.array(available,dtype=float),pr,sfc)
 for k,a in d.items():
  if k in ['cape','cin'] and k in flds:continue
  flds[k]=a[UPSAMPLE];meta[k]=dict(kind='derived',method=f'Pressure-level calculation on a {PGRID["step"]}° sampled profile grid; native model and vertical resolution still apply.')
 for k in ['lclp','storm_u','storm_v']:sfc[k]=d[k]
 return dict(grid=PGRID,pressure=available,**pr,parcel=parcel,surface=sfc)

def clean(a):
 a=np.asarray(a);rounded=np.round(a,2);return np.where(np.isfinite(rounded),rounded, None).tolist()
def serial(obj):
 if isinstance(obj,np.ndarray):return clean(obj)
 if isinstance(obj,dict):return {k:serial(v) for k,v in obj.items()}
 if isinstance(obj,list):return [serial(v) for v in obj]
 return obj

def build_frame(model,run,hour,member='det'):
 path,source=fetch(model,run,hour,member);print(f'Decoding {model} {member} F{hour:03}',flush=True)
 flds,meta,pr,inventory=decode(path)
 if model=='rrfs':
  from sources import indexed
  surface_path=path.with_name(path.stem+'-surface.grib2')
  surface_source=source.replace('.prslev.','.2dfld.')
  if not surface_path.exists():indexed(surface_source,output=surface_path)
  sf,sm,sp,si=decode(surface_path);flds.update(sf);meta.update(sm);inventory.extend(si)
  for lev,values in sp.items():pr.setdefault(lev,{}).update(values)
 profile=derive(flds,meta,pr)
 return dict(model=model,member=member,run=run.isoformat()+'Z',hour=hour,valid=(run+dt.timedelta(hours=hour)).isoformat()+'Z',grid=GRID,sampling=dict(field_step_degrees=GRID['step'],profile_step_degrees=PGRID['step'],method='Nearest native model cell; no added model resolution'),fields=flds,fieldMeta=meta,profile=profile,inventory=inventory,source=source)

def save_frame(frame):
 key=f"{frame['model']}-{frame['run'][:13].replace('-','').replace('T','')}-{frame['hour']:03}-{frame['member']}"
 out=ROOT/'dist'/'data'/f'{key}.json.gz';out.parent.mkdir(exist_ok=True)
 with gzip.open(out,'wt',encoding='utf-8')as f:json.dump(serial(frame),f,separators=(',',':'),allow_nan=False)
 return dict(model=frame['model'],member=frame['member'],run=frame['run'],hour=frame['hour'],valid=frame['valid'],file='data/'+out.name,available=sorted(k for k in frame['fields'] if k in PRODUCTS),keys=sorted(frame['fields']),profile=bool(frame['profile']))

def update_catalog(entries):
 path=ROOT/'dist'/'catalog.json';old=json.loads(path.read_text()) if path.exists() else dict(schema=2,models=MODELS,products=PRODUCTS,frames=[])
 for e in entries:old['frames']=[f for f in old['frames'] if not (f['model']==e['model'] and f['member']==e['member'] and f['valid']==e['valid'])]+[e]
 old.update(models=MODELS,products={**old.get('products',{}),**PRODUCTS},updated=dt.datetime.now(dt.timezone.utc).isoformat())
 temp=path.with_suffix('.tmp');temp.write_text(json.dumps(old,separators=(',',':')));temp.replace(path)
 return old

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--models',nargs='+',default=list(MODELS));p.add_argument('--run',default='2026-09-14T12');p.add_argument('--hours',nargs='+',type=int,default=[12]);p.add_argument('--eps-members',nargs='+',default=['1','2']);p.add_argument('--sref-members',nargs='+',default=['arw-ctl','nmb-ctl']);args=p.parse_args();run=dt.datetime.strptime(args.run,'%Y-%m-%dT%H')
 for model in args.models:
  if model not in MODELS:p.error('Unknown model '+model)
  members=args.eps_members if model=='eps' else args.sref_members if model=='sref' else ['det']
  for member in members:
   for hour in args.hours:
    try:
     r=run-dt.timedelta(hours=3) if model=='sref' else run;h=hour+3 if model=='sref' else hour
     frame=build_frame(model,r,h,member);entry=save_frame(frame);update_catalog([entry]);print(f"Saved {model} {member}: {len(entry['available'])} fields; profile={entry['profile']}",flush=True)
    except Exception as e:print(f'UNAVAILABLE {model} {member} F{hour}: {e}',flush=True)
