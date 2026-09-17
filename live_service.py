"""Persistent, read-only forecast feed and supervised model ingestion.

Run `python live_service.py serve` on one Linux instance with a persistent
MODELSCOPE_DATA_DIR volume. Never run the poller from a web request.
"""
from __future__ import annotations
import argparse, datetime as dt, gzip, hashlib, json, logging, os, pathlib
import signal, subprocess, sys, tempfile, threading, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
from catalog import MODELS, PRODUCTS
from sources import source_url

ROOT = pathlib.Path(__file__).resolve().parent
UTC = dt.timezone.utc
LOG = logging.getLogger('modelscope.live')

def now(): return dt.datetime.now(UTC)
def stamp(value): return value.astimezone(UTC).isoformat().replace('+00:00','Z')
def parse(value): return dt.datetime.fromisoformat(value.replace('Z','+00:00'))
def atomic(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as f:
        tmp=pathlib.Path(f.name)
        try:
            f.write(payload); f.flush(); os.fsync(f.fileno())
            tmp.replace(path)
        finally: tmp.unlink(missing_ok=True)
def json_write(path, obj): atomic(path, json.dumps(obj,separators=(',',':'),allow_nan=False).encode())

def configuration():
    config=json.loads((ROOT/'live-settings.json').read_text())
    extra=os.environ.get('MODELSCOPE_SETTINGS')
    if extra: config.update(json.loads(extra))
    if not 60 <= config['poll_seconds'] <= 3600: raise ValueError('poll_seconds must be 60–3600')
    if config['valid_step_hours'] not in (3,6,12): raise ValueError('Use 3, 6 or 12 hour valid-time steps')
    if not 6 <= config['horizon_hours'] <= 48: raise ValueError('Display horizon must be 6–48 hours')
    if any(m not in MODELS for m in config['models']): raise ValueError('Unknown model')
    if not 1 <= config['workers'] <= 2: raise ValueError('Use 1–2 ingestion workers')
    if any(not str(m).isdigit() or not 1 <= int(m) <= 50 for m in config['eps_members']): raise ValueError('Invalid EPS member')
    allowed={c+'-'+m for c in ('arw','nmb') for m in ['ctl']+[s+str(i) for s in ('n','p') for i in range(1,7)]}
    if any(m not in allowed for m in config['sref_members']): raise ValueError('Invalid SREF member')
    return config

def cycles(model, at):
    """Candidate cycles only; existence is always checked in upstream inventories."""
    at=at.replace(minute=0,second=0,microsecond=0)
    runs=[at-dt.timedelta(hours=h) for h in range(37)]
    hours={'sref':(3,9,15,21),'hiresw':(0,12)}.get(model,(0,6,12,18))
    if model in ('hrrr','rap','rrfs'):
        # Recent hourly cycles plus preceding extended cycles for longer leads.
        extended=(3,9,15,21) if model=='rap' else (0,6,12,18)
        return [r for i,r in enumerate(runs) if i<4 or r.hour in extended][:10]
    return [r for r in runs if r.hour in hours]

def leads(model, run, at, config):
    limit={'hrrr':48 if run.hour%6==0 else 18,
           'rap':51 if run.hour in (3,9,15,21) else 21,
           'hiresw':48,'nam':84,'sref':87,'rrfs':60}.get(model,144)
    step=config['valid_step_hours']
    start=at.replace(minute=0,second=0,microsecond=0)
    start-=dt.timedelta(hours=start.hour%step)
    for offset in range(0,config['horizon_hours']+1,step):
        valid=start+dt.timedelta(hours=offset)
        hour=int((valid-run).total_seconds()/3600)
        if 0 <= hour <= limit: yield hour

def inventory_exists(model, run, hour, member):
    url=source_url(model,run,hour,member)
    url=url.removesuffix('.grib2')+'.index' if model in ('ecmwf','eps') else url+'.idx'
    request=urllib.request.Request(url,headers={'User-Agent':'Modelscope/3.0'})
    with urllib.request.urlopen(request,timeout=20) as response:
        raw=response.read(8_000_001)
    if len(raw)>8_000_000: raise ValueError('Inventory exceeds size limit')
    if model in ('ecmwf','eps'):
        rows=[json.loads(line) for line in raw.splitlines() if line.strip()]
        return any('_offset' in r and (model!='eps' or str(r.get('number'))==str(member)) for r in rows)
    return any(len(line.split(b':'))>=6 and line.split(b':')[0].isdigit() for line in raw.splitlines())

class Store:
    def __init__(self, directory, config):
        self.root=pathlib.Path(directory); self.root.mkdir(parents=True,exist_ok=True)
        (self.root/'data').mkdir(exist_ok=True)
        self.config=config
        path=self.root/'catalog.json'
        self.catalog=json.loads(path.read_text()) if path.exists() else dict(schema=2,models=MODELS,products=PRODUCTS,frames=[])
        # Keep the site's expanded parameter definitions when initializing a feed.
        seed=ROOT/'dist'/'catalog.json'
        if seed.exists(): self.catalog['products']={**json.loads(seed.read_text())['products'],**PRODUCTS}
        self.catalog['live']=dict(poll_seconds=config['poll_seconds'],valid_step_hours=config['valid_step_hours'],horizon_hours=config['horizon_hours'],history=config['history'])
        self.state=dict(state='starting',started_at=stamp(now()),heartbeat=stamp(now()),models={})
        self.save()
    def save(self): json_write(self.root/'catalog.json',self.catalog)
    def status(self):
        self.state['heartbeat']=stamp(now()); self.state['updated']=self.catalog.get('updated')
        json_write(self.root/'status.json',self.state)
    def publish(self, frame):
        model=frame['model']; member=frame['member']; valid=frame['valid']
        grid=frame['grid']; size=grid['width']*grid['height']
        if model not in MODELS or not frame['fields'] or any(len(a)!=size for a in frame['fields'].values()): raise ValueError('Invalid decoded frame')
        if parse(valid)!=parse(frame['run'])+dt.timedelta(hours=frame['hour']): raise ValueError('Run / lead mismatch')
        payload=json.dumps(frame,separators=(',',':'),allow_nan=False).encode()
        digest=hashlib.sha256(payload).hexdigest()[:20]
        filename=f'{model}-{parse(frame["run"]):%Y%m%d%H}-{frame["hour"]:03}-{member}-{digest}.json.gz'
        atomic(self.root/'data'/filename,gzip.compress(payload,mtime=0))
        e=dict(model=model,member=member,run=frame['run'],hour=frame['hour'],valid=valid,file='data/'+filename,
               keys=sorted(frame['fields']),available=sorted(k for k in frame['fields'] if k in self.catalog['products']),profile=bool(frame.get('profile')))
        older=[f for f in self.catalog['frames'] if (f['model'],f['member'],f['valid'])==(model,member,valid)]
        if older and any(parse(f['run'])>parse(frame['run']) for f in older): return
        cutoff=now()-dt.timedelta(hours=self.config['retention_hours'])
        self.catalog['frames']=[f for f in self.catalog['frames'] if f not in older and parse(f['valid'])>=cutoff]+[e]
        # Bound referenced forecast storage to fit small persistent volumes.
        budget=int(os.environ.get('MODELSCOPE_FRAME_BUDGET_MB','320'))*1_000_000
        entries=self.catalog['frames']
        sizes={e['file']:(self.root/e['file']).stat().st_size for e in entries}
        total=sum(sizes.values())
        while total>budget and len(entries)>1:
            # Retain nearer valid times; older initialization loses ties.
            victim=max(entries,key=lambda e:(abs((parse(e['valid'])-now()).total_seconds()),-parse(e['run']).timestamp()))
            entries.remove(victim);total-=sizes[victim['file']]
        self.catalog['updated']=stamp(now()); self.save()
        self.cleanup()
        self.state['models'].setdefault(model,{})['last_success']=stamp(now())
        self.state['models'][model]['last_run']=frame['run']
        self.state['models'][model].pop('error',None)
        self.status()
    def cleanup(self):
        keep={f['file'] for f in self.catalog['frames']}
        # Grace period protects clients holding the previous manifest.
        cutoff=time.time()-min(self.config['retention_hours']*3600,1800)
        for path in (self.root/'data').glob('*.json.gz'):
            if 'data/'+path.name not in keep and path.stat().st_mtime<cutoff: path.unlink(missing_ok=True)
        files=list((self.root/'data').glob('*.json.gz'))
        total=sum(p.stat().st_size for p in files)
        ceiling=int(os.environ.get('MODELSCOPE_FRAME_BUDGET_MB','320'))*1_000_000+60_000_000
        for path in sorted(files,key=lambda p:p.stat().st_mtime):
            if total<=ceiling: break
            if 'data/'+path.name not in keep:
                total-=path.stat().st_size;path.unlink(missing_ok=True)
        cache=self.root/'grib-cache'
        if cache.exists():
            for path in cache.iterdir():
                if path.is_file() and path.stat().st_mtime<cutoff: path.unlink(missing_ok=True)

def choose_job(model, store, checked, clock=now, probe=inventory_exists):
    at=clock(); config=store.config
    members=config['eps_members'] if model=='eps' else config['sref_members'] if model=='sref' else ['det']
    existing={(f['member'],f['valid']):parse(f['run']) for f in store.catalog['frames'] if f['model']==model}
    for run in cycles(model,at):
        for member in members:
            for hour in leads(model,run,at,config):
                key=(model,stamp(run),hour,member); valid=stamp(run+dt.timedelta(hours=hour))
                if existing.get((member,valid),dt.datetime.min.replace(tzinfo=UTC))>=run: continue
                previous=checked.get(key)
                if previous and previous[0]>time.monotonic():
                    if not previous[1]: break
                    continue
                try:
                    ready=probe(model,run,hour,member)
                    checked[key]=(time.monotonic()+config['poll_seconds'],ready)
                    if ready: return key
                except Exception as exc:
                    checked[key]=(time.monotonic()+config['poll_seconds'],False)
                    LOG.info('%s inventory pending/unavailable: %s',model,type(exc).__name__)
                # Forecast files normally arrive in lead-time order. Retry a gap
                # next polling interval instead of hammering every later index.
                break

    return None

def process_job(job, store):
    model,run,hour,member=job
    with tempfile.TemporaryDirectory(dir=store.root,prefix='job-') as directory:
        output=pathlib.Path(directory)/'frame.json.gz'
        cmd=[sys.executable,str(ROOT/'live_service.py'),'frame','--model',model,'--run',run,'--hour',str(hour),'--member',member,'--output',str(output)]
        if store.config['history']: cmd.append('--history')
        env={**os.environ,'MODELSCOPE_CACHE_DIR':str(pathlib.Path('/tmp')/('modelscope-'+pathlib.Path(directory).name)),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
        # Isolated decoder releases GRIB/native memory after every frame.
        try:
            from regional_service import decoder_slot
            with decoder_slot(store.root):
                subprocess.run(cmd,check=True,timeout=store.config['frame_timeout_seconds'],env=env)
        finally:
            import shutil
            shutil.rmtree(env['MODELSCOPE_CACHE_DIR'],ignore_errors=True)
        with gzip.open(output,'rt') as f: return json.load(f)

def run_poller(config, directory):
    import fcntl
    directory=pathlib.Path(directory); directory.mkdir(parents=True,exist_ok=True)
    lock=(directory/'updater.lock').open('w')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    store=Store(directory,config); checked={}; stopped=threading.Event()
    signal.signal(signal.SIGTERM,lambda *_:stopped.set())
    def heartbeat():
        while not stopped.wait(30): store.status()
    threading.Thread(target=heartbeat,daemon=True).start()
    with ThreadPoolExecutor(max_workers=config['workers']) as pool:
        while not stopped.is_set():
            store.state['state']='checking'; store.status(); jobs=[]
            for model in config['models']:
                job=choose_job(model,store,checked)
                store.state['models'].setdefault(model,{})['last_checked']=stamp(now())
                if job: jobs.append(job)
                if stopped.is_set(): break
            store.state['state']='processing' if jobs else 'waiting'; store.state['queued']=len(jobs); store.status()
            futures=[(job,pool.submit(process_job,job,store)) for job in jobs]
            for job,future in futures:
                try: store.publish(future.result())
                except Exception as exc:
                    LOG.exception('Frame failed: %s',job)
                    store.state['models'].setdefault(job[0],{})['error']='Latest attempt failed; retry scheduled ('+type(exc).__name__+')'
                checked[job]=(time.monotonic()+config['poll_seconds'],True)
                store.state['queued']-=1; store.status()
            store.cleanup()
            checked={key:until for key,until in checked.items() if until[0]>time.monotonic()}
            if not jobs: stopped.wait(config['poll_seconds'])
    store.state['state']='stopped'; store.status()

def create_feed_app(directory=None):
    from flask import Flask,abort,jsonify,send_from_directory,request
    app=Flask(__name__,static_folder=None)
    root=pathlib.Path(directory or os.environ.get('MODELSCOPE_DATA_DIR','/data/modelscope'))
    @app.after_request
    def headers(response):
        origin=os.environ.get('MODELSCOPE_SITE_ORIGIN','https://modelscope-weather.wilcob9139.chatgpt.site')
        response.headers['Access-Control-Allow-Origin']=origin
        response.headers['Vary']='Origin'
        response.headers['X-Content-Type-Options']='nosniff'
        response.headers['Cache-Control']='public, max-age=86400, immutable' if request.path.startswith('/data/') else 'no-store'
        return response
    @app.get('/healthz')
    def health():
        try:
            state=json.loads((root/'status.json').read_text())
            healthy=(now()-parse(state['heartbeat'])).total_seconds()<120 and state['state']!='stopped'
            return jsonify(ok=healthy),200 if healthy else 503
        except (OSError,ValueError,KeyError): return jsonify(ok=False),503
    @app.get('/catalog.json')
    @app.get('/status.json')
    def manifest(): return send_from_directory(root,request.path.lstrip('/'))
    @app.get('/data/<filename>')
    def frame_file(filename):
        if not filename.endswith('.json.gz'): abort(404)
        return send_from_directory(root/'data',filename,mimetype='application/gzip',conditional=True)
    from regional_service import install_regional_routes
    install_regional_routes(app,root)
    return app

def serve(config, directory):
    # Web process and updater have separate lifetimes; failure restarts the service.
    poller=subprocess.Popen([sys.executable,str(ROOT/'live_service.py'),'poll'],env=os.environ)
    web=subprocess.Popen([sys.executable,'-m','gunicorn','--bind','0.0.0.0:'+os.environ.get('PORT','8000'),'--workers','1','--threads','4','live_service:create_feed_app()'],cwd=ROOT)
    def stop(*_):
        for child in (poller,web):
            if child.poll() is None: child.terminate()
    signal.signal(signal.SIGTERM,stop); signal.signal(signal.SIGINT,stop)
    try:
        while poller.poll() is None and web.poll() is None: time.sleep(1)
    finally:
        stop()
        for child in (poller,web):
            try: child.wait(timeout=15)
            except subprocess.TimeoutExpired: child.kill(); child.wait()
    raise SystemExit(1)

def main():
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('command',choices=['serve','poll','frame'])
    p.add_argument('--model'); p.add_argument('--run'); p.add_argument('--hour',type=int); p.add_argument('--member',default='det'); p.add_argument('--output'); p.add_argument('--history',action='store_true')
    args=p.parse_args()
    if args.command=='frame':
        from ingest import build_frame,serial
        frame=build_frame(args.model,parse(args.run).replace(tzinfo=None),args.hour,args.member)
        if args.history:
            from history import augment
            frame=augment(frame)
        with gzip.open(args.output,'wt') as f: json.dump(serial(frame),f,separators=(',',':'),allow_nan=False)
        return
    config=configuration(); directory=os.environ.get('MODELSCOPE_DATA_DIR','/data/modelscope')
    if args.command=='poll': run_poller(config,directory)
    else: serve(config,directory)
if __name__=='__main__': main()
