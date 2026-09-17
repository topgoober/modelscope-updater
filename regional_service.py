"""Bounded, on-demand regional forecast jobs; no decoding in request threads."""
import concurrent.futures, contextlib, fcntl, hashlib, json, logging, os
import pathlib, re, shutil, subprocess, sys, tempfile, threading, time
from region_grid import regional_grid

@contextlib.contextmanager
def decoder_slot(root):
    root=pathlib.Path(root);root.mkdir(parents=True,exist_ok=True)
    with (root/'decoder.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        try:yield
        finally:fcntl.flock(lock,fcntl.LOCK_UN)

def install_regional_routes(app, root):
    from flask import request, jsonify, send_from_directory, abort
    root=pathlib.Path(root);cache=root/'regional';cache.mkdir(parents=True,exist_ok=True)
    pool=concurrent.futures.ThreadPoolExecutor(max_workers=1)
    jobs={};mutex=threading.Lock();code=pathlib.Path(__file__).resolve().parent
    budget=80_000_000

    def cleanup(keep=None):
        paths=sorted(cache.glob('*.json.gz'),key=lambda p:p.stat().st_mtime)
        total=sum(p.stat().st_size for p in paths)
        for path in paths:
            size=path.stat().st_size
            if path.name!=keep and (total>budget or time.time()-path.stat().st_mtime>7200):
                path.unlink(missing_ok=True);total-=size

    def work(key,entry,lat,lon):
        try:
            with decoder_slot(root), tempfile.TemporaryDirectory(prefix='modelscope-region-') as directory:
                with mutex:jobs[key]['state']='processing'
                output=pathlib.Path(directory)/'frame.json.gz'
                env={**os.environ,'MODELSCOPE_REGION_CENTER':json.dumps([lat,lon]),
                     'MODELSCOPE_CACHE_DIR':str(pathlib.Path(directory)/'grib'),
                     'MODELSCOPE_EPHEMERAL_HISTORY':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
                cmd=[sys.executable,str(code/'live_service.py'),'frame','--model',entry['model'],
                     '--run',entry['run'],'--hour',str(entry['hour']),'--member',entry['member'],
                     '--output',str(output),'--history']
                subprocess.run(cmd,env=env,check=True,timeout=1800)
                if output.stat().st_size>budget:raise ValueError('Regional frame exceeds cache budget')
                # Publish complete gzip files atomically on the persistent filesystem.
                target=cache/(key+'.json.gz');temporary=target.with_suffix('.tmp')
                shutil.copyfile(output,temporary);temporary.replace(target);cleanup(target.name)
            with mutex:jobs[key]={'state':'ready','time':time.time()}
        except Exception:
            logging.exception('Regional decoder failed: %s',key)
            with mutex:jobs[key]={'state':'failed','time':time.time()}

    @app.get('/regional/capabilities.json')
    def capabilities():
        return jsonify(schema=1,available=True,step=.03,width=161,height=121,products='all_available',profiles=True)

    @app.get('/regional/request')
    def request_region():
        try:
            lat=float(request.args['lat']);lon=float(request.args['lon']);grid=regional_grid(lat,lon)
            model=request.args['model'];member=request.args.get('member','det');run=request.args['run'];hour=int(request.args['hour'])
            catalog=json.loads((root/'catalog.json').read_text())
            entry=next((e for e in catalog['frames'] if (e['model'],e['member'],e['run'],e['hour'])==(model,member,run,hour)),None)
            if entry is None:return jsonify(state='unavailable',message='This model/run/time is not in the live catalog.'),404
        except (KeyError,ValueError,TypeError,OSError):return jsonify(state='unavailable',message='Invalid region or forecast request.'),400
        key=hashlib.sha256(json.dumps([model,member,run,hour,grid],sort_keys=True).encode()).hexdigest()[:32]
        with mutex:
            path=cache/(key+'.json.gz')
            if path.exists():
                path.touch();return jsonify(state='ready',file='regional/'+path.name,grid=grid)
            previous=jobs.get(key)
            if previous and previous['state'] in ('queued','processing'):return jsonify(state=previous['state'],retry_seconds=10),202
            if previous and previous['state']=='failed' and time.time()-previous['time']<300:
                return jsonify(state='failed',message='Regional download or decoding failed; retry in five minutes.'),503
            for old,item in list(jobs.items()):
                if item['state'] not in ('queued','processing') and time.time()-item['time']>300:jobs.pop(old,None)
            if sum(j['state'] in ('queued','processing') for j in jobs.values())>=3:
                return jsonify(state='busy',message='Regional queue is full; retry shortly.',retry_seconds=20),429
            jobs[key]={'state':'queued','time':time.time()};pool.submit(work,key,entry,lat,lon)
        return jsonify(state='queued',retry_seconds=10),202

    @app.get('/regional/<filename>')
    def regional_file(filename):
        if not re.fullmatch(r'[0-9a-f]{32}\.json\.gz',filename):abort(404)
        return send_from_directory(cache,filename,mimetype='application/gzip',conditional=True)
