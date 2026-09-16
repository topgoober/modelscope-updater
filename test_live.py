import copy, datetime as dt, json, pathlib, tempfile, unittest
from unittest.mock import patch
from live_service import Store,configuration,now,stamp,choose_job,create_feed_app,cycles,leads

class LiveTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(); self.root=pathlib.Path(self.tmp.name)
  self.config=configuration(); self.store=Store(self.root,self.config)
  self.run=now().replace(minute=0,second=0,microsecond=0)-dt.timedelta(hours=6)
 def tearDown(self): self.tmp.cleanup()
 def frame(self,run=None,hour=12):
  run=run or self.run
  return dict(model='hrrr',member='det',run=stamp(run),hour=hour,valid=stamp(run+dt.timedelta(hours=hour)),grid=dict(width=2,height=2),fields=dict(cape=[0,100,None,500]),profile=None)
 def test_atomic_publish_and_newer_run(self):
  old=self.frame(); self.store.publish(old); first=json.loads((self.root/'catalog.json').read_text())
  self.assertTrue((self.root/first['frames'][0]['file']).exists())
  new=self.frame(self.run+dt.timedelta(hours=1),11);self.store.publish(new);self.store.publish(old)
  self.assertEqual(len(self.store.catalog['frames']),1)
  self.assertEqual(self.store.catalog['frames'][0]['run'],new['run'])
  reopened=Store(self.root,self.config);self.assertEqual(reopened.catalog['frames'],self.store.catalog['frames'])
 def test_bad_frame_keeps_catalog(self):
  self.store.publish(self.frame());before=(self.root/'catalog.json').read_bytes()
  bad=self.frame();bad['fields']['cape']=[1]
  with self.assertRaises(ValueError):self.store.publish(bad)
  self.assertEqual(before,(self.root/'catalog.json').read_bytes())
 def test_cycles_and_utc_alignment(self):
  at=dt.datetime(2026,9,16,1,tzinfo=dt.timezone.utc)
  self.assertEqual(cycles('sref',at)[0].hour,21)
  self.assertEqual(cycles('hiresw',at)[0].hour,0)
  for model in self.config['models']:
   for run in cycles(model,at):
    for hour in leads(model,run,at,self.config):self.assertEqual((run+dt.timedelta(hours=hour)).hour%3,0)
 def test_discovery_fallback_and_no_duplicate(self):
  at=now().replace(minute=0,second=0,microsecond=0)
  candidates=cycles('hrrr',at);previous=candidates[1]
  def probe(model,run,hour,member):return run==previous
  job=choose_job('hrrr',self.store,{},lambda:at,probe)
  self.assertEqual(job[1],stamp(previous))
  frame=self.frame(previous,job[2]);self.store.publish(frame)
  another=choose_job('hrrr',self.store,{},lambda:at,probe)
  self.assertNotEqual(job,another)
 def test_readonly_feed_health_and_cors(self):
  self.store.publish(self.frame());self.store.status();client=create_feed_app(self.root).test_client()
  response=client.get('/catalog.json');self.assertEqual(response.status_code,200)
  self.assertEqual(response.headers['Cache-Control'],'no-store')
  self.assertEqual(response.headers['Access-Control-Allow-Origin'],'https://modelscope-weather.wilcob9139.chatgpt.site')
  path='/'+self.store.catalog['frames'][0]['file'];self.assertEqual(client.get(path).status_code,200)
  self.assertEqual(client.post('/catalog.json').status_code,405)
  self.assertEqual(client.get('/data/../../live-settings.json').status_code,404)
  self.assertEqual(client.get('/healthz').status_code,200)
  self.store.state['state']='stopped';self.store.status();self.assertEqual(client.get('/healthz').status_code,503)
if __name__=='__main__':unittest.main()
