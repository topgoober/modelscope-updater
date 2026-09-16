"""Explicit product definitions; missing inputs never become synthetic forecasts."""
MODELS = {
'hrrr':dict(name='HRRR',resolution='3 km',color='#bbeb76',family='deterministic'),
'gfs':dict(name='GFS',resolution='0.25°',color='#65c7f2',family='deterministic'),
'nam':dict(name='NAM',resolution='12 km',color='#db9be9',family='deterministic'),
'ecmwf':dict(name='ECMWF IFS',resolution='0.25° open data',color='#ffbb70',family='deterministic'),
'rap':dict(name='RAP',resolution='13 km',color='#ff7e96',family='deterministic'),
'hiresw':dict(name='HRW WRF-ARW',resolution='5 km product',color='#55dfce',family='deterministic'),
'rrfs':dict(name='RRFS',resolution='3 km · parallel',color='#a8a2ff',family='deterministic'),
'sref':dict(name='SREF',resolution='40 km product',color='#edd96d',family='ensemble'),
'eps':dict(name='ECMWF EPS',resolution='0.25° open data',color='#76a8ff',family='ensemble'),
}
PRODUCTS={}
def add(key,name,group,unit,domain,levels=None,base=None,wind=None,**kw):
 PRODUCTS[key]=dict(name=name,group=group,unit=unit,domain=domain,levels=levels or [round(domain[0]+(domain[1]-domain[0])*v,2) for v in [.25,.5,.75]],base=base or key,wind=wind,**kw)
for p,d in [(200,[11000,13000]),(300,[8500,10000]),(500,[5100,6000]),(700,[2700,3300]),(850,[1200,1700])]:
 add(f'z{p}',f'{p} mb Height, Wind','Upper air','m',d,wind=f'{p}')
for p in [700,850]:add(f't{p}',f'{p} mb Temperature, Height, Wind','Upper air','°C',[-30,30],wind=str(p),height=f'z{p}')
add('rh','2 m AGL Relative Humidity','Surface','%',[0,100])
add('tmp','2 m AGL Temperature','Surface','°C',[-20,45])
add('tmp_wind','2 m AGL Temperature, Wind Barbs','Surface','°C',[-20,45],base='tmp',wind='10')
add('apparent','2 m AGL Wind Chill / Heat Index','Surface','°C',[-30,50],note='NWS heat index / wind chill where applicable; air temperature otherwise.')
add('dpt','2 m AGL Dew Point','Surface','°C',[-20,30])
add('dpt_wind','2 m AGL Dew Point, Wind Barbs','Surface','°C',[-20,30],base='dpt',wind='10')
add('thetae','2 m AGL Theta-e, Wind Barbs','Surface','K',[270,380],wind='10')
add('mslp','MSLP, 10 m AGL Wind','Surface','hPa',[980,1040],wind='10')
add('gust','10 m AGL Wind Gusts','Surface','kt',[0,70])
add('smoke_sfc','Near-Surface Smoke Density','Fire weather','µg/m³',[0,100])
add('smoke_col','Vertically Integrated Smoke','Fire weather','mg/m²',[0,300])
add('ptype','Precipitation Type','Precipitation type','category',[0,4],categorical=['None','Rain','Snow','Ice pellets','Freezing rain'])
add('ref_ptype','1 km AGL Reflectivity with Precip Type','Precipitation type','dBZ',[0,70],base='ref1',overlay='ptype')
for n in [1,6,12,24]:add(f'qpf{n}',f'{n}-h QPF','Quantitative precipitation','mm',[0,75],note=f'Exact preceding {n}-hour accumulation. Missing history is unavailable, never prorated.')
add('qpf_total','Total QPF','Quantitative precipitation','mm',[0,150],note='Accumulation since the displayed model initialization.')
add('cloud','Cloud Cover','Moisture & satellite','%',[0,100])
add('cloud_levels','Cloud Cover, Levels','Moisture & satellite','%',[0,100],base='cloud',overlay='cloudlevels')
for key,label in [('cloud_low','Low cloud cover'),('cloud_mid','Middle cloud cover'),('cloud_high','High cloud cover')]:add(key,label,'Moisture & satellite','%',[0,100])
add('lightning','Lightning Flash Density','Moisture & satellite','flashes/km²/5 min',[0,5],note='Only records with matching density units and time window are accepted.')
add('pwat','Precipitable Water','Moisture & satellite','mm',[0,65])
add('ir','Simulated IR Satellite','Moisture & satellite','°C',[-80,35],note='GOES-11 channel 4 (10.7 µm) simulated brightness temperature from the selected feed; see native record metadata.')
add('visibility','Visibility','Moisture & satellite','km',[0,30])
add('ref1','1 km AGL Reflectivity','Radar products','dBZ',[0,70])
add('refc','Composite Reflectivity','Radar products','dBZ',[0,70])
for p in [700,850]:add(f'adv{p}',f'{p} mb Temperature Advection','Temperature advection','°C/h',[-3,3],diverging=True,note='Horizontal advection on the common display grid; terrain-masked.')
for key,name,dom in [('cape03','0–3 km AGL CAPE',[0,300]),('mlcape','Mixed-Layer CAPE',[0,4000]),('mlcin','Mixed-Layer CIN',[-400,0]),('mucape','Most Unstable CAPE',[0,5000]),('cape','Surface-Based CAPE',[0,4000]),('cin','Surface-Based CIN',[-400,0]),('lcl','Surface-Based LCL Height',[0,3000]),('li','Surface-Based Lifted Index',[-12,12])]:
 add(key,name,'Instability','m' if key=='lcl' else '°C' if key=='li' else 'J/kg',dom,note='100-hPa mixed-layer parcel.' if key.startswith('ml') else 'Most-unstable parcel in lowest 300 hPa.' if key=='mucape' else 'Surface parcel.' )
add('lr03','Lapse Rate: 0–3 km AGL, CAPE > 500','Instability','°C/km',[3,10])
add('lr75','Lapse Rate: 700–500 mb','Instability','°C/km',[3,10])
for h in [1,6]:add(f'shear{h}',f'Bulk Shear: 0–{h} km AGL','Wind shear','kt',[0,80])
for h in [1,3]:add(f'srh{h}',f'Storm Relative Helicity: 0–{h} km AGL','Wind shear','m²/s²',[-100,500],note='Bunkers right-moving storm motion; pressure-level profile interpolation.')
for h in [1,3]:add(f'ehi{h}',f'Energy Helicity Index: 0–{h} km AGL','Composite parameters','index',[0,8],note='SBCAPE × signed SRH / 160,000.')
add('scp','Supercell Composite','Composite parameters','index',[0,15],note='Requires effective-layer SRH and effective bulk wind difference. Not replaced with fixed-layer shear.')
add('stp','Significant Tornado Parameter (fixed-layer)','Composite parameters','index',[0,8],note='Fixed-layer STP: SBCAPE, surface LCL, 0–1 km SRH and 0–6 km shear; not effective-layer STP.')
add('ref_uh','Reflectivity, UH > 75','Explicit convection','dBZ',[0,70],base='refc',overlay='uh')
for k,n in [('uh03_run','Updraft Helicity: 0–3 km AGL (run max)'),('uh25_1','Updraft Helicity: 2–5 km AGL (1 h max)'),('uh25_3','Updraft Helicity: 2–5 km AGL (3 h max)'),('uh25_run','Updraft Helicity: 2–5 km AGL (run max)')]:add(k,n,'Explicit convection','m²/s²',[0,300],note='Maxima require complete matching hourly history.')
add('cape_cross','SBCAPE, 500–850–Sfc Crossover','Combination plots','J/kg',[0,4000],base='cape',overlay='crossover')
add('cape_hodo','SBCAPE, Hodographs','Combination plots','J/kg',[0,4000],base='cape',overlay='hodographs')
