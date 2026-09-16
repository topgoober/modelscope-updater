"""Vectorized pressure-profile diagnostics on the sampled display grid.

Thermodynamic approximations: Bolton saturation vapor pressure/LCL, pseudoadiabatic
liquid-water parcel, hydrostatic virtual-temperature buoyancy integration. This is
not native-model or observational analysis. Pressure-level sampling limits detail.
"""
import numpy as np
RD,CP,EPS,G=287.05,1004.,.62195691,9.80665
KAPPA=RD/CP

def es(t):return 6.112*np.exp(17.67*(t-273.15)/(t-29.65))
def mixing(p,td):
 e=np.minimum(es(td),p*.98);return EPS*e/(p-e)
def dewpoint(p,r):
 e=np.maximum(p*r/(EPS+r),1e-5);l=np.log(e/6.112);return 273.15+243.5*l/(17.67-l)
def thetae(p,t,td):
 r=mixing(p,td);tl=1/(1/(td-56)+np.log(t/td)/800)+56
 return t*(1000/(p-es(td)))**(.2854*(1-.28*r))*np.exp((3376/tl-2.54)*r*(1+.81*r))
def lcl(p,t,td):
 td=np.minimum(td,t);tl=1/(1/(td-56)+np.log(t/td)/800)+56
 return p*(tl/t)**(1/KAPPA),tl

def lapse_derivative(p,t):
 r=mixing(p,t);lv=2.5e6
 return (RD*t+lv*r)/(CP+lv*lv*r*EPS/(RD*t*t))/p

def parcel(p,ps,ts,tds):
 pl,tl=lcl(ps,ts,tds);tcur=tl.copy();pcur=pl.copy();out=[]
 for target in p:
  # RK4 in pressure, <= 5 hPa steps, only for saturated parcel columns.
  loops=int(np.nanmax(np.maximum(0,pcur-target))/5)+1 if np.any(np.isfinite(pcur)) else 0
  for _ in range(loops):
   step=-np.minimum(5,np.maximum(0,pcur-target));k1=lapse_derivative(pcur,tcur);k2=lapse_derivative(pcur+step/2,tcur+step*k1/2);k3=lapse_derivative(pcur+step/2,tcur+step*k2/2);k4=lapse_derivative(pcur+step,tcur+step*k3)
   tcur=tcur+step*(k1+2*k2+2*k3+k4)/6;pcur=pcur+step
  out.append(np.where(target>ps,np.nan,np.where(target>=pl,ts*(target/ps)**KAPPA,tcur)))
 return np.array(out),pl

def buoyancy(p,z,t,td,ps,ts,tds,parcel_t):
 r0=mixing(ps,tds);rp=np.minimum(r0[None,:],mixing(p[:,None],parcel_t));re=mixing(p[:,None],td)
 tvp=parcel_t*(1+rp/EPS)/(1+rp);tve=t*(1+re/EPS)/(1+re)
 b=tvp-tve
 # Insert actual parcel start, interpolating environmental temperature at launch.
 cape=np.zeros(ps.shape);cin=np.zeros(ps.shape);cape03=np.zeros(ps.shape);lfc=np.full(ps.shape,np.nan);el=np.full(ps.shape,np.nan)
 prevp=ps.copy();prevb=np.zeros(ps.shape);prevz=np.zeros(ps.shape);found=np.zeros(ps.shape,bool)
 for i,pp in enumerate(p):
  ok=(pp<ps)&np.isfinite(b[i])&np.isfinite(z[i])&(z[i]>=0)
  bb=b[i];positive=ok&(bb>0);new=positive&(~found);frac=np.clip(-prevb/np.where(bb!=prevb,bb-prevb,1),0,1)
  lfc=np.where(new,prevz+(z[i]-prevz)*frac,lfc)
  # Split crossing layers instead of dropping their signed contributions.
  depth=RD*np.log(np.maximum(prevp,pp)/pp)
  cross=(prevb*bb<0)
  pos=np.where(cross,.5*np.maximum(prevb,bb)*(np.where(prevb>0,frac,1-frac)),.5*np.maximum(prevb+bb,0))*depth
  neg=np.where(cross,.5*np.minimum(prevb,bb)*(np.where(prevb<0,frac,1-frac)),.5*np.minimum(prevb+bb,0))*depth
  cape+=np.where(ok,pos,0);cin+=np.where(ok&(~found),neg,0)
  below=np.clip((3000-prevz)/np.maximum(z[i]-prevz,1),0,1);cape03+=np.where(ok,pos*below,0)
  ended=ok&found&(prevb>0)&(bb<=0);el=np.where(ended,prevz+(z[i]-prevz)*frac,el)
  found|=positive;prevp=np.where(ok,pp,prevp);prevb=np.where(ok,bb,prevb);prevz=np.where(ok,z[i],prevz)
 valid=np.isfinite(ts)&np.isfinite(tds)&np.any(np.isfinite(t),axis=0)
 return {'cape':np.where(valid,cape,np.nan),'cin':np.where(valid&found,cin,np.where(valid,0,np.nan)), 'cape03':np.where(valid,cape03,np.nan),'lfc':lfc,'el':el}

def interp_height(z,values,target,surface):
 out=np.full(surface.shape,np.nan);prevz=np.zeros(surface.shape);prevv=surface.copy()
 target=np.broadcast_to(target,surface.shape)
 for zz,v in zip(z,values):
  ok=np.isfinite(zz)&np.isfinite(v)&(zz>prevz)
  take=ok&np.isnan(out)&(target>=prevz)&(target<=zz)
  f=(target-prevz)/np.maximum(zz-prevz,1e-6);out=np.where(take,prevv+(v-prevv)*f,out)
  prevz=np.where(ok,zz,prevz);prevv=np.where(ok,v,prevv)
 return out

def diagnostics(p,pr,sfc):
 t,td,z,u,v=[pr[k] for k in ['t','td','z','u','v']]
 ps,ts,tds,u0,v0=[sfc[k] for k in ['ps','tmp','dpt','u10','v10']]
 ts=ts+273.15;tds=tds+273.15;t=t+273.15;td=td+273.15
 pt,pl=parcel(p,ps,ts,tds);d=buoyancy(p,z,t,td,ps,ts,tds,pt)
 d['lcl']=interp_height(z,np.broadcast_to(p[:,None],z.shape),0,ps) # replaced below with log-pressure interpolation
 d['lcl']=np.full(ps.shape,np.nan);prevp=ps.copy();prevz=np.zeros(ps.shape)
 for i,pp in enumerate(p):
  valid=np.isfinite(z[i])&(pp<prevp);take=valid&np.isnan(d['lcl'])&(pl<=prevp)&(pl>=pp)
  f=np.log(prevp/pl)/np.maximum(np.log(prevp/pp),1e-9);d['lcl']=np.where(take,prevz+(z[i]-prevz)*f,d['lcl']);prevp=np.where(valid,pp,prevp);prevz=np.where(valid,z[i],prevz)
 d['thetae']=thetae(ps,ts,tds)
 # Pressure-weighted 100-hPa mixed layer, including interpolated top endpoint.
 theta=t*(1000/p[:,None])**KAPPA;r=mixing(p[:,None],td)
 prevp=ps.copy();prevth=ts*(1000/ps)**KAPPA;prevr=mixing(ps,tds);sumth=np.zeros(ps.shape);sumr=np.zeros(ps.shape);mass=np.zeros(ps.shape)
 for i,pp in enumerate(p):
  ok=(pp<prevp)&np.isfinite(theta[i])&np.isfinite(r[i]);end=np.maximum(pp,ps-100);dp=np.maximum(prevp-end,0);frac=(prevp-end)/np.maximum(prevp-pp,1e-9)
  sumth+=np.where(ok,(prevth+prevth+(theta[i]-prevth)*frac)*.5*dp,0);sumr+=np.where(ok,(prevr+prevr+(r[i]-prevr)*frac)*.5*dp,0);mass+=np.where(ok,dp,0)
  prevth=np.where(ok,theta[i],prevth);prevr=np.where(ok,r[i],prevr);prevp=np.where(ok,pp,prevp)
 mt=sumth/np.maximum(mass,1)*(ps/1000)**KAPPA;md=dewpoint(ps,sumr/np.maximum(mass,1));mt=np.where(mass>=99,mt,np.nan)
 mpt,_=parcel(p,ps,mt,md);ml=buoyancy(p,z,t,td,ps,mt,md,mpt);d['mlcape'],d['mlcin']=ml['cape'],ml['cin']
 # Most-unstable parcel in the lowest 300 hPa, with the surface included.
 eth=thetae(p[:,None],t,td);eth=np.where((p[:,None]<=ps)&(p[:,None]>=ps-300),eth,-np.inf)
 eth=np.where(np.isfinite(eth),eth,-np.inf);idx=np.argmax(eth,axis=0);cols=np.arange(ps.size);better=eth[idx,cols]>d['thetae']
 mp=np.where(better,p[idx],ps);mt=np.where(better,t[idx,cols],ts);md=np.where(better,td[idx,cols],tds)
 mupt,_=parcel(p,mp,mt,md);mu=buoyancy(p,z,t,td,mp,mt,md,mupt);d['mucape']=mu['cape']
 if 500 in p:d['li']=t[list(p).index(500)]-pt[list(p).index(500)]
 t3=interp_height(z,t,3000,ts);d['lr03']=np.where(d['cape']>500,(ts-t3)/3,np.nan)
 if 700 in p and 500 in p:
  a,b=list(p).index(700),list(p).index(500);d['lr75']=(t[a]-t[b])/(z[b]-z[a])*1000
 us={h:interp_height(z,u,h,u0) for h in [500,1000,3000,5500,6000]};vs={h:interp_height(z,v,h,v0) for h in us}
 for h in [1,6]:d[f'shear{h}']=np.hypot(us[h*1000]-u0,vs[h*1000]-v0)*1.943844
 # Height-weighted Bunkers 0–6 km mean wind and 0–0.5 / 5.5–6 km shear.
 levels=np.arange(0,6001,250);uu=np.array([interp_height(z,u,h,u0) for h in levels]);vv=np.array([interp_height(z,v,h,v0) for h in levels])
 meanu=np.trapezoid(uu,levels,axis=0)/6000;meanv=np.trapezoid(vv,levels,axis=0)/6000
 du=np.mean(uu[-3:],axis=0)-np.mean(uu[:3],axis=0);dv=np.mean(vv[-3:],axis=0)-np.mean(vv[:3],axis=0);mag=np.hypot(du,dv)
 rm_u=meanu+7.5*dv/np.where(mag>.01,mag,np.nan);rm_v=meanv-7.5*du/np.where(mag>.01,mag,np.nan)
 d['storm_u'],d['storm_v']=rm_u,rm_v
 for h in [1,3]:
  n=h*4+1;su,sv=uu[:n]-rm_u,vv[:n]-rm_v;srh=np.sum(su[1:]*sv[:-1]-su[:-1]*sv[1:],axis=0);d[f'srh{h}']=srh;d[f'ehi{h}']=d['cape']*srh/160000
 sh=d['shear6']/1.943844;d['stp']=d['cape']/1500*np.clip((2000-d['lcl'])/1000,0,1)*d['srh1']/150*np.where(sh<12.5,0,np.minimum(sh,30)/20)
 # Effective inflow layer: contiguous sampled launch levels with CAPE >=100 and CIN >=-250.
 eff_base=np.where((d['cape']>=100)&(d['cin']>=-250),0,np.nan)
 eff_top=eff_base.copy();closed=np.zeros(ps.shape,bool)
 for i,pp in enumerate(p):
  if pp<500:break
  launch=np.isfinite(t[i])&(pp<ps)&(z[i]>=0)&(z[i]<=4000)
  lp=np.where(launch,pp,np.nan);lt=np.where(launch,t[i],np.nan);ld=np.where(launch,td[i],np.nan)
  lifted,_=parcel(p,lp,lt,ld);test=buoyancy(p,z,t,td,lp,lt,ld,lifted)
  qualifies=launch&(test['cape']>=100)&(test['cin']>=-250)
  start=qualifies&np.isnan(eff_base)&(~closed);eff_base=np.where(start,z[i],eff_base)
  closed|=launch&(~qualifies)&np.isfinite(eff_base)
  eff_top=np.where(qualifies&(~closed),z[i],eff_top)
 end=eff_base+.5*(mu['el']-eff_base)
 ub=interp_height(z,u,eff_base,u0);vb=interp_height(z,v,eff_base,v0)
 ue=interp_height(z,u,end,u0);ve=interp_height(z,v,end,v0);ebwd=np.hypot(ue-ub,ve-vb)
 path=np.linspace(0,1,15)[:,None]*(eff_top-eff_base)+eff_base
 eu=np.array([interp_height(z,u,h,u0) for h in path])-rm_u;ev=np.array([interp_height(z,v,h,v0) for h in path])-rm_v
 esrh=np.sum(eu[1:]*ev[:-1]-eu[:-1]*ev[1:],axis=0)
 d['scp']=mu['cape']/1000*esrh/50*np.where(ebwd<10,0,np.minimum(ebwd,20)/20)
 d['scp']=np.where(np.isnan(eff_base)&np.isfinite(d['cape']),0,d['scp'])
 d['lclp']=pl
 return d,pt-273.15
