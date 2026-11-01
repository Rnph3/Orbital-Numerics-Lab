"""Compute reproducible orbital comparisons and figures."""
import argparse, csv, json
from math import pi
from pathlib import Path
from .physics import integrate,diagnostics,circular_initial,eccentric_initial


def run(out:Path, *, eccentricity:float=0.0)->dict:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out.mkdir(parents=True,exist_ok=True)
    initial=eccentric_initial(eccentricity) if eccentricity else circular_initial()
    steps=160
    dt=2*pi/steps
    rows=[]
    fig,ax=plt.subplots(figsize=(7,6))
    ax.scatter([0],[0],marker='*',s=110,label='Central body')
    for method in ('euler','verlet','rk4'):
        states=integrate(initial,dt=dt,steps=steps,method=method)
        d=diagnostics(states)
        rows.append({'method':method,'steps':steps,'dt':dt,**d})
        ax.plot([s.x for s in states],[s.y for s in states],label=method,linewidth=1.6)
    ax.set_aspect('equal');ax.set_title(f'1 orbital period: e={eccentricity} (dimensionless)')
    ax.set_xlabel('x');ax.set_ylabel('y');ax.legend();ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'orbit_comparison.png',dpi=160);plt.close(fig)
    with (out/'orbit_results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    # Fixed time interval and grid refinement across all methods
    conv=[]
    for method in ('euler','verlet','rk4'):
        for n in (40,80,160,320):
            series=integrate(initial,dt=2*pi/n,steps=n,method=method)
            conv.append({'method':method,'steps':n,'endpoint_position_error':diagnostics(series)['endpoint_position_error']})
    with (out/'convergence.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(conv[0]));w.writeheader();w.writerows(conv)
    fig,ax=plt.subplots(figsize=(7,4))
    for method in ('euler','verlet','rk4'):
        part=[r for r in conv if r['method']==method]
        ax.loglog([r['steps'] for r in part],[r['endpoint_position_error'] for r in part],'o-',label=method)
    ax.set_xlabel('Time steps per period (more = smaller dt)');ax.set_ylabel('Endpoint position error')
    ax.set_title('Grid-refinement comparison, 1 period');ax.grid(alpha=.25,which='both');ax.legend();fig.tight_layout()
    fig.savefig(out/'convergence.png',dpi=160);plt.close(fig)
    result={'model':'dimensionless two-body point-mass test particle',
            'eccentricity':eccentricity,'period':2*pi,'comparison':rows,
            'convergence':conv,'provenance':'analytical initial conditions, no observational/mission data'}
    (out/'results.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=Path('artifacts'))
    p.add_argument('--eccentricity',type=float,default=0.0)
    a=p.parse_args();result=run(a.out,eccentricity=a.eccentricity)
    for x in result['comparison']:
        print(f"{x['method']}: energy drift {x['max_relative_energy_drift']:.5g}; endpoint err {x['endpoint_position_error']:.5g}")
