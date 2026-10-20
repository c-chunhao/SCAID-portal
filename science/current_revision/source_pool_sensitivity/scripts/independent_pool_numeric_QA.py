"""Independent least-squares/group-score recheck of additional pool sensitivity.
Uses saved individual model rows locally; contains and prints no individual IDs.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np,pandas as pd,patsy
from scipy import linalg,stats

def manual_cr1(X,y,groups):
    N,K=X.shape;levels=np.unique(groups);G=len(levels)
    beta,_,rank,_=linalg.lstsq(X,y,lapack_driver='gelsy')
    if rank!=K:raise ValueError('Full-rank design required')
    residual=y-X@beta
    scores=np.vstack([X[groups==g].T@residual[groups==g] for g in levels])
    # QR-based inverse information, independent of statsmodels covariance helpers.
    _,R=linalg.qr(X,mode='economic');Ri=linalg.solve_triangular(R,np.eye(K));bread=Ri@Ri.T
    cov=bread@scores.T@scores@bread*(G/(G-1))*((N-1)/(N-K))
    return beta,cov,G

def self_test():
    rng=np.random.default_rng(512);g=np.repeat(np.arange(20),10);x=rng.normal(size=200);c=np.tile([0,1],100)
    X=np.column_stack([np.ones(200),x,c]);y=.3*c+.1*x+np.repeat(rng.normal(0,.3,20),10)+rng.normal(0,.5,200)
    beta,cov,G=manual_cr1(X,y,g)
    assert G==20 and np.linalg.eigvalsh(cov).min()>0 and np.all(np.isfinite(beta))
    return {'passed':True,'input':'200 fully synthetic observations in 20 groups','coefficient':float(beta[2])}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--directory',type=Path,default=Path(__file__).resolve().parents[1]);a.add_argument('--self-test',action='store_true');args=a.parse_args()
    if args.self_test:print(json.dumps(self_test()));raise SystemExit
    W=args.directory;models=pd.read_csv(W/'tables/source_pool_sensitivity_models.tsv',sep='\t');rows=pd.read_csv(W/'tables/PRIVATE_model_input_rows.tsv',sep='\t')
    checks=[];all_values={}
    for _,r in models[models.status.eq('estimated')].iterrows():
        d=rows[rows.model_key.eq(r.model_key)].copy();d['cohort']=d.cohort.astype(str);d['condition']=pd.Categorical(d.condition,categories=['HD','SLE'])
        y,X=patsy.dmatrices(r.formula,d,return_type='dataframe');i=list(X.columns).index('C(condition)[T.SLE]')
        b,V,G=manual_cr1(X.to_numpy(),y.to_numpy().ravel(),d.batch_cov.to_numpy());se=np.sqrt(V[i,i]);critical=stats.t.ppf(.975,G-1)
        actual=np.array([b[i],se,b[i]-critical*se,b[i]+critical*se,2*stats.t.sf(abs(b[i]/se),G-1)])
        expected=r[['coefficient','CR1_SE','CI_low','CI_high','p_raw']].to_numpy(float)
        error=float(abs(actual-expected).max());assert np.allclose(actual,expected,atol=1e-9,rtol=1e-7),(r.model_key,error)
        checks.append({'model_key':r.model_key,'variant':r.variant,'max_abs_error':error,'N':len(d),'G':G})
    # Recompute BH by sorted-rank formula, including unavailable tests as P=1.
    bh_checks=[]
    for (variant,family),d in models.groupby(['variant','family']):
        p=d.p_raw.fillna(1).to_numpy();order=np.argsort(p);ordered=p[order];adjusted=np.minimum.accumulate((ordered*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1);q=np.empty(len(p));q[order]=adjusted
        available=d.status.eq('estimated').to_numpy();err=float(abs(q[available]-d.BH_q.to_numpy()[available]).max());assert err<1e-12
        bh_checks.append({'variant':variant,'family':family,'size':len(p),'max_abs_error':err})
    result={'passed':True,'method':'Independent scipy pivoted-QR least squares, QR inverse-information, grouped-score CR1 covariance, scipy t(G-1), manual BH sorted-rank formula.','n_models':len(checks),'max_abs_error':max(c['max_abs_error'] for c in checks),'BH_families':bh_checks,'checks':checks,'synthetic_self_test':self_test(),'input_models_sha256':hashlib.sha256((W/'tables/source_pool_sensitivity_models.tsv').read_bytes()).hexdigest()}
    (W/'provenance/independent_numeric_QA.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['checks','BH_families']}))
