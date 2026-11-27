import numpy as np, sys
fermi=float(sys.argv[1])
def read_unfold(fn):
    with open(fn) as f:
        nq,ne,emin,emax,ef=f.readline().split(); nq=int(nq); ne=int(ne); emin=float(emin); emax=float(emax)
        W=np.zeros((nq,ne)); Q=np.zeros((nq,3))
        for i in range(nq):
            Q[i]=[float(x) for x in f.readline().split()[:3]]
            for j in range(ne): W[i,j]=float(f.readline())
    E=np.linspace(emin,emax,ne)
    return Q,E,W
def read_band(fn):
    ev=[];cur=None
    for l in open(fn):
        t=l.split()
        if not t: continue
        if t[0]=='KPT': cur=[]; ev.append(cur)
        else: cur.append(float(t[1]))
    return np.array(ev)
Q,E,W=read_unfold('run2/unfold.refoldedBands')
prim=read_band('prim/bands/band.out')  # (46, nbands)
Eabs=E+fermi
print("nq",len(Q)," bands in primitive:",prim.shape)
print(" q#  totalWeight  frac_within_0.15eV_of_prim_bands   ")
for iq in range(0,len(Q),5):
    tot=W[iq].sum()
    near=np.zeros(len(E),bool)
    for e in prim[iq]: near |= abs(Eabs-e)<0.15
    print(f"{iq:3d}  {tot:8.4f}   {W[iq][near].sum()/tot:8.4f}")
# global
tot=W.sum(); 
nearfrac=0
for iq in range(len(Q)):
    near=np.zeros(len(E),bool)
    for e in prim[iq]: near |= abs(Eabs-e)<0.15
    nearfrac+=W[iq][near].sum()
print("global total weight per q:",tot/len(Q)," fraction near primitive bands:",nearfrac/tot)
# energy-weighted centroid error: for each prim band find weighted mean energy in +-0.4 eV window
err=[]
for iq in range(len(Q)):
    for e in prim[iq]:
        m=abs(Eabs-e)<0.4
        if W[iq][m].sum()>0.5: err.append(( (W[iq][m]*Eabs[m]).sum()/W[iq][m].sum())-e)
print("centroid deviation (eV): mean %.4f  max abs %.4f  n=%d"%(np.mean(err),np.max(np.abs(err)),len(err)))
