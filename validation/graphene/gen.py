import numpy as np
a=2.46; L=10.0
a1=np.array([a,0,0]); a2=np.array([a/2,a*np.sqrt(3)/2,0]); a3=np.array([0,0,L])
def gen(n,fn):
    A=np.array([n*a1,n*a2,a3])
    atoms=[]
    for i in range(n):
        for j in range(n):
            for b in [(1/3,1/3,0),(2/3,2/3,0)]:
                frac=np.array([ (i+b[0])/n,(j+b[1])/n,0.0])
                atoms.append(frac)
    with open(fn,'w') as f:
        f.write(f"  {len(atoms)} F\n  C\n")
        for k,fr in enumerate(atoms): f.write(f"  {k+1} 1 {fr[0]:.12f} {fr[1]:.12f} {fr[2]:.12f}\n")
        f.write("  0 0 0\n")
        for v in A: f.write(f"  {v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n")
gen(3,'sc.gen'); gen(1,'prim.gen')
