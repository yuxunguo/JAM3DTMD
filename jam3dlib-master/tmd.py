import os,sys
import numpy as np
from tools.config  import load_config, conf
from tools.tools   import load
from fitlib.resman import RESMAN
from obslib.sidis import upol0 as upol
from obslib.sidis import sivers0 as sivers
from obslib.sidis import collins0 as collins
from obslib.sidis import aUTsPs0 as sinphiS
from obslib.AN_pp import AN_theory0 as AN_theory

class TMD:

    def __init__(self,tag):

        load_config('data/input-%s.py'%tag)
        resman=RESMAN(nworkers=1,parallel=False,datasets=False)
        self.parman=resman.parman
        self.jar=np.load('data/jar-%s.npy'%tag,allow_pickle=True,encoding = 'latin1').item()
        self.parman.order=self.jar['order']
        self.par=self.jar['par']
        self.nrep=len(self.par)

    def eval(self,x,Q2,kT,had,dist,irep,icol=False,deriv=False):

        self.parman.set_new_params(self.par[irep],initial=True)

        if  dist=='pdf' and had=='p':
            return conf['pdf'].get_tmd(x,Q2,kT,'p','pdf',icol=icol,deriv=False)
        elif dist=='ff' and had=='pi':
            return conf['ffpi'].get_tmd(x,Q2,kT,'pi','ffpi',icol=icol,deriv=False)
        elif dist=='sivers' and had=='p':
            return conf['sivers'].get_tmd(x,Q2,kT,'p','sivers',icol=icol,deriv=False)
        elif dist=='transversity' and had=='p':
            return conf['transversity'].get_tmd(x,Q2,kT,'p','transversity',icol=icol,deriv=False)
        elif dist=='collinspi' and had=='pi':
            return conf['collinspi'].get_tmd(x,Q2,kT,'pi','collinspi',icol=icol,deriv=False)
        elif dist=='Htildepi' and had=='pi':
            return conf['Htildepi'].get_tmd(x,Q2,kT,'pi','Htildepi',icol=True,deriv=False)
        else:
            print('dist and had not available')
            sys.exit()

    def eval_stfunc(self,stfunc,x,z,Q2,pT,tar,had,irep,icol=False):

        self.parman.set_new_params(self.par[irep],initial=True)

        if stfunc=='FUU':
            return upol.get_FUU(x,z,Q2,pT,tar,had)
        elif stfunc=='FUTSiv':
            return sivers.get_FUT(x,z,Q2,pT,tar,had)
        elif stfunc=='FUTCol':
            return collins.get_FUT(x,z,Q2,pT,tar,had)
        elif stfunc=='FUTsinphiS':
            return sinphiS.get_FX(x,z,Q2,None,tar,had)
        else:
            print(stfunc,'is not available')
            sys.exit()

    def eval_TMDasym(self,asym,x,z,Q2,pT,tar,had,irep,icol=False):

        self.parman.set_new_params(self.par[irep],initial=True)

        if asym=='Siv': # sin(phi_s-phi_h)
            return self.eval_stfunc('FUTSiv',x,z,Q2,pT,tar,had,irep,icol)/self.eval_stfunc('FUU',x,z,Q2,pT,tar,had,irep,icol)
        elif asym=='Col':  # sin(phi_s+phi_h) NO DEPOLARIZATION FACTOR ADDED
            return self.eval_stfunc('FUTCol',x,z,Q2,pT,tar,had,irep,icol)/self.eval_stfunc('FUU',x,z,Q2,pT,tar,had,irep,icol)
        elif asym=='sinphiS':
            return self.eval_stfunc('FUTsinphiS',x,z,Q2,None,tar,had,irep,icol)/self.eval_stfunc('FUU',x,z,Q2,None,tar,had,irep,icol)
        else:
            print(stfunc,'is not available')
            sys.exit()

    def eval_ANasym(self,asym,xF,pT,rs,tar,had,irep):

        self.parman.set_new_params(self.par[irep],initial=True)

        if asym=='ANpp':
            den = AN_theory.get_sig(xF, pT, rs, tar, had, mode='gauss', nx=10, nz=10)
            num = AN_theory.get_sigST(xF, pT, rs, tar, had,mode='gauss', nx=10, nz=10)
            return num/den


if __name__ == '__main__':


    x = 0.25
    z = 0.5
    Q2 = 2.4
    pT = 0.3
    tar = 'p'
    had = 'pi+'
    stfunc='FUTCol'

    tag='JAM3D_2022'
    tmd=TMD(tag)

    print(tmd.eval_stfunc(stfunc, x, z, Q2, pT, tar, had, 0, icol=False))
    print(tmd.eval_asymmetry('Siv', x, z, Q2, pT, tar, had, 0, icol=False))
    print(tmd.eval_asymmetry('sinphiS', x, z, Q2, None, tar, had, 0, icol=False))
