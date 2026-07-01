#!/usr/bin/env python
import os,sys
import numpy as np
from tmd import TMD


def main1():

    tag='JAM3D_2020'
    x=0.1
    Q2=10.0
    kT=0.1
    tmd=TMD(tag)
    print(tmd.eval(x,Q2,kT,'p','pdf',0))


if __name__=="__main__":

    main1()
