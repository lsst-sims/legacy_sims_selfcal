import pylab
from lsst.sims.catalogs.measures.photometry.Bandpass import Bandpass
import lsst.sims.atmosphere.transmission.atmoComp as ac


def build_atmos(atmocmp=None, X=1.0, t0=(3.9/100.0), t1=(0.02/100.0), t2=(-0.03/100.0), alpha=-1.7,
                mol=0.96, BP=782, O3=0.9, H2O=0.9, doPlot=True):
    if atmocmp==None:
        atmocmp = ac.AtmoComp()
    # examples
    # max atmo = build_atmos(atmocmp, X=X, t0=5.6/100.0, alpha=-1.8, O3=1.5, H2O=1.3)
    # min atmo = build_atmos(atmocmp, X=X, t0=0.2/100.0, alpha=-0.5, O3=0.6, H2O=0.5)
    # 30p atmo = build_atmos(atmocmp, X=2.5, t0=(0.8/100), alpha=-1.0, O3=0.9, H2O=0.8)
    # 30p atmo = build_atmos(atmocmp, X=2.5, t0=(2.4/100.0), alpha=-1.4, O3=1.17, H2O=1.04)
    # 10p/30 atmo = build_atmos(atmocmp, X=X, t0=(0.8/100.0), alpha=-1.0, O3=0.9, H2O=0.8)
    # 10p/30 atmo = build_atmos(atmocmp, X=X, t0=(1.3/100.0), alpha=-1.13, O3=0.99, H2O=1.04)
    atmocmp.setCoefficients(t0=t0, t1=t1, t2=t2, alpha=alpha, mol=mol, BP=BP, O3=O3, H2O=H2O)
    atmocmp.buildAtmos(secz=X, doPlot=doPlot)
    atmos_bp = Bandpass(wavelen=atmocmp.wavelen, sb=atmocmp.trans_total)    
    return atmocmp, atmos_bp


if __name__ == '__main__':
    
    atmocomp, atmos1 = build_atmos(X=1.2, t0=(3.9/100.0), t1=(0.02/100.0), t2=(-0.03/100.0), alpha=-1.7,
                                   mol=0.96, BP=782, O3=1, H2O=1)
    atmocomp, atmos2= build_atmos(atmocmp=atmocomp, X=1.2, t0=(3.9/100.0), t1=(0.02/100.0),
                                   t2=(-0.03/100.0), alpha=-1.7,
                                   mol=0.96, BP=782, O3=1, H2O=1.5)
    atmocomp, atmos3 = build_atmos(atmocmp=atmocomp, X=1.2, t0=(3.9/100.0), t1=(0.02/100.0), \
                                   t2=(-0.03/100.0), alpha=-1.7, \
                                   mol=0.96, BP=782, O3=1, H2O=0.5)
    atmocomp, atmos4 = build_atmos(atmocmp=atmocomp, X=1.55, t0=(3.9/100.0), t1=(0.02/100.0), \
                                   t2=(-0.03/100.0), alpha=-1.7, \
                                   mol=0.96, BP=782, O3=1, H2O=0.5)

    pylab.show()

