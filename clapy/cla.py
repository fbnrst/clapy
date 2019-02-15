from . import ccla as _cla


def run(mode=2,samples = 1000,nCells=100,mCells=100,times=[0.0,0.2],GF=1,G1=0.2,S=0.3,G2M=0.5,sSamples=0.01,sPopulation=0.01,seed = 42,**args):
    '''
    mode             ... 1 daughter : 1 (int)
                         2 daughter : 2 (int)
    samples          ... number of samples per timepoint (int)
    nCells           ... number of cells per samples (int)
    mCells           ... number of measured cells (only 2 daughter)
    times            ... list of timepoints (floats)
    GF               ... growth fraction (float)
    G1               ... length of G1-phase (float)
    S                ... length of S-phase (float)
    G2M              ... length of G2- and M-phase (float)
                     ... cell cycle Tc = G1 + S + G2M
    sPopulations     ... std for Sample level (float)
    sSamples         ... std for Cell level (float)
    seed             ... seed

    '''
    return _cla.run(seed,samples,nCells,mCells,times,GF,G1,S,G2M,sSamples,sPopulation,mode)
