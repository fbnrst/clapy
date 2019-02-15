import numpy as np
import scipy as sp
import scipy.stats
import pandas as pd
idx = pd.IndexSlice

def calc_statistic(dataframe,parameter):
    tmp = dataframe.copy()
    try:
        tmp.rename({'GF_error' : 'GFf_error'}, axis="columns",inplace=True)
    except:
        pass
    for i in ['Tc','r','ss','sc','GFf']:
        try:
            tmp[i+'_2std'] = tmp[i+'_error'].mul(2)
            tmp[i+'_mse'] = tmp[i].sub(1).pow(2)
            tmp[i+'_d'] = tmp[i].sub(parameter[i])
            tmp[i+'_p'] = tmp[i+'_d'].abs() < tmp[i+'_2std']
            tmp[i+'_p2'] = np.logical_and( tmp[i+'_hpd46_l'] < parameter[i] ,parameter[i] < tmp[i+'_hpd46_u']    )
        except:
            print(i+' not found')
    sdata = tmp.mean(level=['model','sPopulation','sSample','M','GF'])
    for i in ['Tc','r','ss','sc','GFf']:
        try:
            sdata[i+'_mse'] = sdata[i+'_mse'].pow(0.5)
        except:
            print(i+' not found')
    return tmp,sdata

def ecdf(x):
    xs = np.sort(x)
    ys = np.arange(1, len(xs)+1)/float(len(xs))
    return xs, ys


class asym_lh:#cdef is here to reduce name lookup for __call__

    def __init__(self,data,times,ncell):
        ''' data = labeling fraction
            times = time for labeling fraction
            ncell = number of cells
        '''
        self.data = np.round(data*ncell)
        self.datalen = np.size(data)
        self.times = times
        if np.size(ncell) !=  self.datalen:
            self.ncell = np.ones_like(data,dtype=np.int32)*ncell
        else:
            self.ncell = ncell

    def compute(self, Tc,r,GF,sigma_cell,sigma_sample):
        ''' compute log liklyhood for parameters given
        '''
        pmf = self.pmf_f(Tc,r,GF,sigma_cell,sigma_sample) 
        pmf[np.abs(pmf) < 1e-300] = 1e-300 #fix nan in log
        return np.sum(-np.log( pmf) )
        #return np.sum(-np.log( self.pmf_f(Tc,r,GF,sigma_cell,sigma_sample) ) )


    def log_params(self, mu, sigma):
        """ A transformation of paramteres such that mu and sigma are the
            mean and variance of the log-normal distribution and not of
            the underlying normal distribution.
        """
        s2 = np.log(1.0 + sigma**2/mu**2)
        m = np.log(mu) - 0.5 * s2
        s = np.sqrt(s2)
        return m, s

    def pdf_LN(self, X, mu, sigma):
        ''' lognormal pdf with actual miu and sigma
        '''
        mu_tmp, sigma_tmp = self.log_params(mu, sigma)
        return sp.stats.lognorm.pdf(X, s=sigma_tmp, scale = np.exp(mu_tmp))

    def logn(self, sigma_cell,Tc,r,n):
        '''  def logn(t,sigma_cell,Tc,r):
             analytic solution for cell only noise
        '''
        mean,sigma = self.log_params(Tc,sigma_cell)
        la = np.log(self.times[n]/(1-r))
        idf = 0.5 * ( 1 + sp.special.erf(  (mean - la) / (np.sqrt(2) * sigma) )  )
        int2 = 0.5 * np.exp( -mean + 0.5 * sigma * sigma ) * sp.special.erfc(  ( -mean + sigma*sigma + la ) / (np.sqrt(2)*sigma) )
        return 1-1*idf+r*idf+self.times[n]*int2



    def pmf_f(self, Tc,r, GF, sigma_cell,sigma_sample):
        """ pmf for the number of labelled cells
            to test: using epsabs=0.1 and epsrel=0.1 in quad might significantly
            speed up the computation without loosing too much precision
        """
        self.Tc = Tc
        self.r = r
        self.GF = GF
        self.sigma_cell = sigma_cell
        self.sigma_sample = sigma_sample

        #P =[sp.integrate.quadrature(self.f, 0.0001, Tc+(sigma_sample*10),args=([i]),tol=1.48e-08, rtol=1.48e-08)[0] for i in range(self.datalen)]
        #P = [sp.integrate.fixed_quad(self.f, 0.0001, Tc+(sigma_sample*10),n=100,args=([i]))[0] for i in range(self.datalen)]
        P = [sp.integrate.fixed_quad(self.f, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200,args=([i]))[0] for i in range(self.datalen)]
        return np.array(P)

    def f(self, TC_,n):
        return sp.stats.binom.pmf(self.data[n], self.ncell[n], self.GF*self.logn(self.sigma_cell,TC_,self.r,n) ) * self.pdf_LN(TC_, self.Tc, self.sigma_sample)
    
    
class dist:#cdef is here to reduce name lookup for __call__
    
    def __init__(self,):
        ''' data = labeling fraction
            times = time for labeling fraction
            ncell = number of cells
        '''
        
        
    
    
    def log_params(self, mu, sigma):
        """ A transformation of paramteres such that mu and sigma are the 
            mean and variance of the log-normal distribution and not of
            the underlying normal distribution.
        """
        s2 = np.log(1.0 + sigma**2/mu**2)
        m = np.log(mu) - 0.5 * s2
        s = np.sqrt(s2)
        return m, s

    def pdf_LN(self, X, mu, sigma):
        ''' lognormal pdf with actual miu and sigma
        '''
        mu_tmp, sigma_tmp = self.log_params(mu, sigma)
        return sp.stats.lognorm.pdf(X, s=sigma_tmp, scale = np.exp(mu_tmp))

    def logn(self, sigma_cell,Tc,r):
        '''  def logn(t,sigma_cell,Tc,r):
             analytic solution for cell only noise
        '''
        mean,sigma = self.log_params(Tc,sigma_cell)
        la = np.log(self.t/(1-r))
        idf = 0.5 * ( 1 + sp.special.erf(  (mean - la) / (np.sqrt(2) * sigma) )  )
        int2 = 0.5 * np.exp( -mean + 0.5 * sigma * sigma ) * sp.special.erfc(  ( -mean + sigma*sigma + la ) / (np.sqrt(2)*sigma) )
        return 1-1*idf+r*idf+self.t*int2



    def pmf_f(self,ncell,Tc,r, GF, sigma_cell,sigma_sample, t, x):
        """ pmf for the number of labelled cells
            to test: using epsabs=0.1 and epsrel=0.1 in quad might significantly 
            speed up the computation without loosing too much precision
        """
        self.ncell = ncell
        self.Tc = Tc
        self.r = r
        self.GF = GF
        self.sigma_cell = sigma_cell
        self.sigma_sample = sigma_sample
        self.t = t
        
        #P =  sp.integrate.quadrature(self.f, 0.01, 11,args=([x]))[0] 
        P = sp.integrate.fixed_quad(self.f, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200,args=([x]))[0]
        return P

    def f(self, TC_,x):
        return sp.stats.binom.pmf(x, self.ncell, self.GF*self.logn(self.sigma_cell,TC_,self.r) ) * self.pdf_LN(TC_, self.Tc, self.sigma_sample)

    def pmf_mean(self,ncell,Tc,r, GF, sigma_cell,sigma_sample, t):
        """ pmf for the number of labelled cells
            to test: using epsabs=0.1 and epsrel=0.1 in quad might significantly 
            speed up the computation without loosing too much precision
        """
        self.ncell = ncell
        self.Tc = Tc
        self.r = r
        self.GF = GF
        self.sigma_cell = sigma_cell
        self.sigma_sample = sigma_sample
        self.t = t
        
        #P =  sp.integrate.quadrature(self.fm, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10))[0] 
        P = sp.integrate.fixed_quad(self.fm, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200)[0]
        return P


    def fm(self, TC_):
        return  self.GF*self.logn(self.sigma_cell,TC_,self.r) * self.pdf_LN(TC_, self.Tc, self.sigma_sample)
    

@np.vectorize
def brdu_model(t,G1, S,G2M , GF,**args):
    S 
    TC = S + G2M + G1
    if t < TC - S:
        return GF * (t + S) / TC
    else:
        return GF

