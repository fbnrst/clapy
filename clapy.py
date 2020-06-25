import numpy as np
from scipy.stats import lognorm,binom
from scipy.special import erf,erfc
from scipy import integrate

import matplotlib.pyplot as plt
from iminuit import Minuit
import io


def log_params( mu, sigma):
    """ A transformation of paramteres such that mu and sigma are the
        mean and variance of the log-normal distribution and not of
        the underlying normal distribution.
    """
    s2 = np.log(1.0 + sigma**2/mu**2)
    m = np.log(mu) - 0.5 * s2
    s = np.sqrt(s2)
    return m, s

def pdf_LN(X, mu, sigma):
    ''' lognormal pdf with actual miu and sigma
    '''
    mu_tmp, sigma_tmp = log_params(mu, sigma)
    return lognorm.pdf(X, s=sigma_tmp, scale = np.exp(mu_tmp))

def logn(sigma_cell,Tc,r,t):
    '''  def logn(t,sigma_cell,Tc,r):
         analytic solution for cell only noise
         sigma_cell : $\hat{\sigma}_c$
         r          : $f_S$
         t          : time
    '''
    s2 = np.log(1.0 + sigma_cell**2/Tc**2)
    mean = np.log(Tc) - 0.5 * s2
    sigma = np.sqrt(s2)
    la = np.log(t/(1-r))
    idf = 0.5 * ( 1 + erf(  (mean - la) / (np.sqrt(2) * sigma) )  )
    int2 = 0.5 * np.exp( -mean + 0.5 * sigma * sigma ) * erfc(  ( -mean + sigma*sigma + la ) / (np.sqrt(2)*sigma) )
    return 1-1*idf+r*idf+t*int2


class asym_lh:
    ''' class with the likelihood for the asymetric cell labelling assays 
        usable for minuit
    '''

    def __init__(self,data,times,ncell):
        ''' data = number of labeld cells
            times = time for labeling fraction
            ncell = number of cells 
        '''
        self.data = np.round(data)
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

        #P =[integrate.quadrature(self.f, 0.0001, Tc+(sigma_sample*10),args=([i]),tol=1.48e-08, rtol=1.48e-08)[0] for i in range(self.datalen)]
        #P = [integrate.fixed_quad(self.f, 0.0001, Tc+(sigma_sample*10),n=100,args=([i]))[0] for i in range(self.datalen)]
        P = [integrate.fixed_quad(self.f, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200,args=([i]))[0] for i in range(self.datalen)]
        return np.array(P)

    def f(self, TC_,n):
        return binom.pmf(self.data[n], self.ncell[n], self.GF*logn(self.sigma_cell,TC_,self.r,self.times[n]) ) * pdf_LN(TC_, self.Tc, self.sigma_sample)

class dist:
    ''' distribution for the asymetric labelling assay   '''

    def __init__(self):
        pass


    def pmf_f(self,ncell,Tc,r, GF, sigma_cell,sigma_sample, t, x):
        """ pmf for the number of labelled cells
              ncell      : number of cells counted
              Tc         : cell cycle length $\tau$
              r          : $f_S$
              GF         : growth fraction $g$
              sigma_cell : $\hat{\sigma}_c$
              t          : time
              x          : number of labelled cells

        """
        self.ncell = ncell
        self.Tc = Tc
        self.r = r
        self.GF = GF
        self.sigma_cell = sigma_cell
        self.sigma_sample = sigma_sample
        self.t = t

        #P =  sp.integrate.quadrature(self.f, 0.01, 11,args=([x]))[0] 
        P = integrate.fixed_quad(self.f, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200,args=([x]))[0]
        return P

    def f(self, TC_,x):
        return binom.pmf(x, self.ncell, self.GF*logn(self.sigma_cell,TC_,self.r,self.t) ) * pdf_LN(TC_, self.Tc, self.sigma_sample)

    def pmf_mean(self,ncell,Tc,r, GF, sigma_cell,sigma_sample, t):
        """ mean number of labelled cells
        """
        self.ncell = ncell
        self.Tc = Tc
        self.r = r
        self.GF = GF
        self.sigma_cell = sigma_cell
        self.sigma_sample = sigma_sample
        self.t = t

        #P =  sp.integrate.quadrature(self.fm, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10))[0] 
        P = integrate.fixed_quad(self.fmean, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200)[0]
        return ncell*P


    def fmean(self, TC_):
        return  self.GF*logn(self.sigma_cell,TC_,self.r,self.t) * pdf_LN(TC_, self.Tc, self.sigma_sample)

    def pmf_std(self,ncell,Tc,r, GF, sigma_cell,sigma_sample, t):
        """ mean number of labelled cells
        """
        self.ncell = ncell
        self.Tc = Tc
        self.r = r
        self.GF = GF
        self.sigma_cell = sigma_cell
        self.sigma_sample = sigma_sample
        self.t = t

        mean = self.pmf_mean(ncell,Tc,r, GF, sigma_cell,sigma_sample, t)
        #P =  sp.integrate.quadrature(self.fm, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10))[0] 
        P = integrate.fixed_quad(self.fstd, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200)[0]
        return np.sqrt(P - mean*mean)

    def fstd(self, TC_):
        p = self.GF*logn(self.sigma_cell,TC_,self.r,self.t)
        n = self.ncell
        res =   (n**2*p**2 - n*p**2 + n*p)* pdf_LN(TC_, self.Tc, self.sigma_sample)
        return np.nan_to_num(res)


    def pmf_skw(self,ncell,Tc,r, GF, sigma_cell,sigma_sample, t):
        """ mean number of labelled cells
        """
        self.ncell = ncell
        self.Tc = Tc
        self.r = r
        self.GF = GF
        self.sigma_cell = sigma_cell
        self.sigma_sample = sigma_sample
        self.t = t
        mean = self.pmf_mean(ncell,Tc,r, GF, sigma_cell,sigma_sample, t)
        std = self.pmf_std(ncell,Tc,r, GF, sigma_cell,sigma_sample, t)
        #P =  sp.integrate.quadrature(self.fm, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10))[0] 
        P = integrate.fixed_quad(self.fskw, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200)[0]
        return (P - 3*mean*std**2 - mean**3)/std**3

    def fskw(self, TC_):
        p = self.GF*logn(self.sigma_cell,TC_,self.r,self.t)
        n = self.ncell
        res =  (n*p - 3*(1 - n)*n*p**2 + (1 - n)*(2 - n)*n*p**3 ) * pdf_LN(TC_, self.Tc, self.sigma_sample)
        return np.nan_to_num(res)

def calc_sigma_true(sigma,fg1,fg2m):
    a = np.sqrt(-(-1 + 3*fg1 - 3*fg1*fg1 + 3*fg2m - 6*fg1*fg2m - 3*fg2m*fg2m + 3*fg1*fg2m*fg2m))
    b = np.sqrt(-(-1 + 3*fg1 - 3*fg1*fg1 + 3*fg2m - 6*fg1*fg2m - 3*fg2m*fg2m))
    return sigma*a/b
    
@np.vectorize
def cla_det_model(t, G1=0.2, S=0.3, G2M=0.5, GF=1, mode=1, **kwargs):
    """ Model for labeling assays in vivo.
        Based on Lefevre et al., 2013 and extended with an initial
        growth fraction.
        
        t    ... time after start of labeling
        S    ... absolute length of S-Phase
        G1   ... absolute length of G1-Phase
        G2M  ... absolute length of G2-Phase and M-Phase
        rmode    ... mean number of daughter cells after cell division remaining
                 in the population
        GF   ... initial growth fraction
        
        Lefevre, J., Marshall, D. J., Combes, A. N., Ju, A. L., Little, M. H.
        & Hamilton, N. A. (2013). Modelling cell turnover in a complex tissue
        during development. Journal of Theoretical Biology, 338, 66-79.
    """
    r = mode
    TC = S + G1 + G2M
    if G2M < 0:
        return sp.nan
    if S + G2M > TC:
        return sp.nan
    else:
        if r==1:
            if t < TC - S:
                return GF * (t + S) / TC
            else:
                return GF
        else:
            # calculate the growth fraction at time t
            g = ( ( GF * r ** (t / TC) ) / ( GF * r ** (t / TC) + (1 - GF) ) )
            if t < G2M:
                return  g * ((r ** ( ( G2M + S ) / TC ) - r ** (( G2M - t ) / TC) ) / (r - 1.0) )
            elif t < TC - S:
                return g * (1.0 - ( r ** ( ( TC + G2M - t ) / TC ) - r ** ( ( G2M  + S) / TC ) ) / (r - 1.0) )
            else:
                return g


def web_fit(times,datas,ncells):
    dat = np.array(datas)
    tim = np.array(times)
    ncell = np.array(ncells)

    Tc_init = np.max(tim)*0.5
    r_init = 0.5
    GF_init = np.mean(np.sort(dat)[-len(dat)//10:])
    Tc_lower = np.min(tim)
    Tc_upper = np.max(tim)
    error = Tc_init*0.1

    lh = asym_lh(dat,tim,ncell)
    mi = Minuit(lh.compute, Tc=Tc_init, r=r_init,GF=GF_init,sigma_sample=Tc_init*0.2,sigma_cell=Tc_init*0.2, \
               error_Tc=error,error_r=0.1,error_GF=0.1,error_sigma_sample=error,error_sigma_cell=error,\
               limit_Tc=(Tc_lower,Tc_upper), limit_r=(0.00001,1),limit_GF=(0,1),limit_sigma_sample=(0.00001,Tc_init),limit_sigma_cell=(0.00001,Tc_init),\
               errordef=0.5,print_level=0)
    mi.migrad();

    fit = dict()
    for i in  mi.values:
        fit.update( {i : {'value' : mi.values[i], '2sigma' : 2*mi.errors[i]}})

    fig = plt.figure(1,figsize=(5,4))


    tf2 = np.linspace(0.01,np.max(tim)*1.1,100)
    d = dist()
    prob = np.zeros(len(tf2))
    nc = np.mean(ncell)
    for t_n,t in enumerate(tf2):
        prob[t_n] = d.pmf_mean(nc,fit['Tc']['value'],fit['r']['value'],fit['GF']['value'],fit['sigma_cell']['value'],fit['sigma_sample']['value'],t)
    plt.plot(tim,dat/ncell,'k.',label='Measurements',zorder=4)
    colorp =  np.array([0.5647058823529412, 0.9333333333333333, 0.5647058823529412]) - np.array([0.4,0.1,0.4])
    plt.plot(tf2,prob/nc,label='probabilistic model',color=colorp,lw=2,zorder=2)
    plt.ylim(0,1.1)
    plt.legend()
    plt.xlabel('time [original units]')
    plt.ylabel('labeling fraction')

    image = io.BytesIO()
    fig.savefig(image, format='png', bbox_inches='tight')
    plt.close(fig)
    return fit,image.getvalue()




from scipy.stats import gamma as gammad
from scipy.special import gamma,gammaincc


def pdf_Gamma(x,mu,sigma):
    beta  = sigma*sigma/mu
    alpha = mu*mu/(sigma*sigma)
    return gammad.pdf(x, a=alpha, scale=beta) 


def loggamma(sigma_cell,Tc,r,t):
    '''  def logn(t,sigma_cell,Tc,r):
         analytic solution for cell only noise
         sigma_cell : $\hat{\sigma}_c$
         r          : $f_S$
         t          : time
    '''
    alpha = Tc*Tc/(sigma_cell*sigma_cell)
    beta  = sigma_cell*sigma_cell/Tc
    ga = gamma(alpha)
    frac = t/(beta-r*beta)
    igaf = ga*gammaincc(alpha,frac)
    
    f1 = beta*ga - r*beta*ga + t*gamma(-1+alpha)*gammaincc(-1+alpha,frac) - beta*igaf + r*beta*igaf
    return f1/(beta*ga) + r



class asym_lhgamma:
    ''' class with the likelihood for the asymetric cell labelling assays 
        usable for minuit
    '''

    def __init__(self,data,times,ncell):
        ''' data = number of labeld cells
            times = time for labeling fraction
            ncell = number of cells 
        '''
        self.data = np.round(data)
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
        #P =[integrate.quadrature(self.f, 0.0001, Tc+(sigma_sample*10),args=([i]),tol=1.48e-08, rtol=1.48e-08)[0] for i in range(self.datalen)]
        #P = [integrate.fixed_quad(self.f, 0.0001, Tc+(sigma_sample*10),n=100,args=([i]))[0] for i in range(self.datalen)]
        P = [integrate.fixed_quad(self.f, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200,args=([i]))[0] for i in range(self.datalen)]
        return np.array(P)

    def f(self, TC_,n):
        res = binom.pmf(self.data[n], self.ncell[n], self.GF*loggamma(self.sigma_cell,TC_,self.r,self.times[n]) ) * pdf_Gamma(TC_, self.Tc, self.sigma_sample)
        return np.nan_to_num(res)

class dist_gamma:
    ''' distribution for the asymetric labelling assay   '''

    def __init__(self):
        pass


    def pmf_f(self,ncell,Tc,r, GF, sigma_cell,sigma_sample, t, x):
        """ pmf for the number of labelled cells
              ncell      : number of cells counted
              Tc         : cell cycle length $\tau$
              r          : $f_S$
              GF         : growth fraction $g$
              sigma_cell : $\hat{\sigma}_c$
              t          : time
              x          : number of labelled cells

        """
        self.ncell = ncell
        self.Tc = Tc
        self.r = r
        self.GF = GF
        self.sigma_cell = sigma_cell
        self.sigma_sample = sigma_sample
        self.t = t

        #P =  sp.integrate.quadrature(self.f, 0.01, 11,args=([x]))[0] 
        P = integrate.fixed_quad(self.f, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200,args=([x]))[0]
        return np.array(P)

    def f(self, TC_,x):
        res = binom.pmf(x, self.ncell, self.GF*loggamma(self.sigma_cell,TC_,self.r,self.t) ) * pdf_Gamma(TC_, self.Tc, self.sigma_sample)
        return np.nan_to_num(res)

    def pmf_mean(self,ncell,Tc,r, GF, sigma_cell,sigma_sample, t):
        """ mean number of labelled cells
        """
        self.ncell = ncell
        self.Tc = Tc
        self.r = r
        self.GF = GF
        self.sigma_cell = sigma_cell
        self.sigma_sample = sigma_sample
        self.t = t

        #P =  sp.integrate.quadrature(self.fm, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10))[0] 
        P = integrate.fixed_quad(self.fmean, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200)[0]
        return ncell*P


    def fmean(self, TC_):
        res =  self.GF*loggamma(self.sigma_cell,TC_,self.r,self.t) * pdf_Gamma(TC_, self.Tc, self.sigma_sample)
        return np.nan_to_num(res)


    def pmf_std(self,ncell,Tc,r, GF, sigma_cell,sigma_sample, t):
        """ mean number of labelled cells
        """
        self.ncell = ncell
        self.Tc = Tc
        self.r = r
        self.GF = GF
        self.sigma_cell = sigma_cell
        self.sigma_sample = sigma_sample
        self.t = t

        mean = self.pmf_mean(ncell,Tc,r, GF, sigma_cell,sigma_sample, t)
        #P =  sp.integrate.quadrature(self.fm, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10))[0] 
        P = integrate.fixed_quad(self.fstd, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200)[0]
        return np.sqrt(P - mean*mean)

    def fstd(self, TC_):
        p = self.GF*loggamma(self.sigma_cell,TC_,self.r,self.t)
        n = self.ncell
        res =   (n**2*p**2 - n*p**2 + n*p)* pdf_Gamma(TC_, self.Tc, self.sigma_sample)
        return np.nan_to_num(res)

    def pmf_skw(self,ncell,Tc,r, GF, sigma_cell,sigma_sample, t):
        """ mean number of labelled cells
        """
        self.ncell = ncell
        self.Tc = Tc
        self.r = r
        self.GF = GF
        self.sigma_cell = sigma_cell
        self.sigma_sample = sigma_sample
        self.t = t
        mean = self.pmf_mean(ncell,Tc,r, GF, sigma_cell,sigma_sample, t)
        std = self.pmf_std(ncell,Tc,r, GF, sigma_cell,sigma_sample, t)
        #P =  sp.integrate.quadrature(self.fm, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10))[0] 
        P = integrate.fixed_quad(self.fskw, max(0.0001,Tc-(sigma_sample*5)), Tc+(sigma_sample*10),n=200)[0]
        return (P - 3*mean*std**2 - mean**3)/std**3

    def fskw(self, TC_):
        p = self.GF*loggamma(self.sigma_cell,TC_,self.r,self.t)
        n = self.ncell
        res =  (n*p - 3*(1 - n)*n*p**2 + (1 - n)*(2 - n)*n*p**3 ) * pdf_Gamma(TC_, self.Tc, self.sigma_sample)
        return np.nan_to_num(res)
