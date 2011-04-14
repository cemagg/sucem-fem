from FenicsCode import Consts

mu0 = Consts.mu0
eps0 = Consts.eps0

class MaterialProperties(object):
    eps_r = 1.0
    mu_r = 1.0

    def __init__(self, eps_r=None, mu_r=None):
        if eps_r is not None: self.set_eps_r(eps_r)
        if mu_r is not None: self.set_mu_r(mu_r)
        
    def set_eps_r(self, eps_r):
        self.eps_r = eps_r

    def set_mu_r(self, mu_r):
        self.mu_r = mu_r

    def get_eps(self):
        return eps0*self.eps_r

    def get_mu(self):
        return mu0*self.mu_r


