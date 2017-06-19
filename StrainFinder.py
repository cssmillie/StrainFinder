import argparse, copy, cPickle, itertools, os.path, random, sys, time, uuid
import numpy as np
import scipy.spatial.distance as ssd
from openopt import NLP, MINLP
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)

bps = {'A':[1,0,0,0], 'C':[0,1,0,0], 'G':[0,0,1,0], 'T':[0,0,0,1]}
nts = np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]])
msg = False
t0 = time.time()
stdout = sys.stdout

class DummyFile(object):
    def write(self, x):
        pass
    


def quiet():
    # Turn off stdout
    sys.stdout = DummyFile()


def loud():
    # Turn on stdout
    sys.stdout = stdout


def message(self, x):
    # Print message to screen
    global msg
    if msg == True:
        sys.stderr.write('(%s) %s\n' %(self.__class__.__name__, x))


def rselect(x):
    # Select random index using weights in x
    if sum(x) == 0:
        x = [1.] * len(x)
    if sum(x) != 1:
        x = norm(x)
    ci = 0
    r = random.random()
    for i in range(len(x)):
        ci += x[i]
        if r <= ci:
            return i


def norm(x):
    # Normalize a vector x so that sum(x) = 1
    sum_x = sum(x)
    return np.array([1.*xi/sum_x for xi in x])


def error(nt, e):
    # Simulate error with probability e
    if random.random() < e:
        return random.choice(nts)
    return nt


def discretize_genotypes(a):
    b = a.max(-1)
    condition = a == b[..., np.newaxis]
    a[condition] = 1
    a[np.logical_not(condition)] = 0
    return a


def gdist(u, v):
    dist = 0
    for i in range(len(u)):
        if u[i] == 'N' or v[i] == 'N' or u[i] == v[i]:
            continue
        else:
            dist += 1
    dist = 1.*dist/len(u)
    return dist


fdist = ssd.cosine

def parse_args():
    
    # Initialize parser
    parser = argparse.ArgumentParser()
    
    # Input options
    group1 = parser.add_argument_group('General')
    group1.add_argument('--em', help='Input EM object', default=None)
    group1.add_argument('--sim', help='Simulate data?', action='store_true', default=False)
    group1.add_argument('--aln', help='Input alignment (numpy)', default=None)
    group1.add_argument('--data', help='Input data object', default=None)
    group1.add_argument('--msg', help='Print messages?', action='store_true', default=False)
    
    # Simulation options
    group2 = parser.add_argument_group('Simulation')
    group2.add_argument('-m', help='Number of samples', type=int, default=None)
    group2.add_argument('-n', help='Number of strains', type=int, default=None)
    group2.add_argument('-l', help='Alignment length', type=int, default=None)
    group2.add_argument('-d', help='Sequencing depth', type=int, default=None)
    group2.add_argument('-e', help='Sequencing error', type=float, default=1e-2)
    group2.add_argument('-u', help='Mutation rate', type=float, default=1.0)
    group2.add_argument('--sparse', help='Sparse frequencies?', action='store_true', default=False)
    group2.add_argument('--phylo', help='Phylo genotypes?', action='store_true', default=False)
    group2.add_argument('--noise', help='Fraction of alignment to disrupt', type=float, default=0)
    
    # Search options
    group3 = parser.add_argument_group('Search')
    group3.add_argument('-N', help='Number of strains to estimate', type=int, default=None)
    group3.add_argument('--random', help='Use random strain genotypes (default = dominant SNPs)', action='store_true', default=False)
    group3.add_argument('--s_reps', help='Number of searches (shallow)', type=int, default=sys.maxint)
    group3.add_argument('--s_iter', help='Number of iterations (shallow)', type=int, default=sys.maxint)
    group3.add_argument('--d_reps', help='Number of searches (deep)', type=int, default=0)
    group3.add_argument('--d_iter', help='Number of iterations (deep)', type=int, default=sys.maxint)
    group3.add_argument('--n_keep', help='Number of searches to keep', type=int, default=0)
    group3.add_argument('--converge', help='Search until convergence', action='store_true', default=False)
    group3.add_argument('--robust', help='Robust EM (pay penalty to use uniform frequencies)?', action='store_true', default=False)
    group3.add_argument('--penalty', help='Penalty for robust EM', type=float, default=1.25)
    group3.add_argument('--exhaustive', help='Exhaustively search genotypes?', action='store_true', default=False)
    group3.add_argument('--reset', help='Reset reps and time', action='store_true', default=False)
    
    # Stop options
    group4 = parser.add_argument_group('Stop')
    group4.add_argument('--dtol', help='Stop when d(log-likelihood) < dtol for ntol iterations', type=float, default=np.nan)
    group4.add_argument('--ftol', help='Stop when f(log-likelihood) < ftol for ntol iterations', type=float, default=np.nan)
    group4.add_argument('--ntol', help='See --dtol and --ftol options', type=int, default=np.nan)
    group4.add_argument('--detect_limit', help='', type=float, default=0.)
    group4.add_argument('--max_reps', help='Max number of searches', type=int, default=np.nan)
    group4.add_argument('--max_time', help='Max time (seconds)', type=float, default=np.nan)
    group4.add_argument('--min_reps', help='Min number of searches (for convergence)', type=int, default=0)
    group4.add_argument('--min_gdist', help='Min distance between estimated genotypes (p)', type=float, default=np.nan)
    group4.add_argument('--min_fdist', help='Min distance between estimated frequencies (z)', type=float, default=np.nan)
    
    # Write options
    group5 = parser.add_argument_group('Write')
    group5.add_argument('--log', help='Log file', default=None)
    group5.add_argument('--em_out', help='Output EM object', default=None)
    group5.add_argument('--aln_out', help='Output alignment (numpy array)', default=None)
    group5.add_argument('--data_out', help='Output data object', default=None)
    group5.add_argument('--otu_out', help='Output OTU table', default=None)
    group5.add_argument('--merge_out', help='Merge estimates with output file?', action='store_true', default=False)
    group5.add_argument('--force_update', help='Force update?', action='store_true', default=False)
    
    # Parse arguments
    args = parser.parse_args()
    
    global msg
    msg = args.msg
    
    return args


class Data():
    
    def __init__(self, sim=None, m=None, n=None, l=None, d=None, u=None, x=None, p=None, z=None, e=0., sparse=None, phylo=None):
        
        # initialize
        self.sim = sim
        self.m = m
        self.n = n
        self.l = l
        self.d = d
        self.u = u
        self.x = x
        self.p = p
        self.z = z
        self.e = e
        self.nt_freq = None
        self.shuffled = None
        self.sparse = sparse
        self.phylo = phylo
        
        # get dimensions
        if self.x is not None:
            self.m, self.l = np.shape(self.x)[:2]
        
        # load data
        if self.sim:
            self = self.simulate()
        
        # nucleotide frequencies
        self = self.calc_nt_freq()
    
    
    def simulate(self):
        message(self, 'Simulating random dataset')
        
        # draw random frequencies (m x n)
        if self.sparse == False:
            self = self.random_z()
        elif self.sparse == True:
            self = self.sparse_z()
        
        # draw random genotypes (n x l x 4)
        if self.phylo == False:
            self = self.random_p()
        elif self.phylo == True:
            self = self.phylo_p()
        
        # generate alignment (m x l x 4)
        self.x = self.random_x()
        
        return self
    
    
    def random_z(self):
        message(self, 'Generating random strain frequencies (%d x %d)' %(self.n, self.m))
        self.z = np.array([np.random.dirichlet(np.ones(self.n)) for i in range(self.m)])
        return self
    
    
    def sparse_z(self):
        message(self, 'Generating random strain frequencies (%d x %d)' %(self.n, self.m))
        z = []
        for i in range(self.m):
            # get number of non-zero genotypes
            n = random.choice(range(2, self.n))
            # draw frequencies from dirichlet
            d = np.random.dirichlet([1]*n)
            # randomly select zeros
            q = [random.choice(range(n+1)) for j in range(self.n-n)]
            # insert zeros in frequency matrix
            z.append(np.insert(d, sorted(q), [0]*(self.n-n)))
        self.z = np.array(z)
        return self
    
    
    def random_p(self):
        message(self, 'Generating random genotypes (%d x %d)' %(self.n, self.l))
        self.p = np.array([[random.choice(nts) for j in range(self.l)] for i in range(self.n)])
        return self
    
    
    def majority_p(self, k=None):
        
        # select random number of strains
        if k is None:
            k = random.randint(1, min(self.m, self.n))
        message(self, 'Guessing initial strain genotypes (%d x %d) from dominant SNPs in %d samples' %(self.n, self.l, k))
        
        # generate random genotypes (n x l x 4)
        self.p = np.array([[random.choice(nts) for j in range(self.l)] for i in range(self.n)])
        
        # get dominant snps in each sample (m x l x 4)
        p = nts[self.data.x.argmax(axis=2)]
        
        # select random strain indices and replace
        i = random.sample(range(self.n), k)
        j = random.sample(range(self.m), k)
        self.p[i,:,:] = p[j,:,:]
        
        return self
    

    def weighted_p(self, k=None):

        # select random number of strains
        if k is None:
            k = random.randint(1, min(self.m, self.n))
        message(self, 'Guessing initial strain genotypes (%d x %d) from random SNPs in %d samples' %(self.n, self.l, k))

        # generate random genotypes (n x l x 4)
        self.p = np.array([[random.choice(nts) for j in range(self.l)] for i in range(self.n)])

        # select k random strains from samples
        i = np.array(random.sample(range(self.m), self.m))[[j % self.m for j in range(k)]]
        p = np.apply_along_axis(lambda x: nts[rselect(x)], 2, self.data.x[i,:,:])

        # select random strain indices and replace
        i = random.sample(range(self.n), k)
        self.p[i,:,:] = p
        print p
        
        return self
    
    
    def phylo_p(self):
        message(self, 'Generating random genotypes (%d x %d)' %(self.n, self.l))
        
        # load dendropy
        import dendropy
        
        # make alignment from random tree
        tree = dendropy.treesim.pure_kingman(dendropy.TaxonSet(map(str, range(self.n))))
        seqs = [si for si in dendropy.seqsim.generate_hky_dataset(seq_len = self.l*10, tree_model = tree, mutation_rate = self.u).as_string('fasta').split('\n') if si != '' and si[0] != '>']
        
        # count k-morphic sites
        counts = np.array([len(set([seqs[i][j] for i in range(self.n)])) for j in range(self.l*10)])
        
        # look at 2-3 morphic sites
        sites = [i for i in range(len(counts)) if 1 < counts[i] < 4][:self.l]
        
        # construct p
        self.p = np.array([[bps[seqs[i][j]] for j in sites] for i in range(self.n)])
        
        return self
    
    
    def random_x(self):
        message(self, 'Simulating random alignment (%d x %d)' %(self.m, self.l))
        return np.array([[sum([error(self.p[rselect(self.z[i])][j], self.e) for k in range(self.d)]) for j in range(self.l)] for i in range(self.m)])
    
    
    def resample_x(self):
        message(self, 'Resampling alignment from estimated nucleotide frequencies')
        
        # calculate sequencing depth in each position
        total = self.x.sum(axis=2)
        
        # select k nucleotides from freq at each position
        x = np.array([[np.array([0,0,0,0])+sum([nts[rselect(self.x[i,j,:])] for k in range(int(total[i,j]))]) for j in range(self.l)] for i in range(self.m)])
        
        return x
    
    
    def add_noise(self, f):
        message(self, 'Replacing %.2f%% of the alignment with [A,C,G,T] @ random frequencies' %(100*f))
        self.shuffled = np.zeros([self.m, self.l], dtype=bool)
        for i in range(self.m):
            for j in range(self.l):
                if random.random() <= f:
                    depth = sum(self.x[i,j,:])
                    self.x[i,j,:] = (depth * np.random.dirichlet([1,1,1,1])).round()
                    self.shuffled[i,j] = True
        return self
    
    
    def get_genotypes(self):
        acgt = np.array('A C G T'.split())
        seqs = acgt[np.where(self.p == 1)[2]].reshape(self.n, self.l)
        seqs = map(lambda a: ''.join(a), seqs)
        return seqs   
    
    
    def calc_nt_freq(self):
        total = np.maximum(1., self.x.sum(axis=2))
        self.nt_freq = ((1-self.e)*(np.divide(self.x, 1.*total[:,:,np.newaxis])).clip(1e-10) + self.e/4.).clip(1e-10)
        return self
    
    
    def write_aln(self, out_fn):
        if out_fn:
            message(self, 'Writing alignment to "%s"' %(out_fn))
            cPickle.dump(self.x, open(out_fn, 'wb'), protocol=2)
    
    
    def write_data(self, out_fn):
        if out_fn:
            message(self, 'Writing data object to "%s"' %(out_fn))
            cPickle.dump(self, open(out_fn, 'wb'), protocol=2)
    
    


class Estimate(Data):
    
    def __init__(self, data_obj, n, p=None, z=None, random=False, e=.01, robust=False, penalty=None):
        
        self.data = data_obj # alignment data
        self.x = self.data.x # alignment (M,L,4)
        self.m = self.data.m # number of subjects
        self.l = self.data.l # alignment length
        self.n = n # number of strains
        self.p = p # strain genotypes (N,L,4)
        self.z = z # strain frequencies (M,N)
        self.e = e # error rate
        self.random = random
        self.loglik = None # current log-likelihood
        self.logliks = [] # past log-likelihoods
        self.aic = None # current aic
        self.aics = [] # past aics
        self.bic = None
        self.bics = []
        self.update = True # update estimate? (bool)
        self.robust = robust # robust estimation? (bool)
        self.mask = np.ones([self.m, self.l], dtype=bool) # masked sites (M, L)
        self.penalty = penalty # robust penalty
        self.uid = uuid.uuid4() # unique id
        
        # Random guess
        if self.z is None:
            self = self.random_z()
        if self.p is None:
            if self.random == True:
                self = self.random_p()
            else:
                self = self.majority_p()
        
        # Log-likelihood
        self = self.calc_likelihood()
    
    
    def get_genotypes(self, detect_limit=0):
        
        # Total counts at every alignment site (M,L)
        u = np.einsum('ij,ij->ij', self.mask, self.x.sum(axis=2))
        
        # Expected counts for each strain at each alignment site (N,L)
        v = np.einsum('ij,ik->jk', self.z, u)
        
        # Strain genotypes (N,L)
        acgt = np.array('A C G T'.split())
        w = acgt[np.where(self.p == 1)[2]].reshape(self.n, self.l)
        
        # Mask strain genotypes
        w[v <= detect_limit] = 'N'
        
        # Collapse genotypes
        w = map(lambda a: ''.join(a), w)
        
        return w
    
    
    def calc_likelihood(self):
        # Site likelihoods (M,L)
        l1 = self.calc_site_likelihoods()
        # Penalize masked sites
        if self.robust == True:
            i = np.logical_not(self.mask)
            l1[i] = self.penalty * self.calc_site_likelihoods(maxent=True)[i]
        # Update likelihoods
        self.loglik = np.sum(l1)
        self.logliks.append(self.loglik)
        message(self, 'Log-likelihood is %f' %(self.loglik))
        return self    
    
    
    def calc_site_likelihoods(self, optimal=False, maxent=False):
        # Calculate site likelihoods (M x L)
        if optimal == True:
            self = self.calc_nt_freq()
            return np.einsum('ijk,ijk->ij', self.x, np.log(((1-self.e)*self.nt_freq + (self.e/4.)).clip(1e-10)))
        if maxent == True:
            return (np.log(.25)*self.x).sum(axis=2)
        else:
            return np.einsum('ijk,ijk->ij', self.x, np.log(((1-self.e)*np.einsum('ij...,j...k->i...k', self.z, self.p) + (self.e/4.)).clip(1e-10)))
    
    
    def calc_aic(self):
        # Calculate AIC
        penalty = 1 + self.m*(self.n - 1) + self.n*self.l*3
        self.aic = 2*penalty - 2*self.loglik
        if not hasattr(self, 'aics'):
            self.aics = []
        self.aics.append(self.aic)
        message(self, 'AIC is %f' %(self.aic))
        return self
    
    
    def calc_bic(self):
        # Calculate BIC
        pp = 1 + self.m*(self.n - 1) + self.n*self.l*3
        dd = self.x.sum().sum().sum()
        self.bic = pp*np.log(dd) - 2*self.loglik
        if not hasattr(self, 'bics'):
            self.bics = []
        self.bics.append(self.bic)
        message(self, 'BIC is %f' %(self.bic))
        return self

    def exhaustive_search_p(self, c):
        message(self, 'Running exhaustive search of strain genotypes')
        
        # Calculate nucleotide frequencies @ j, dim = (M,C,4)
        a1 = lambda j: (1-self.e)*np.einsum('ij,jkl->ikl', self.z, c) + (self.e/4.)
        
        # Calculate site likelihoods @ j, dim = (M,C)
        l1 = lambda j: np.einsum('ik,ijk->ij', self.x[:,j,:], np.log(a1(j).clip(1e-10)))
        # Calculate alternative likelihood @ j, dim = (M)
        l2 = lambda j: (np.log(.25)*self.x[:,j,:]).sum(axis=1)
        
        if self.robust == False:
            # Get index of maximum likelihood genotypes
            ci = lambda j: np.argmax(np.sum(l1(j), axis=0))
        
        elif self.robust == True:
            def ci(j):
                # Get index of penalized maximum likelihood genotypes
                l1j = l1(j)
                l2j = l2(j)
                i = np.argmax(np.sum(np.maximum(l1j, self.penalty*l2j[:,np.newaxis]), axis=0))
                # Update mask
                self.mask[:,j] = (l1j[:,i] >= self.penalty*l2j)
                return i
            
        # Exhaustive search function
        self.p = np.swapaxes(np.array([c[:,ci(j),:] for j in range(self.l)]), 0, 1)
        
        return self    
    
    
    def max_loglik_p(self, method='scipy_slsqp'):
        
        message(self, 'Optimizing strain genotypes using %s' %(method))
        
        # Initialize frequencies
        p = []
        
        # Bound genotypes in (0,1)
        lb = np.zeros(4*self.n)
        ub = np.ones(4*self.n)
        
        # Constrain genotypes to sum to 1
        def h(a):
            return np.reshape(a, (self.n, 4)).sum(axis=1) - 1.0
        
        quiet()
        
        # Optimize genotypes at every position
        for j in range(self.l):
            
            # Copy mask
            mj = self.mask[:,j].copy()
            
            # Objective function
            def f(a):
                # Reshape strain genotypes
                pi = np.reshape(a, (self.n,4))
                # Nucleotide frequencies at position j
                a1 = (1-self.e)*(np.dot(self.z, pi)) + (self.e/4.)
                # Site likelihoods
                l1 = (self.x[:,j,:] * np.log(a1.clip(1e-10))).sum(axis=1)
                # Robust estimation
                if self.robust == True:
                    # Alternative likelihoods
                    l2 = (np.log(.25)*self.x[:,j,:]).sum(axis=1)
                    # Pay likelihood penalty to mask sites
                    i = l1 >= self.penalty*l2
                    # Update mask
                    mj[i] = True
                    mj[np.logical_not(i)] = False
                    # Penalized likelihood
                    lf = l1[i].sum() + self.penalty*l2[np.logical_not(i)].sum()
                else:
                    # Normal likelihood
                    lf = l1.sum()
                # L2 penalty
                #return -1.*lf - pi.max(axis=1).sum()
                return -1.*lf - (pi**2).sum()
            
            # Calculate original likelihood
            x0 = self.p[:,j,:].flatten()
            l0 = f(x0)
            
            # Optimize genotypes
            g = [.25] * 4 * self.n
            soln = NLP(f, g, h=h, lb=lb, ub=ub, gtol=1e-5, contol=1e-5, name='NLP1').solve(method, plot=0)
            
            # Update genotypes
            if soln.ff <= l0 and soln.isFeasible == True:
                # Discretize results
                xf = discretize_genotypes(np.reshape(soln.xf.clip(0,1), (self.n, 4)))
                lf = f(xf.flatten())
                # Validate likelihood
                if lf < l0:
                    # Update genotypes and mask
                    p.append(xf)
                    if self.robust == True:
                        self.mask[:,j] = mj
                    continue
            # If likelihood not improved, use original genotypes
            p.append(self.p[:,j,:])
        
        loud()
        
        # Fix shape
        self.p = np.swapaxes(np.array(p), 0, 1)
        return self
    
    
    def max_loglik_z(self, method='scipy_slsqp'):
        message(self, 'Optimizing strain frequencies using %s' %(method))
        
        # Initialize frequencies
        z = []
        
        # Bound frequencies in (0,1)
        lb = np.zeros(self.n)
        ub = np.ones(self.n)
        
        # Constrain frequencies to sum to 1
        def h(a):
            return sum(a) - 1
        
        quiet()
        
        # For every subject
        for i in range(self.m):
            
            # Objective function
            def f(a):
                # Get frequencies (N) and error rate
                zi = a; ei = self.e;
                # Get nucleotide frequencies at every position (L,4)
                a1 = np.einsum('i...,i...', zi, self.p)[self.mask[i]]
                # Remove masked sites from alignment (L,4)
                xi = self.x[i][self.mask[i]]
                # Error correct and weight by counts
                a2 = np.einsum('ij,ij', xi, np.log(((1-ei)*a1 + ei/4.).clip(1e-10)))
                # Return negative log-likelihood
                return -1.*a2
            
            # Run optimization
            g = self.z[i,:]
            soln = NLP(f, g, lb=lb, ub=ub, h=h, gtol=1e-5, contol=1e-5, name='NLP1').solve(method, plot=0)
            if soln.ff <= f(g) and soln.isFeasible == True:
                zi = soln.xf
                z.append(zi.clip(0,1))
            else:
                z.append(g)
            
        loud()
        
        # Update frequencies and error rate
        self.z = np.array(z)
        
        return self
    
    
    def run_em(self, n_iter=None, c=None, dtol=None, ftol=None, ntol=None, max_time=None, exhaustive=False):
        
        # Run EM algorithm and quit on (1) n_iter, (2) max_time, or (3) tol/ntol
        message(self, 'Running %d iterations of EM algorithm' %(n_iter))
        
        for i in xrange(n_iter):
            
            # Check time
            if time.time() - t0 >= max_time:
                break
            
            # Check convergence
            if self.local_convergence(dtol, ftol, ntol):
                break
            
            # Optimize genotypes and frequencies
            self = self.max_loglik_z()
            if exhaustive == True:
                self = self.exhaustive_search_p(c)
            else:
                self = self.max_loglik_p()
            
            # Update likelihoods
            self = self.calc_likelihood()
            self = self.calc_aic()
            self = self.calc_bic()
        
        return self
    
    
    def local_convergence(self, dtol=np.nan, ftol=np.nan, ntol=np.nan):
        
        # Check if tolerance is set
        if (np.isnan(dtol) and np.isnan(ftol)) or np.isnan(ntol):
            return False
        
        # Check if the estimate has run for ntol iterations
        if len(self.logliks) <= ntol:
            return False
        
        # If the absolute change is less than tol, estimate has converged
        if not np.isnan(dtol):
            l0 = self.logliks[-ntol]
            l1 = self.logliks[-1]
            if l1 - l0 <= dtol:
                message(self, 'Estimate has converged (loglik0 = %f, loglik1 = %f)' %(l0, l1))
                return True
        
        # If the percent change is less than tol, estimate has converged
        if not np.isnan(ftol):
            d0 = self.logliks[-1] - self.logliks[0]
            d1 = self.logliks[-1] - self.logliks[-ntol]
            if d0 > 0 and 1.*d1/d0 <= ftol:
                message(self, 'Estimate has converged (dloglik0 = %f, dloglik1 = %f, ratio = %f)' %(d0, d1, 1.*d1/d0))
                return True
        
        return False
    
    


class EM():
    
    def __init__(self, data, estimates=[]):
        
        # Initialize data
        self.data = data
        self.estimates = np.array(estimates)
        
        # Progress tracking
        self.r0 = None
        self.t0 = None
        self.total_reps = 0
        self.total_time = 0
    
    
    def current_reps(self):
        return self.total_reps + self.r0
    
    
    def current_time(self):
        return time.time() - self.t0
    
    
    def add_estimate(self, estimate, i=None):
        if i is None:
            message(self, 'Appending estimate')
            self.estimates = np.append(self.estimates, estimate)
        else:
            message(self, 'Replacing estimate %d' %(i))
            self.estimates[i] = estimate
        return self
    
    
    def clear_estimates(self, keep_n=None):
        n0 = len(self.estimates)
        i = [index for index, estimate in enumerate(self.estimates) if estimate.n == keep_n]
        self.estimates = self.estimates[i]
        n1 = len(self.estimates)
        message(self, 'Cleared %d estimates' %(n0 - n1))
        return self
    
    
    def select_best_estimates(self, n_keep=None):
        if n_keep is None:
            n_keep = len(self.estimates)
        message(self, 'Selecting %d best estimates' %(n_keep))
        l = np.array([estimate.loglik for estimate in self.estimates])
        i = np.array([a for a in l.argsort() if not np.isnan(l[a])])[-n_keep:]
        return self.estimates[i]
    
    
    def update_best_estimates(self, n_keep=None):
        
        # Collapse estimates with the same uid
        estimates = {}
        for estimate in self.estimates:
            uid = estimate.uid
            if uid not in estimates:
                estimates[uid] = estimate
            else:
                if estimate.loglik >= estimates[uid].loglik:
                    estimates[uid] = estimate
        self.estimates = np.array(estimates.values())
        
        # Select best estimates
        self.estimates = self.select_best_estimates(n_keep)
        return self
    
    
    def shallow_search(self, n, n_reps=sys.maxint, n_iter=sys.maxint, n_keep=None, c=None, exhaustive=False, random=False, robust=False, penalty=None, dtol=None, ftol=None, ntol=None, max_reps=sys.maxint, max_time=sys.maxint, log_fn=None, out_fn=None):
        message(self, 'Running shallow search')
        
        # Quickly search initial conditions
        for i in xrange(n_reps):
            
            # Check reps
            if self.current_reps() >= max_reps:
                self.write_log(message='%s\treps=%d' %(out_fn, self.current_reps()), log_fn=log_fn)
                break

            # Check time
            if time.time() - t0 >= max_time:
                break
            
            # Initialize estimate and run EM
            estimate = Estimate(self.data, n=n, robust=robust, penalty=penalty, random=random)
            estimate = estimate.run_em(n_iter, c=c, dtol=dtol, ftol=ftol, ntol=ntol, max_time=max_time, exhaustive=exhaustive)
            self = self.add_estimate(estimate)
            self.r0 += 1
        
        # Select best estimates by log-likelihood
        self = self.update_best_estimates(n_keep)
        
        return self
    
    
    def deep_search(self, n, n_reps=1, n_iter=sys.maxint, n_keep=None, c=None, exhaustive=False, dtol=None, ftol=None, ntol=None, max_time=sys.maxint, log_fn=None, out_fn=None):
        
        # Get indices of estimates for deep search
        if len(self.estimates) == 0:
            return self
        if n_reps < 0:
            order = sorted(random.sample(range(len(self.estimates)), abs(n_reps)))
        else:
            order = range(len(self.estimates)) * n_reps
        random.shuffle(order)
        
        message(self, 'Running deep search from %d initial conditions' %(len(order)))
        # For each estimate, run EM and update estimates
        for i in order:
            
            # Check time
            if time.time() - t0 >= max_time:
                break

            estimate = self.estimates[i]
            new_estimate = copy.copy(self.estimates[i])
            new_estimate.update = True
            new_estimate = new_estimate.run_em(n_iter, c, dtol=dtol, ftol=ftol, ntol=ntol, max_time=max_time, exhaustive=exhaustive)
            if new_estimate.loglik > estimate.loglik:
                message(self, 'Updating estimate (loglik0 = %f, loglik1 = %f)' %(estimate.loglik, new_estimate.loglik))                
                self = self.add_estimate(new_estimate, i)
        
        # Select best estimates by log-likelihood
        self = self.update_best_estimates(n_keep)
        
        return self
    
    
    def converge_search(self, n, n_keep=None, c=None, exhaustive=False, robust=False, penalty=None, random=False, dtol=None, ftol=None, ntol=None, max_reps=None, max_time=None, log_fn=None, out_fn=None):
        # Refine estimates (run until convergence)
        self = self.deep_search(n=n, c=c, exhaustive=exhaustive, dtol=dtol, ftol=ftol, ntol=ntol, max_time=max_time, log_fn=log_fn, out_fn=out_fn)
        # New estimates (run until convergence)
        self = self.shallow_search(n=n, n_keep=n_keep, c=c, exhaustive=exhaustive, random=random, robust=robust, penalty=penalty, dtol=dtol, ftol=ftol, ntol=ntol, max_reps=max_reps, max_time=max_time, log_fn=log_fn, out_fn=out_fn)
        return self
    
    
    def genetic_distances(self, detect_limit=0, use_true=False):
        if use_true == True:
            ps = [self.data.get_genotypes(), self.select_best_estimates(1)[0].get_genotypes(detect_limit=detect_limit)]
        else:
            ps = [e.get_genotypes(detect_limit=detect_limit) for e in self.estimates]
        dists = []
        for [p1, p2] in itertools.combinations(ps, 2):
            n = len(p1)
            dist = np.zeros([n, n])
            dist[:] = np.nan
            for i in range(n):
                for j in range(n):
                    dist[i][j] = gdist(p1[i], p2[j])
            dists.append(dist.min(axis=1))
            dists.append(dist.min(axis=0))
        return dists
    
    
    def frequency_distances(self, use_true=False):
        if use_true == True:
            zs = [self.data.z, self.select_best_estimates(1)[0].z]
        else:
            zs = [estimate.z for estimate in self.estimates]    
        dists = []
        for [z1, z2] in itertools.combinations(zs, 2):
            n = z1.shape[1]
            dist = np.zeros([n, n])
            dist[:] = np.nan
            for i in range(n):
                for j in range(n):
                    dist[i][j] = abs(fdist(z1[:,i], z2[:,j]))
            dists.append(dist.min(axis=1))
            dists.append(dist.min(axis=0))
        return dists
    
    
    def global_convergence(self, min_reps=0, min_gdist=None, min_fdist=None, detect_limit=0, log_id=None, log_fn=None):
        
        # Check arguments
        if np.isnan(min_gdist) or np.isnan(min_fdist) or len(self.estimates) == 0:
            return False

        # Number of searches
        if self.current_reps() < min_reps:
            return False
        
        # Distances between estimates
        u = max(map(np.mean, self.frequency_distances()))
        v = max(map(np.mean, self.genetic_distances(detect_limit=detect_limit)))
        
        # Convergent -> return True
        if u <= min_fdist and v <= min_gdist:
            self.write_log(message='%s\treps=%d,fdist=%f,gdist=%f' %(log_id, self.current_reps(), u, v), log_fn=log_fn)
            return True
        
        # Divergent -> return False
        return False
    
        
    def write_log(self, message, log_fn):
        # Write message to log file
        if message and log_fn:
            os.system('touch %s' %(log_fn))
            fh = open(log_fn, 'a')
            fh.write('%s\n' %(message))
            fh.close()
    
    
    def fix_references(self):
        # Make sure all references point to the same object to reduce file size
        message(self, 'Updating references to data object')
        for estimate in self.estimates:
            estimate.data = self.data
            estimate.x = self.data.x
        return self
    
    
    def merge_estimates(self, em, n_keep=1, reset=False):
        # Select the n_keep best estimates from 2 EM objects and fix references
        message(self, 'Merging EM objects')
        self.estimates = np.concatenate([em.estimates, self.estimates])
        self = self.update_best_estimates(n_keep)
        if reset == False:
            self.total_reps = em.total_reps + self.r0
            self.total_time = em.total_time + self.t0
        else:
            self.total_reps = self.r0
            self.total_time = self.t0
        self = self.fix_references()
        return self
    
    
    def check_update(self):
        updates = [estimate.update for estimate in self.estimates]
        if True in updates:
            return True
        else:
            return False        
    
    
    def write_em(self, out_fn, n_keep=1, merge_out=False, force_update=False, reset=False):
        
        if merge_out == True and os.path.exists(out_fn):
            em = cPickle.load(open(out_fn, 'rb'))
            self = self.merge_estimates(em, n_keep=n_keep, reset=reset)
        else:
            self.total_reps += self.r0
            self.total_time += self.t0
                
        if out_fn and (self.check_update() or force_update == True):
            for estimate in self.estimates:
                estimate.update = False
            cPickle.dump(self, open(out_fn, 'wb'), protocol=2)
    
    
    def write_otu_table(self, out_fn, detect_limit=0.):
        # Write otu table (no rownames, colnames = genotypes)
        np.savetxt(out_fn, self.z, delimiter='\t', comments='', header='\t'.join(self.get_genotypes(detect_limit=detect_limit)))
    
    


def load_em(args):
    
    # Load existing EM object
    if args.em and os.path.exists(args.em):
        message(str, 'Loading EM object from file "%s"' %(args.em))
        em = cPickle.load(open(args.em, 'rb'))
    
    # Make new EM object
    else:
        if args.data and os.path.exists(args.data):
            data = cPickle.load(open(args.data, 'rb'))
        elif args.aln and os.path.exists(args.aln):
            data = Data(x=cPickle.load(open(args.aln, 'rb')))
        elif args.sim:
            data = Data(sim=args.sim, m=args.m, n=args.n, l=args.l, d=args.d, u=args.u, e=args.e, sparse=args.sparse, phylo=args.phylo)
            data = data.add_noise(f=args.noise)
        else:
            quit()
        em = EM(data=data)
    
    # Reset reps and time
    if args.reset == True:
        em.total_reps = 0
        em.total_time = 0
    em.r0 = 0
    em.t0 = time.time()
    
    return em


def write_results(args, em, detect_limit=0, reset=False):
    # Write alignment, data, and EM object
    if args.aln_out:
        em.data.write_aln(out_fn=args.aln_out)
    if args.data_out:
        em.data.write_data(out_fn=args.data_out)
    if args.em_out:
        n_keep = args.n_keep
        em.write_em(out_fn=args.em_out, merge_out=args.merge_out, n_keep=n_keep, force_update=args.force_update, reset=reset)
    if args.otu_out:
        em.write_otu_table(out_fn=args.otu_out, detect_limit=detect_limit)
    



def run():
    
    # Parse command line arguments
    args = parse_args()
        
    # Get EM object
    em = load_em(args)

    # Clear estimates
    em = em.clear_estimates(keep_n=args.N)
    
    # Enumerate nucleotide search space
    if args.exhaustive == True:
        c = np.swapaxes(np.array([combo for combo in itertools.product(nts, repeat=args.N)]), 0, 1)
    else:
        c = None
    
    # Check convergence
    if em.global_convergence(min_reps=args.min_reps, min_gdist=args.min_gdist, min_fdist=args.min_fdist, detect_limit=args.detect_limit, log_id=args.em_out, log_fn=args.log):
        quit()
    
    # Shallow search
    if args.s_reps != 0 and args.converge == False:
        em = em.shallow_search(n=args.N, n_reps=args.s_reps, n_iter=args.s_iter, n_keep=args.n_keep, c=c, robust=args.robust, penalty=args.penalty, random=args.random, exhaustive=args.exhaustive, dtol=args.dtol, ftol=args.ftol, ntol=args.ntol, max_reps=args.max_reps, max_time=args.max_time, log_fn=args.log, out_fn=args.em_out)
    
    # Deep search
    if args.d_reps != 0 and args.converge == False:
        em = em.deep_search(n=args.N, n_reps=args.d_reps, n_iter=args.d_iter, n_keep=args.n_keep, c=c, exhaustive=args.exhaustive, dtol=args.dtol, ftol=args.ftol, ntol=args.ntol, max_time=args.max_time, log_fn=args.log, out_fn=args.em_out)
    
    # Converge search
    if args.converge == True:
        em = em.converge_search(n=args.N, c=c, exhaustive=args.exhaustive, robust=args.robust, penalty=args.penalty, random=args.random, dtol=args.dtol, ftol=args.ftol, ntol=args.ntol, max_reps=args.max_reps, max_time=args.max_time, log_fn=args.log, out_fn=args.em_out)
    
    # Check convergence
    em.global_convergence(min_reps=args.min_reps, min_gdist=args.min_gdist, min_fdist=args.min_fdist, detect_limit=args.detect_limit, log_id=args.em_out, log_fn=args.log)
    
    # Write results
    write_results(args, em, detect_limit=args.detect_limit, reset=args.reset)


if __name__ == '__main__':
    run()

