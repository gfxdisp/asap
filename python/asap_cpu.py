import numpy as np
import scipy
import random
import networkx as nx
import scipy.stats
import copy

class ASAP():

    def __init__(self,N, selective_eig = False, approx = False):
        self.approx = approx
        # Initialize trueskill solver with the number of conditions
        self.ts_solver = TrueSkillSolver(N, approx)
        self.selective_eig = selective_eig


    def get_scores(self):
        '''
        Get current score estimation
        '''
        
        return self.ts_solver.Ms, np.sqrt(self.ts_solver.Vs)
        
    def unroll_mat(self,M):
        '''
        Function to convert matrix M with M[ii][jj] = (number of time ii was chosen over jj)
        to a matrix G with number of rows equal to the total number of comparisons performed.
        Each row of this matrix contains three elements, [ii, jj, 1/0], i.e. conditions participating
        in a single comparison and the outcome of this comparison.
        '''
        N = np.shape(M)[0]
        G = np.empty((0,2))
        for ii in range(0,N):
            for jj in range(0,N):
                G = np.append(G, np.tile(np.array([ii,jj]),(M[ii][jj],1)),0)
                M[ii][jj] = 0
        G = G.astype(int)
        return G

    def compute_minimum_spanning_tree(self,inf_mat):
        '''
        Given an information gain matrix, we want to extract a set of comparisons to perform that
        would have the largest total information gain and would form a connected graph of comparisons.
        The function takes as input an information gain matrix, then computes its reciprocal (1/inf_gain)
        and extracts a minimum spanning tree from it.
        '''
        
        inf_mat = inf_mat+inf_mat.T
        inf_mat[inf_mat<=0] = np.inf
        inf_mat = 1/inf_mat

        GrMST=nx.from_numpy_matrix(inf_mat)
        T=nx.minimum_spanning_tree(GrMST)

        pairs_to_compare = np.asarray(T.edges())
        
        edges=sorted(T.edges(data=True), key=lambda t: t[2].get('weight', 1))
        pairs_to_compare = np.array([t[0:2] for t in edges ])
        
        
        return pairs_to_compare

    def run_asap(self, M, mst_mode=True):
        '''
        The main function to generate pairs for comparisons from the pairwise comparison matrix M
        ''' 
        
        N = np.shape(M)[0]
        
        G = self.unroll_mat(M.copy())
        
        
        # Returns an information gain matrix and the next pair of comparisons with the largest information gain
        inf_mat, pairs_to_compare = self.compute_information_gain_mat(N,G)
        # If we are interested in the set of (N-1) pairs to compare instead of a single comparison, we extract
        # this set with a minimum spanning tree of the inverse of the information gain matrix
        if mst_mode:
            
            # In the first iteration populate information matrix at random
            if M.sum()<(N-1):
                
                inf_mat = np.tril(np.random.rand(N,N))
                np.fill_diagonal(inf_mat,0)
            
                
            pairs_to_compare = self.compute_minimum_spanning_tree(inf_mat)
        
        
        return pairs_to_compare

    
    
    def compute_prob_cmps(self):
        '''
        prob: matrix with probability of one condition chosen over another with prob[ii][jj] computed from 
        CDF(mu_ii-mu_jj, 1+var_ii+var_jj). Used to weight the computation of the expected information gain
        prob_cmp: matrix with probability of performing an evaluation of expected information gain for each 
        pair of conditions
        '''
        means, vrs = self.ts_solver.Ms, self.ts_solver.Vs
        N = np.shape(means)[0]

        diff_means = np.tile(means,(N,1))-np.tile(means,(N,1)).T
        vars_sum = 1+np.tile(vrs,(N,1)).T+np.tile(vrs,(N,1))

        Nd = scipy.stats.norm(0, 1)
        prob = Nd.cdf(diff_means/np.sqrt(vars_sum))
        np.fill_diagonal(prob,0)
        
        if self.selective_eig:
            prob_cmp = prob.copy()
            prob_cmp = np.minimum(prob_cmp, 1-prob_cmp)
            prob_cmp = np.tril(prob_cmp)
            prob_cmp = np.divide(prob_cmp.T,np.amax(prob_cmp, 1), out=np.zeros_like(prob_cmp.T), where=np.amax(prob_cmp, 1)!=0).T
        else:
            prob_cmp = np.ones(np.shape(prob))
            
        return prob, prob_cmp

    def get_maximum(self,gain_mat):
        '''
        Function to find the pair of conditions, which, compared would attain a maximum in the information gain matrix. 
        '''
        result = np.where(gain_mat == np.amax(gain_mat))
        result = np.stack((result[0],result[1]), axis=1)
        pair_to_compare = np.expand_dims(result[random.randint(0,np.shape(result)[0]-1),:],0)
        return pair_to_compare

    def compute_information_gain_mat(self,N,G):
        '''
        Given the number of conditions (N) and comparisons performed (G) the function returns the information
        gain matrix and the pair of conditions maximizing the information gain
        '''
        
        kl_divs = np.zeros((N,N))
        
        # Estimate the mean and the standard deviation of the scores obtained from the comparisons collected so far
        Ms_curr, Vs_curr = self.ts_solver.solve(G)
        
        # Save the computed matrices during the intermediate steps of the ts_solver to speed up computations in the 
        # next steps
        self.ts_solver.add_cmps()
        
        # Compute prob to weight entries in the expected information gain matrix and prob_cmps for the selective 
        # expected information gain evaluations
        prob, prob_cmps = self.compute_prob_cmps()
        # Iterate over all possible pairs of conditions

        for ii in range(1,N):
            for jj in range(0,ii):
                # only calculate the the expected information gain for the pairs that are close in the scale (selective evaluations)
                if prob_cmps[ii][jj]>= random.random():

                    if self.approx:
                        Ms, Vs = self.ts_solver.solve_approx(np.array([[0,1]]), Ms_curr[[ii,jj]], Vs_curr[[ii,jj]])
                        kl1 = self.kl_divergence_approx(Ms,Vs,Ms_curr[[ii,jj]],Vs_curr[[ii,jj]])
                    else:
                        Ms, Vs = self.ts_solver.solve(np.vstack((G,np.array([ii,jj]))), num_iters=4, save = False)
                        kl1 = self.kl_divergence_approx(Ms,Vs,Ms_curr,Vs_curr)

                    if self.approx:
                        Ms, Vs = self.ts_solver.solve_approx(np.array([[1,0]]), Ms_curr[[ii,jj]], Vs_curr[[ii,jj]])
                        kl2 = self.kl_divergence_approx(Ms,Vs,Ms_curr[[ii,jj]],Vs_curr[[ii,jj]])
                    else:
                        Ms, Vs  = self.ts_solver.solve(np.vstack((G,np.array([jj,ii]))), num_iters=4, save = False)
                        kl2 = self.kl_divergence_approx(Ms,Vs,Ms_curr,Vs_curr)
                    
                    # Compute expected information gain by weighting the kl divergence by the probability of one condition selected over another
                    kl_divs[ii][jj] = prob[jj][ii]*kl1+prob[ii][jj]*kl2
                else:
                    kl_divs[ii][jj] = -1

        pair_to_compare = self.get_maximum(kl_divs)

        return kl_divs, pair_to_compare

    
    
    def kl_divergence_approx(self,mean_1, var_1, mean_2, var_2):
        '''
        Aproximation of the multivariate normal KL divergence: 
        https://stats.stackexchange.com/questions/60680/kl-divergence-between-two-multivariate-gaussians
        '''
        total = np.sum(np.log(var_2)) - np.sum(np.log(var_1))+sum(var_1/var_2)+np.dot(1/var_2, (mean_1-mean_2)**2)#
        return total


class TrueSkillSolver():
    '''
    Implementation of the TrueSkill:
    http://mlg.eng.cam.ac.uk/teaching/4f13/1920/message%20in%20TrueSkill.pdf
    '''
    
    def __init__(self,N, approx=False):
        self.approx = approx
        self.N = N
        self.Ms = np.zeros(shape=(N))
        self.Vs = 0.5*np.ones(shape = (N))
        self.Mgs = np.empty((0,2))
        self.Pgs = np.empty((0,2))
        self.pv = 0.5
            
        
    def add_cmps(self, numb_cmps=1):
        if numb_cmps>0:
            self.Pgs = np.vstack((self.Pgs, np.zeros((numb_cmps,2))))
            self.Mgs = np.vstack((self.Mgs, np.zeros((numb_cmps,2))))
        else:
            self.Pgs = self.Pgs[:np.shape(self.Pgs)[0]-numb_cmps][:]
            self.Mgs = self.Mgs[:np.shape(self.Mgs)[0]-numb_cmps][:]
        self.Pgs = np.zeros(np.shape(self.Pgs))
        self.Mgs = np.zeros(np.shape(self.Mgs))
        
    def psi(self, x):
        return scipy.stats.norm(0, 1).pdf(x)/scipy.stats.norm(0, 1).cdf(x)

    def lamb(self, x):
        ps = self.psi(x)
        return ps*(ps + x)


    def solve_approx(self, G, mv, pv):

        # https://www.microsoft.com/en-us/research/project/trueskill-ranking-system/
        N = np.shape(G)[0]
        for ii in range(0,N):
            c = np.sqrt(2+pv[G[ii,0]]+pv[G[ii,1]])

            term = (mv[G[ii,0]]-mv[G[ii,1]])/c
            #normpdf_term = scipy.stats.norm(0, 1).pdf(term)
            #normcdf_term = scipy.stats.norm(0, 1).cdf(term)
            #print(normcdf_term)
            fact_1 =self.psi(term)
            fact_2 =self.lamb(term)
            mv[G[ii,0]] = mv[G[ii,0]]+pv[G[ii,0]]/c*fact_1
            mv[G[ii,1]] = mv[G[ii,1]]-pv[G[ii,1]]/c*fact_1

            pv[G[ii,0]] = pv[G[ii,0]]*(1-pv[G[ii,0]]*fact_2/c**2)
            pv[G[ii,1]] = pv[G[ii,1]]*(1-pv[G[ii,1]]*fact_2/c**2)

        return mv, pv
    
    def solve(self, G, num_iters = 6, save = True):

        if np.shape(self.Mgs)[0]!=np.shape(G)[0]:
            self.add_cmps(np.shape(G)[0]-np.shape(self.Mgs)[0])

        Ps = 1./self.Vs
        Ms = self.Ms.copy()
        for ii in range(0, num_iters):
            Psg = np.reshape(Ps[G],np.shape(G)) - self.Pgs
            Msg = (np.reshape(Ps[G]*Ms[G],np.shape(G)) - self.Pgs*self.Mgs)/Psg

            vgt = 1+np.sum(1/Psg,1)
            mgt = Msg[:,0] - Msg[:,1]

            Mt = mgt + np.sqrt(vgt)*self.psi(mgt/np.sqrt(vgt))
            Pt = 1/( vgt*( 1-self.lamb(mgt/np.sqrt(vgt)) ) )

            ptg = Pt -1/vgt
            mtg = (Mt*Pt - mgt/vgt)/ptg

            self.Pgs = 1/(1+np.tile(1/ptg,(2,1)).T+1/Psg[:,[1,0]])
            self.Mgs = np.stack((mtg, -mtg),axis=-1) + Msg[:,[1,0]]

            for p in range(0,self.N):
                Ps[p] = 1/self.pv+np.sum(self.Pgs[G==p])
                Ms[p] = np.sum(self.Pgs[G==p]*self.Mgs[G==p])/Ps[p]

        if save:
            self.Vs = 1/Ps
            self.Ms = Ms.copy()

        return Ms, 1/Ps
