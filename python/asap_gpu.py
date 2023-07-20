import torch
import torch.distributions as dist
import numpy as np
import scipy.sparse as ssp
import networkx as nx
import random
torch.set_printoptions(precision=3,sci_mode=False)
np.set_printoptions(precision=6, suppress=True)

def Psi(x):
    normal = dist.Normal(
        torch.tensor(0.0, device=x.device),
        torch.tensor(1.0, device=x.device))
    return normal.log_prob(x).exp() / normal.cdf(x)

def Phi(x):
    psi = Psi(x)
    return psi*(psi+x)

# Figure 1 in the paper
def true_skill(G: torch.Tensor, M: int, num_iters=4):
    '''
    Implementation of the TrueSkill from http://mlg.eng.cam.ac.uk/teaching/4f13/1920/message%20in%20TrueSkill.pdf
    '''
    _, *dim, N = G.size()
    I, J = 0, 1
    idx, count = G[:2,...], G[2,...].float()
    pv = 0.5

    Ms = torch.zeros(*dim, M).to(G.device)
    Ps = torch.full((*dim, M), 1/pv).to(G.device)

    Mgs = torch.zeros(2, *dim, N).to(G.device)
    Pgs = torch.zeros(2, *dim, N).to(G.device)

    Msg = torch.zeros(2, *dim, N).to(G.device)
    Psg = torch.zeros(2, *dim, N).to(G.device)

    n2m = lambda An, value=0: (torch.full((*dim, M), value,dtype =torch.float)
                               .to(An.device)
                               .scatter_add_(-1, idx[I], An[I])
                               .scatter_add_(-1, idx[J], An[J]))
    m2n = lambda Am: torch.stack((
        Am.gather(-1, idx[I]),
        Am.gather(-1, idx[J])))

    for i in range(num_iters):
        # step2 compute skill to game message
        Psg = m2n(Ps) - Pgs
        Msg = (m2n(Ps*Ms) - Pgs*Mgs)/Psg
        # step3 compute to performance message
        vgt = 1 + 1/Psg[I] + 1/Psg[J]
        mgt = Msg[I] - Msg[J]
        # step4 compute marginal performance
        Mt = mgt + vgt.sqrt() * Psi(mgt/vgt.sqrt())
        Pt = 1/(vgt * (1 - Phi(mgt/vgt.sqrt())))
        # step5 compute performance to game message
        ptg = Pt - 1/vgt
        mtg = (Mt*Pt - mgt/vgt)/ptg
        # step6 compute game to skills messages
        Pgs[I] = 1/(1 + 1/ptg + 1/Psg[J])
        Pgs[J] = 1/(1 + 1/ptg + 1/Psg[I])
        Mgs[I] = Msg[J] + mtg
        Mgs[J] = Msg[I] - mtg
        # step1 compute marginal variance
        Ps = n2m(Pgs*count, 1/pv)
        Ms = n2m(Pgs*Mgs*count) / Ps
    return dist.Normal(Ms, 1/Ps.sqrt()), Ms, 1/Ps.sqrt()


def prob_cmps(normal: dist.Normal):
    '''
    prob: matrix with probability of one condition chosen over another with 
    prob[ii][jj] computed from 
        
    '''
    m, v = normal.mean, normal.variance
    mi, mj = torch.meshgrid((m, m))
    vi, vj = torch.meshgrid((v, v))
    N = dist.Normal(0, 1)
    prob = N.cdf((mi - mj)/(1+vi+vj).sqrt()).fill_diagonal_(0)

    return prob


def compute_minimum_spanning_tree(inf_mat):
    '''
    Given an information gain matrix, we want to extract a set of comparisons to perform that
    would have the largest total information gain and would form a connected graph of comparisons.
    The function takes as input an information gain matrix, then computes its reciprocal (1/inf_gain)
    and extracts a minimum spanning tree from it.
    '''

    inf_mat[inf_mat<=0.0] = np.inf
    inf_mat = 1/inf_mat
    
    GrMST=nx.from_numpy_array(inf_mat)
    T=nx.minimum_spanning_tree(GrMST)

    pairs_to_compare = np.asarray(T.edges())
    
    edges=sorted(T.edges(data=True), key=lambda t: t[2].get('weight', 1))

    pairs_to_compare = np.array([t[0:2] for t in edges ])
    return pairs_to_compare


def kl_divergence_approx(mean_1, std_1, mean_2, std_2):
    normal0 = dist.Normal(mean_1, std_1)
    normal = dist.Normal(mean_2, std_2)
    total = (1 + 2*dist.kl.kl_divergence(normal0, normal)).sum()
    return total


def get_maximum(gain_mat):
    '''
    Function to find the pair of conditions, which, compared would attain a maximum in the information gain matrix. 
    '''
    result = np.where(gain_mat == np.amax(gain_mat))
    result = np.stack((result[0],result[1]), axis=1)
    pair_to_compare = np.expand_dims(result[random.randint(0,np.shape(result)[0]-1),:],0)
    return pair_to_compare

def ASAP(cmp_matrix: np.ndarray, mst_mode=True, cuda=False, get_scores=False):
    '''

    Function to compute the next batch of comparisons to perform. 

    Note: in current implementation selective EIG evaluations are not implemented.
    This does not impact the accuracy of the algorithm.

    '''

    # Create a sparse matrix
    cmp_matrix = ssp.coo_matrix(cmp_matrix)

    # Check that the matrix is square
    M, M2 = cmp_matrix.shape
    assert M == M2, "The comparsion must be square matrix"

    # convert the matrix to tensor to have G(rows: condition1, condition2, comparison_outcomes, cols: pairing_combination) 
    G0 = torch.stack((
        torch.tensor(cmp_matrix.row).long(),
        torch.tensor(cmp_matrix.col).long(),
        torch.tensor(cmp_matrix.data).long(),
    ))

    # Put on GPU if available
    if cuda:
        G0 = G0.cuda()

    # Compute the scores
    normal0, Ms0, StD0 = true_skill(G0, M)


    # tensor to hold each possible pairwise outcome in the next comparison
    I, J = torch.nonzero((1 - torch.eye(M)).to(G0.device), as_tuple=False).unbind(-1)
    G = torch.zeros(G0.size(0),I.size(0),G0.size(1)+1).to(G0)
    G[:,:,:-1] = G0.unsqueeze(-2)
    G[0,:,-1] = I
    G[1,:,-1] = J
    G[2,:,-1] = 1

    # Compute score distributions for each of the possible outcome combinations
    normal,_,_ = true_skill(G, M, num_iters=4)

    
    # Compute expected information gain for all posible outcomes
    kl_divs = prob_cmps(normal0)
    for cc in range(0,len(I)):
        kl_divs[I[cc],J[cc]] *= kl_divergence_approx(normal0.loc, normal0.scale,
                                                     normal.loc[cc], normal.scale[cc])

    info_gain = kl_divs + kl_divs.T
    info_gain = info_gain.cpu().detach().numpy()
    pairs_to_compare = get_maximum(info_gain)
    if mst_mode:    
        # minial spanning tree for batch mode
        pairs_to_compare = compute_minimum_spanning_tree(info_gain)
    
    
    if get_scores:
        return pairs_to_compare, Ms0.cpu().detach().numpy(), StD0.cpu().detach().numpy()
        
    return pairs_to_compare
