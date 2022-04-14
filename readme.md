# ASAP: Active Sampling for Pairwise Comparisons
A toolkit for actively pairing conditions in pairwise comparison preference aggregation.

Pairwise comparison data arise in many domains with subjective assessment experiments. In these experiments participants are asked to express a preference between two conditions. However, many pairwise comparison protocols require a large number of comparisons to infer accurate scores, which may be unfeasible when each comparison is time-consuming or expensive. To address this problem we propose ASAP, an active sampling algorithm, offering the highest accuracy of inferred scores compared to the existing methods. Unlike most existing methods, which rely on partial updates of the posterior distribution, we are able to perform full updates and therefore much improve the accuracy of the inferred scores.

If using the code please cite:

A. Mikhailiuk, C. Wilmot, M. Perez-Ortiz, D. Yue and R. K. Mantiuk, 2020. “[Active Sampling for Pairwise Comparisons via Approximate Message Passing and Information Gain Maximization](https://arxiv.org/abs/2004.05691)”, International Conference on Pattern Recognition (ICPR)

```
@INPROCEEDINGS{mikhailiuk2020active,
    title={Active Sampling for Pairwise Comparisons via Approximate Message Passing and Information Gain Maximization},
    author={Aliaksei Mikhailiuk and Clifford Wilmot and Maria Perez-Ortiz and Dingcheng Yue and Rafal Mantiuk},
    booktitle={2020 IEEE International Conference on Pattern Recognition (ICPR)}, 
    year={2021},
    month={Jan},
}
```


The active sampling package includes:

* Two versions of the algorithm: (i) ASAP - an accurate version with full posterior update; (ii) ASAP-approx. - version with an online posterior update for reduced computation cost;
* ASAP and ASAP-approx. both in python and matlab;
* ASAP and ASAP-approx. with support for batch and sequential modes
* Pytorch-GPU implementation of ASAP for large scale experiments;
* Matlab implementation of 9 commonly used sampling strategies.

### Video explanation

Click on the image below to be re-directed to youtube.

[![IMAGE ALT TEXT](http://img.youtube.com/vi/Yt3-3zOR9u4/0.jpg)](http://www.youtube.com/watch?v=Yt3-3zOR9u4 "Video Title")


### Usage

ASAP and ASAP-approx return either a pair or a batch of pairs with the highest expected information gain. By default batch mode is used.

Both python and matlab code use run_asap(M) function to produce pairs to be compared. Where M is a pairwise comparison matrix (NxN), with each element M[ii][jj] - the number of times condition ii was selected over condition jj in comparisons performed so far.

run_asap(M) function returns an either a ((N-1)x2) array of indeces ii jj of elements to be compared, or a single pair maximising information gain. Note that the output might be different from the one given below due to the use of selective expected information gain evaluation.

Python:

```
# import asap from the file
import asap_cpu
import numpy as np

# Matrix with pairwise comparisons
M = np.array([[0,1,2,3,1],
              [1,0,2,3,1],
              [1,2,0,3,1],
              [1,2,3,0,1],
              [1,2,3,1,0]])
              
N = np.shape(M)[0]

# Create an object of class passing the number of conditions
asap = asap_cpu.ASAP(N)

# Run active sampling algorithm on the matrix of comaprisons
pairs_to_compare = asap.run_asap(M)

# Calling print
print(pairs_to_compare)

Output: 

[[0 4]
 [2 4]
 [0 1]
 [0 3]]

```

Matlab:


```
% Create matrix with pairwise comparisons
M = [0,1,2,3,1;
     1,0,2,3,1;
     1,2,0,3,1;
     1,2,3,0,1;
     1,2,3,1,0];
              
% Run active sampling
pairs_to_compare = run_asap(M, 'mst');

% Calling print
display(pairs_to_compare)

Output: 

pairs_to_compare =

     1     5
     1     2
     1     3
     1     4

```


## Literature

We also include the code for the compared in the paper methods.

Swiss-system: L. Csato, 2017. Ranking in Swiss system chess team tournaments. Annals of Operations Research. 254. 10.1007/s10479-017-2440-4. 

Quicksort: L.  Maystre  and  M.  Grossglauser, 2017.  “Just  sort  it!  A  simple  and  effective approach  to  active  preference  learning”,  in Proceedings  of  the 34th International Conference on Machine Learning, vol. 70. pp. 2344–2353.

Adaptive-rectangular-design (ARD): ITU-R Recomendations, 2016.  “Subjective  assessment  methods  for  3D  video  quality”,  P.915.

Matchmaking: R.  Herbrich,  T.  Minka,  and  T.  Graepel, 2019.  “Trueskill TM:  A  bayesian  skill rating system”, in Advances in Neural Information Processing Systems

YeDoermann (computationally very slow, included for completeness): P. Ye and D. Doermann, 2014. “Active sampling for subjective image quality assessment”, IEEE Conference on Computer Vision and Pattern Recognition (CVPR), pp. 4249–4256.

API (Numerically unstable, included for completeness): T. Pfeiffer, X. Gao, Y. Chen, A. Mao, and D. Rand, 2012. “Adaptive polling for information aggregation”, AAAI Conference on Artificial Intelligence.

In the paper we used official implemenation from:

HR-Active: Q.  Xu,  J.  Xiong,  X.  Chen,  Q.  Huang,  and  Y.  Yao, 2018  “Hodgerank  with information  maximization  for  crowdsourced  pairwise  ranking  aggregation”, AAAI  Conference  on  Artificial Intelligence, pp. 4326–4334. 

AKG: X.  Chen,  K.  Jiao,  and  Q.  Lin, 2016.  “Bayesian  decision  process  for  cost-efficient  dynamic  ranking  via  crowdsourcing”, Journal  of  MachineLearning Research, vol. 17, no. 216, pp. 1–40.

Crowd-BT: X. Chen, P. N. Bennett, K. Collins-Thompson, and E. Horvitz, 2013. “Pairwise ranking aggregation in a crowdsourced setting”, Proceedings of the Sixth ACM  International  Conference  on  Web  Search  and  Data  Mining,  pp.193–202

Hybrid-MST: J.  Li,  R.  Mantiuk,  J.  Wang,  S.  Ling,  and  P.  Le  Callet, 2018.  “Hybrid-mst: A hybrid active sampling strategy for pairwise preference aggregation,” NIPS, 31st Conference on Neural Information Processing Systems. [LINK](https://github.com/jingnantes/hybrid-mst)



## License

MIT
