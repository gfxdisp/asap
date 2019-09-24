#include "mex.h"
//compute AUC for ranking
// [auc, total_cnt]=calc_auc(true_label, pred_label)
#define sign(X) (((X) > 0) ? 1 : (((X) < 0) ? -1 : 0))

void CheckInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Check for proper number of arguments. */
  if (nrhs < 2)
     mexErrMsgTxt("calc_auc(true_label, pred_label).\n");
  
  if(mxIsSparse(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[1]) || mxIsComplex(prhs[1]))
    mexErrMsgTxt("Input must be a dense real vector\n");
  
  if ((mxGetM(prhs[0])*mxGetN(prhs[0]))!= (mxGetM(prhs[1])*mxGetN(prhs[1])))
    mexErrMsgTxt("The lengths of true label and test label are not the same\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
  
  double auc, *t_label, *pred_label;
  size_t n,m;
  int i,j,cnt=0,total_cnt=0;
  
  CheckInput(nlhs, plhs, nrhs, prhs);
  
  n=mxGetM(prhs[0])*mxGetN(prhs[0]);
  
  t_label=mxGetPr(prhs[0]);
  
  pred_label=mxGetPr(prhs[1]);
  
  for (i=0; i<n-1; ++i)
      for (j=i+1; j<n; ++j)
          if (t_label[i]==t_label[j])
            continue;
          else{
            total_cnt++;
            cnt+=sign(t_label[i]-t_label[j])==sign(pred_label[i]-pred_label[j])? 1 : 0;
          }
  auc=(double)cnt/(double) total_cnt;
  plhs[0] = mxCreateDoubleScalar(auc);
  plhs[1] = mxCreateDoubleScalar(total_cnt);
}
