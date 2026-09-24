"""Post-hoc stratified audit of the confounded UniProt scale-up test set."""
import csv,json,pathlib
import numpy as np,torch
from sklearn.metrics import roc_auc_score,average_precision_score
from pepdesign.encoding import combined_features
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import predict_scores,predict_scores_shared_adj
ROOT=pathlib.Path(__file__).resolve().parents[1]
def cnn_inputs(ss,n): return (torch.from_numpy(np.stack([combined_features(s,n) for s in ss])),)
def main():
 rows=list(csv.DictReader(open(ROOT/'data/processed_large/test.csv')))
 seqs=[r['sequence'] for r in rows]; y=np.array([int(r['label']) for r in rows]);
 cnn=PepCNN(); cnn.load_state_dict(torch.load(ROOT/'results/study07_cnn.pt',weights_only=False)['state_dict']);cnn.eval()
 gnn=PepGNN(); gnn.load_state_dict(torch.load(ROOT/'results/study07_gnn.pt',weights_only=False)['state_dict']);gnn.eval()
 scores={'PepCNN':predict_scores(cnn,cnn_inputs,seqs,150), 'PepGNN':predict_scores_shared_adj(gnn,seqs,150)}
 strata={'all':np.ones(len(rows),dtype=bool),
         'reviewed_only':np.array([r['header'].startswith('sp|') for r in rows]),
         'unreviewed_positive_vs_reviewed_negative':np.array([r['header'].startswith('tr|') or r['label']=='0' for r in rows])}
 out={}
 for name,mask in strata.items():
  yy=y[mask]; rec={'n':int(mask.sum()),'n_pos':int(yy.sum()),'n_neg':int((yy==0).sum())}
  for model,s in scores.items():
   ss=s[mask]; rec[model]={'auroc':roc_auc_score(yy,ss),'auprc':average_precision_score(yy,ss),
                           'mean_positive_score':float(ss[yy==1].mean()),'mean_negative_score':float(ss[yy==0].mean())}
  out[name]=rec
 # Any difference is descriptive because positive strata differ in families and size.
 (ROOT/'results/study17_review_strata.json').write_text(json.dumps(out,indent=2))
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
