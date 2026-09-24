"""Precursor/signal/mature segment ablation for P0ACW4, post hoc hypothesis check."""
import json,pathlib
import numpy as np, torch
from pepdesign.encoding import combined_features
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import predict_scores,predict_scores_shared_adj
ROOT=pathlib.Path(__file__).resolve().parents[1]

def score(seq,cnn,gnn):
 def feat(ss,n): return (torch.from_numpy(np.stack([combined_features(s,n) for s in ss])),)
 cs=float(predict_scores(cnn,feat,[seq],150)[0]);gs=float(predict_scores_shared_adj(gnn,[seq],150)[0])
 return {'length':len(seq),'cnn':cs,'gnn':gs,'ensemble':(cs+gs)/2}

def main():
 c=PepCNN();c.load_state_dict(torch.load(ROOT/'results/study07_cnn.pt',weights_only=False)['state_dict']);c.eval()
 g=PepGNN();g.load_state_dict(torch.load(ROOT/'results/study07_gnn.pt',weights_only=False)['state_dict']);g.eval()
 s=json.loads((ROOT/'results/study12_candidates.json').read_text())['candidates'][0]['sequence']
 cases={'P0ACW4_precursor':s,'P0ACW4_mature_after20':s[20:],'P0ACW4_signal_1_20':s[:20]}
 out={name: {'sequence':seq,**score(seq,c,g)} for name,seq in cases.items()}
 out['source']='https://ecoliwiki.org/colipedia/index.php/ydcA:Gene_Product(s) (signal peptide residues 1-20)'
 (ROOT/'results/study15_signal_ablation.json').write_text(json.dumps(out,indent=2))
 print(json.dumps(out,indent=2))
if __name__=='__main__': main()
