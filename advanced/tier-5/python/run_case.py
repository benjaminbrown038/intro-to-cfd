"""Run Tier 5. See the tier README for model scope."""
from pathlib import Path
import sys,json,argparse
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE.parents[1]/'common'/'python'))
from flow_examples import run
if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--case',type=Path,default=BASE/'cases'/'adaptive.json')
    p.add_argument('--output',type=Path,default=BASE/'results'/'python')
    args=p.parse_args()
    run(5,json.loads(args.case.read_text()),args.output)
