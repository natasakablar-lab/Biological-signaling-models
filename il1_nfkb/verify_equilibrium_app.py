"""Functional tests for equilibrium/branch actions and boundary-case messages."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from streamlit.testing.v1 import AppTest

def button(at,key): return next(b for b in at.button if b.key==key)

def main():
    at=AppTest.from_file(str(Path(__file__).parent/'app.py'),default_timeout=90).run()
    assert not at.exception,[e.message for e in at.exception]
    button(at,'equilibrium_button').click().run()
    assert not at.exception,[e.message for e in at.exception]
    assert at.session_state.equilibrium_result[2]['stability']=='stabilna'
    branch_button=next(b for b in at.button if b.label=='Prati ravnotežnu granu')
    branch_button.click().run()
    assert not at.exception,[e.message for e in at.exception]
    branch=at.session_state.branch_result[1]
    assert len(branch['rows'])==81 and not branch['hopf_candidates']
    selector=next(s for s in at.selectbox if s.label=='Parametar za analizu')
    selector.select('k19');next(b for b in at.button if b.label=='Prati ravnotežnu granu').click().run()
    assert not at.exception,[e.message for e in at.exception]
    assert at.session_state.branch_result[1]['parameter']=='k19'
    at.sidebar.selectbox[0].select('Bez IL-1, sa početnim aktivnim IKK')
    next(b for b in at.button if b.label=='Izračunaj odzive').click().run()
    button(at,'equilibrium_button').click().run()
    assert not at.exception,[e.message for e in at.exception]
    assert any('bilansnih' in e.value for e in at.error)
    out=Path(__file__).parent/'results_equilibria'
    (out/'ui_verification_equilibria.json').write_text(json.dumps({'passed':['equilibrium action','81-point branch action','parameter change to k19','unsupported zero IL1 total clearly reported']},indent=2))
    print('Equilibrium UI checks passed')

if __name__=='__main__': main()
