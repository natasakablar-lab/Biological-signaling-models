"""Meaningful UI regression checks. Run: python test_streamlit_app.py."""
from pathlib import Path
from streamlit.testing.v1 import AppTest

APP=Path(__file__).resolve().parent/'app.py'
at=AppTest.from_file(str(APP),default_timeout=120).run()
assert not at.exception, at.exception
assert at.metric[0].value=='Ustaljene oscilacije'
assert abs(at.session_state['simulation']['metrics']['period_s']-1291.04)<.1
buttons=lambda label:next(b for b in at.button if b.label==label)
buttons('Izračunaj Hopfove pragove').click().run()
assert not at.exception,at.exception
assert abs(at.session_state['hopf']['points'][0]['V1']-.386227346)<1e-6
buttons('Izračunaj prvi Ljapunovljev koeficijent').click().run()
assert not at.exception,at.exception
assert len(at.session_state['lyapunov'])==2
assert all(row['l1']<0 for row in at.session_state['lyapunov'])
# A preset updates pending widgets but does not relabel previous computed data.
at.selectbox(key='preset').select('Fosfataza KK · 0.5 nM').run()
assert at.session_state['simulation']['parameters']['phosphatase_KK']==1.
buttons('Primeni i simuliraj').click().run()
assert not at.exception,at.exception
assert at.metric[0].value=='Dostignuta ravnoteža'
assert 'hopf' not in at.session_state and 'lyapunov' not in at.session_state
assert at.session_state['simulation']['parameters']['phosphatase_KK']==.5
assert at.session_state['simulation']['metrics']['KPP_half_amplitude_nM']==0.
# Recompute thresholds for changed parameters, not cached baseline values.
buttons('Izračunaj Hopfove pragove').click().run()
assert not at.exception,at.exception
assert abs(at.session_state['hopf']['points'][1]['V1']-4.33302)<1e-4
at.selectbox(key='preset').select('Referentni model · ravnoteža').run()
buttons('Primeni i simuliraj').click().run()
assert not at.exception,at.exception
assert abs(at.session_state['simulation']['equilibrium'][7]-2.573760)<1e-5
print('PASS: initial oscillations, Hopf, Lyapunov, pending preset, equilibrium, cache invalidation, reference equilibrium')
