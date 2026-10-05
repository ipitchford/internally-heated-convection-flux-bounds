"""Local slope d ln Nu / d ln Gamma along a fixed-Ra roll family, and the 'shoulder' (least negative slope beyond
the optimum), located by a quadratic fit through the three nearest slopes."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from paths import ROOT
import numpy as np, json
def slopes(path):
    d = sorted(json.load(open(path)), key=lambda x: x['Gamma'])
    G = np.array([x['Gamma'] for x in d]); N = np.array([x['Nu'] for x in d])
    k = np.unique(G, return_index=True)[1]; G, N = G[k], N[k]
    return np.sqrt(G[1:]*G[:-1]), np.diff(np.log(N))/np.diff(np.log(G))
def shoulder(path):
    d = sorted(json.load(open(path)), key=lambda x: x['Gamma'])
    G = np.array([x['Gamma'] for x in d]); N = np.array([x['Nu'] for x in d]); Gopt = G[np.argmax(N)]
    Gm, s = slopes(path); lg = np.log(Gm)
    i0 = int(np.argmax(np.where(lg > np.log(Gopt)+0.15, s, -9.0)))
    c = np.polyfit(lg[i0-1:i0+2], s[i0-1:i0+2], 2); x = -c[1]/(2*c[0])
    return dict(Gamma_opt=float(Gopt), Nu_opt=float(N.max()), Gamma_shoulder=float(np.exp(x)), slope_shoulder=float(np.polyval(c, x)),
                Gmin=float(G.min()), Gmax=float(G.max()), max_slope_beyond_opt=float(s[lg > np.log(Gopt)+0.01].max()))
if __name__ == '__main__':
    for Ra in ['1e6', '1e7']:
        r = shoulder(f'{ROOT}/data/steady/gscan_rb_{Ra}.json'); r['Gamma_shoulder_Ra^0.2'] = r['Gamma_shoulder']*float(Ra)**0.2
        print(Ra, json.dumps(r))

def local_maxima(path):
    """All interior local maxima of Nu(Gamma) on the scan, refined by a parabola in log Gamma through the
    three points around each discrete maximum."""
    d = sorted(json.load(open(path)), key=lambda x: x['Gamma'])
    G = np.array([x['Gamma'] for x in d]); N = np.array([x['Nu'] for x in d])
    k = np.unique(G, return_index=True)[1]; G, N = G[k], N[k]; out = []
    for i in range(1, len(G)-1):
        if N[i] >= N[i-1] and N[i] >= N[i+1]:
            c = np.polyfit(np.log(G[i-1:i+2]), N[i-1:i+2], 2); x = -c[1]/(2*c[0])
            out.append(dict(Gamma=float(np.exp(x)), Nu=float(np.polyval(c, x))))
    return out
