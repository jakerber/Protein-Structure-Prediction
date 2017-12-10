"""
CS86 Final Project – Main Script

Joshua Kerber
06/01/2107
Compare protein structure prediction methods:
    Fragment substitution
        vs.
    Perturbing Phi Psi angles
"""

import threading
from ProteinPredictSimulation import ProteinPredictSimulation

# methods used in simulation
METHODS = [
    'fragment_substitution',
    'perturb_phi_psi'
]

# proteins we will use in simulation, mapped to respective sequences
PROTEINS = {
    '4hgu':'EAAVCTTEWDPVCGKDGKTYSNLCWLNEAGVGLDHEGECL',
    '4xdx':'KELRCQCIKTYSKPFHPKFIKELRVIESGPHCANTEIIVKLSDGRELCLDPKENWVQRVVEKFLKRAENS',
    '5e4x':'YAVAESVIGKRVGDDGKTIEYLVKWTDMSDATWEPQDNVDSTLVLLYQQQ',
    '5eca':'SHMTFVALYDYVASGETDLSFKKGERLQIVNNTEGDWWLAHSLTTGRTGYIPSNYVAPSD',
    '5eh6':'EPEITLIIFGVIAGVIGTILLISYGIRRLC'
}

def runSim(protein, method):
    """
    Run structure prediction simulaton using ProteinPredictSimulation object
    :param protein: protein used in simulation
    :param method: folding method used
    """
    sequence = PROTEINS[protein]
    newSim = ProteinPredictSimulation(protein, sequence, method)
    newSim.run()


if __name__ == "__main__":
    # run simulation in own thread for all protein and method combinations
    for protein in PROTEINS:
        for method in METHODS:
            t = threading.Thread(target=runSim, args=[protein, method])
            t.start()
