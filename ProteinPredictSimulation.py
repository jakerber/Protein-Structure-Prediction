"""
CS86 Final Project – Simulation Object

Joshua Kerber
06/01/2107
Predict a protein structure from scratch using a pyRosetta simulation
    Monte Carlo method used
    Root-mean-square-deviation (RMSD) data output written to file
"""

from rosetta import *
import random
import math
import datetime

def getCurrentDateTime():
    """
    Get current formatted datetime stamp
    """
    now = datetime.datetime.now()
    return '{0}{1}{2}_{3}{4}{5}'.format(
        now.month,
        now.day,
        now.year,
        now.hour,
        now.minute,
        now.second
    )

class ProteinPredictSimulation:
    """
    Monte Carlo simulation to predict the structure of a protein
    """
    def __init__(self, protein, sequence, method, steps=100000):
        self.protein = protein
        self.sequence = sequence
        self.method = method
        self.steps = steps

        # simulation prep
        self.counter = 0
        self.initSim = False
        self.bestScore = float('inf')
        self.datetime = getCurrentDateTime()
        self.filename = 'output/{0}/{1}/output_{2}.csv'.format(
            self.method,
            self.protein,
            self.datetime)
        self.ktVal = None
        self.fragSet = None
        self.curScore = None
        self.scoreFunc = None

        # unqiue folding method functions
        self.methodFuncs = {
            'fragment_substitution': lambda pose: self.fragmentSub(pose),
            'perturb_phi_psi': lambda pose: self.perturbPhiPsi(pose),
        }

        # poses used in sim
        self.testPose = None
        self.correctPose = None
        self.curPose = None
        self.optimalPose = None

    def perturbPhiPsi(self, pose):
        """
        Perturb phi psi angles step
        :param pose: current protein pose
        """
        # pick random residue number between 1 and 12
        n = random.randrange(1, 13)

        # perturbate phi and psi angles of pose
        pose.set_phi(n, self.testPose.phi(n) + random.randrange(-10, 11))
        pose.set_psi(n, self.testPose.psi(n) + random.randrange(-10, 11))

    def fragmentSub(self, pose):
        """
        Fragment substitution step
        :param pose: current protein pose
        """
        # create mover that will swap in a random fragment
        movemap = MoveMap()
        movemap.set_bb(True)
        mover = ClassicFragmentMover(self.fragSet, movemap)

        # make a move
        mover.apply(pose)

    def init(self):
        """
        Initialize all simulation variables to ready state
        """
        init()  # pyrosetta
        self.ktVal = 1
        self.counter = 0

        # open file used for output
        self.filePointer = open(self.filename, 'w')
        self.filePointer.write('RMSD,score\n')

        # generate protein poses
        self.curPose = pose_from_sequence(self.sequence, 'centroid')
        self.correctPose = pose_from_pdb('proteins/{0}/{0}_clean.pdb'.format(self.protein))
        self.optimalPose = Pose()
        self.testPose = Pose()
        self.testPose.assign(self.curPose)

        # generate fragment library
        self.fragSet = ConstantLengthFragSet(3)
        self.fragSet.read_fragment_file('proteins/{0}/{0}-3frag.txt'.format(self.protein))

        # create scoring function
        self.scoreFunc = create_score_function('score3')
        self.curScore = self.scoreFunc(self.curPose)

    def run(self):
        """
        Initialize and run the protein folding simulation
        """
        self.init()
        for _ in range(self.steps):

            # make move using unique method step
            self.methodFuncs[self.method](self.testPose)
            testScore = self.scoreFunc(self.testPose)

            # accept move based on random chance algorithm (Monte Carlo)
            randomChance = math.exp(-(testScore - self.curScore) / self.ktVal)
            if testScore < self.curScore or random.random() < randomChance:
                if testScore < self.bestScore:
                    self.bestScore = testScore

                    # save pose if best score
                    self.optimalPose.assign(self.testPose)

                # decrease kT energy value
                self.ktVal /= 1.05
                self.curScore = testScore
                self.curPose.assign(self.testPose)
                self.counter += 1

                # write RMSD data to file every 100 accepted moves
                if self.counter == 100:
                    rmsd = CA_rmsd(self.correctPose, self.curPose)
                    output = '{},{}\n'.format(rmsd, str(self.curScore))
                    self.filePointer.write(output)
                    self.counter = 0

            # move not accepted, increase kT value
            else:
                self.testPose.assign(self.curPose)
                self.ktVal *= 1.01

        # write final RMSD data to file
        rmsd = CA_rmsd(self.correctPose, self.optimalPose)
        output = '{},{}\n'.format(rmsd, self.scoreFunc(self.optimalPose))
        self.filePointer.write(output)

        # save best pose into new pdb file
        dump_pdb(self.optimalPose, 'output/{0}/{1}/pose_{2}.pdb'.format(
            self.method,
            self.protein,
            self.datetime))
