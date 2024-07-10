'''
Script for integrative modeling of Notch1-DLL4 complex
using a topology file
'''
# Imports
from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.crosslinking
import IMP.pmi.dof
import IMP.atom

import sys

num_frames=20000 #TWEAK 1

rb_max_trans = float(sys.argv[1])

bead_max_trans = float(sys.argv[2])

rex_max_temp = float(sys.argv[3]) 
# Topology File
topology_file = "/home/arastu/notch/dll4/topology_notch1_dll4_rigid.txt"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# All IMP systems start out with a Model
mdl = IMP.Model()

# Read the topology file for a given state
t = IMP.pmi.topology.TopologyReader(topology_file)

# Create a BuildSystem macro to and add a state from a topology file
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(t)

# executing the macro will return the root hierarchy and degrees of freedom (dof) objects
root_hier, dof = bs.execute_macro(max_rb_trans=rb_max_trans,
                                  max_rb_rot=0.1,
                                  max_bead_trans=bead_max_trans)

# TWEAK 2 the dof bead max trans and rigid body max trans

# It's useful to have a list of the molecules.
molecules = t.get_components()

#############################################################################
##################### Fix the part that is transmembrane ####################
#############################################################################
notch1_tm_bead = IMP.atom.Selection(
        root_hier, molecule="notch1", residue_indexes=range(1736,1757)).get_selected_particles()[0]

IMP.core.XYZ(notch1_tm_bead).set_coordinates(IMP.algebra.Vector3D(0,0,0))

dll4_tm_bead = IMP.atom.Selection(
        root_hier, molecule="dll4", residue_indexes=range(530,551)).get_selected_particles()[0]

IMP.core.XYZ(dll4_tm_bead).set_coordinates(IMP.algebra.Vector3D(0,1200,0))

# First select and gather all particles to fix.
fixed_particles = []
fixed_particles.append(notch1_tm_bead)

fixed_particles.append(dll4_tm_bead)
print(fixed_particles)

# Fix the Corresponding Rigid movers and Super Rigid Body movers using dof
# The flexible beads will still be flexible (fixed_beads is an empty list)!
fixed_beads, fixed_rbs = dof.disable_movers(fixed_particles,
                                            [IMP.core.RigidBodyMover,
                                             IMP.pmi.TransformMover,IMP.core.BallMover])

#####################################################
##################### RESTRAINTS ####################
#####################################################

# Restraints define functions that score the model based on
# input information.
#
# Restraint objects are first created in the definition.
# To be evaluated, the restraint object must be add_to_model().
#
# In some cases, sampled parameters for restraints must be added to the DOF
# object

# The output_objects list is used to collect all restraints
# where we want to log the output in the STAT file.
# Each restraint should be appended to this list.
output_objects = []

# -----------------------------
# %%%%% CONNECTIVITY RESTRAINT
#
# Restrains residues/particles that are collected in sequence
# This should be used for any system without an atomic force field (e.g. CHARMM)
# We apply the restraint to each molecule

for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)

# -----------------------------
# %%%%% EXCLUDED VOLUME RESTRAINT
#
# Keeps particles from occupying the same area in space.
# Here, we pass a list of both molecule chains to included_objects to apply this to every residue.
# We could also have passed root_hier to obtain the same behavior.
#
# resolution=1000 applies this expensive restraint to the lowest resolution for each particle.
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                            included_objects=[root_hier],
                                            resolution=1000)
output_objects.append(evr)

#Added this restraint to keep notch from bending
dsr5 = IMP.pmi.restraints.basic.DistanceRestraint(tuple_selection1=(20,20,"notch1",0),tuple_selection2=(1735,1735,"notch1",0),distancemax=865,distancemin=845,root_hier=root_hier,label='ntcstr')
dsr5.add_to_model()
output_objects.append(dsr5)

#####################################################
###################### SAMPLING #####################
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.

# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=100,
                                    excluded_rigid_bodies=[notch1_tm_bead, dll4_tm_bead])

#TWEAK 3 make MAX_TRANS higher

# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(500)

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
        root_hier=root_hier,                    # pass the root hierarchy
        #replica_exchange_maximum_temperature=rex_max_temp,
        monte_carlo_sample_objects=dof.get_movers(),  # pass all objects to be moved ( almost always dof.get_movers() )
        #global_output_directory='output/',      # The output directory for this sampling run.
        output_objects=output_objects,          # Items in output_objects write information to the stat file.
        monte_carlo_steps=10,                   # Number of MC steps between writing frames
        number_of_best_scoring_models=0,        # set >0 to store best PDB files (but this is slow)
        number_of_frames=num_frames)            # Total number of frames to run / write to the RMF file.
# TWEAK 4

# Ok, now we finally do the sampling!
rex.execute_macro()

# Outputs are then analyzed in a separate analysis script.
