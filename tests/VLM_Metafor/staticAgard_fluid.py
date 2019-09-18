import wing
import pythonVLM.VLM_inputs as inputs
import os

airfoils = ["NACA 0011.dat", "NACA 0011.dat"]
span = [0.762]
taper = [0.659]
twist = [0.0, 0.0]
sweep_LE = [46.74]
dihedral = [0.0]
root_chord = 0.55905

offset = [0.0, 0.0]
if not os.path.exists("models"):
    os.mkdir("models") # Needs to exist to write airfoil data

w = inputs.VLMWing(airfoils, span, taper, sweep_LE, dihedral, twist, root_chord, offset)

v_airfoils = []

v_span = []
v_taper = []
v_sweep = []
v_dihedral = []
v_twist = []
v_root_chord = 0.5
v_offset = [0.0, 0.0]

vtail = inputs.VLMVTail(v_airfoils, v_span, v_taper, v_sweep, v_dihedral, v_twist, v_root_chord, v_offset)

h_airfoils = []
h_span = []
h_taper = []
h_sweep = []
h_dihedral = []
h_twist = []
h_root_chord = 0.5
h_offset = [0., 0.]


htail = inputs.VLMHTail(h_airfoils, h_span, h_taper, h_sweep, h_dihedral, h_twist, h_root_chord, h_offset)

w.write_geofile("models/AGARD445_Wing.arp")

properties = inputs.VLMProperties(w, htail, vtail)
properties.u = 247.09
properties.rho = 0.09411
properties.AoA = 1.0
properties.timesteps = 40
properties.infile = "AGARD445_fluid.arp"

properties.wing.chordwise_panels = 10
properties.wing.spanwise_panels = 25
properties.wing.geometry_file = "models/AGARD445_Wing.arp"

properties.htail.chordwise_panels = 0
properties.htail.spanwise_panels = 0
properties.htail.geometry_file = "models/AGARD445_HTail.arp"

properties.vtail.chordwise_panels = 0
properties.vtail.spanwise_panels = 0
properties.vtail.geometry_file = "models/AGARD445_VTail.arp"


properties.write_infile()
