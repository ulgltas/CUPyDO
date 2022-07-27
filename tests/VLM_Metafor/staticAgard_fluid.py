
def getParams():
    p = {}
    # Geometry of the problem
    p["Airfoil"] = ["NACA 0011.dat", "NACA 0011.dat"] # Symmetric airfoil - thickness is irrelevant for VLM
    p["Span"] = [0.762]
    p["Taper"] = [0.659]
    p["Twist"] = [0.0, 0.0]
    p["SweepLE"] = [46.74]
    p["Dihedral"] = [0.0]
    p["RootChord"] = 0.55905
    p["Offset"] = [0.0, 0.0] # X and Z offset
    p["Geofile"] = "models/AGARD445_Wing.arp"
    # Flow properties
    p["U_inf"] = 247.09
    p["Rho"] = 0.09411
    p["AoA"] = 1.0
    p["TimeSteps"] = 40
    p["TimeDenominator"] = 1.0
    p["FreeWake"] = 0
    p["Infile"] = "AGARD445_fluid.arp"
    # Mesh
    p["ChordwisePanels"] = 10
    p["SpanwisePanels"] = 25

    return p

