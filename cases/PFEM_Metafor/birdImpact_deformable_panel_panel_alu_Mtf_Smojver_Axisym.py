# -*- coding: latin-1; -*-
# bends a simple beam from a gmsh file

from wrap import *

metafor = None

def getMetafor(p={}):
    global metafor
    if metafor: return metafor
    metafor = Metafor()

    domain = metafor.getDomain()
    geometry = domain.getGeometry()
    geometry.setDimAxisymmetric()

    # import .geo
    from toolbox.gmshOld import GmshImport
    f = os.path.join(os.path.dirname(__file__), "birdImpact_deformable_panel_Mtf_Pfem_Smojver_Axisym.msh")
    importer = GmshImport(f, domain)
    importer.execute2D()

    groupset = domain.getGeometry().getGroupSet()
    p['FSI'] = groupset(17) 

    # solid elements / material
    interactionset = domain.getInteractionSet()

    app1 = FieldApplicator(1)
    app1.push( groupset(22) )  # physical group 22: panneau
    interactionset.add( app1 )

    materset = domain.getMaterialSet()
    materset.define( 1, EvpIsoHHypoMaterial)
    mater1 = materset(1)
    mater1.put(MASS_DENSITY,    2713.0)  # [kg/m³]
    mater1.put(ELASTIC_MODULUS, 68.26e9)  # [Pa]
    mater1.put(POISSON_RATIO,   0.33)   # [-] 
    mater1.put(YIELD_NUM, 1)

    materlawset = domain.getMaterialLawSet()
    materlawset.define (1, LinearIsotropicHardening)
    materlawset(1).put(IH_SIGEL, 262.0e6) # [Pa]
    materlawset(1).put(IH_H,    2590.0e6) # [Pa]
    
    prp = ElementProperties(Volume2DElement)
    app1.addProperty(prp)
    prp.put (MATERIAL, 1)
    prp.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_SRIPR)
    
    # boundary conditions
    loadingset = domain.getLoadingSet()

    loadingset.define(groupset(15), Field1D(TY,RE)) # Forbid vertical displacements on part of the panel
    
    loadingset.define(groupset(16), Field1D(TX,RE)) # Clamping
    loadingset.define(groupset(16), Field1D(TY,RE)) # Clamping
    
    loadingset.define(groupset(18), Field1D(TX,RE)) # Axisymmetry


    mim = metafor.getMechanicalIterationManager()
    mim.setResidualTolerance(1.0e-4)
    # mim.setResidualComputationMethod(Method4ResidualComputation(1000.)) 

    ti = AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)

    # visu
    if 0:
        tsm = metafor.getTimeStepManager()
        tsm.setInitialTime(0.0, 1.0)
        tsm.setNextTime(1.0, 1, 1.0)

    # results (DOES NOT WORK PROPERLY!)
    vmgr = metafor.getValuesManager()
    vmgr.add(1, TdFieldValueExtractor(metafor, groupset(22), THERMODYN_TRAV_FINT), 'panel_work_int_F')
    vmgr.add(2, TdFieldValueExtractor(metafor, groupset(22), THERMODYN_EN_CIN), 'panel_kin_En')
    vmgr.add(3, TdFieldValueExtractor(metafor, groupset(22), THERMODYN_EN_DIS), 'panel_diss_En')
    vmgr.add(4, TdFieldValueExtractor(metafor, groupset(22), THERMODYN_TRAV_FEXT), 'panel_work_ext_F')
    
    metafor.getDbTdFieldsValueManager().setComputeEXW(True)
    metafor.getDbTdFieldsValueManager().setComputeINW(True)
    
    p['exporter'] = Extractor()
    return metafor

class Extractor(object):
    def __init__(self):

        self.metafor = metafor
        

    def to_ascii(self,extractor):

        file = open(extractor.buildName()+'.ascii', 'a')
        
        file.write('{0:12.6f}\t'.format(self.metafor.getCurrentTime()))
        file.write('{0:12.6f}\n'.format(extractor.extract()[0]))
        file.close()

    def write(self):

        groupset = self.metafor.getDomain().getGeometry().getGroupSet()
        self.to_ascii(TdFieldValueExtractor(metafor, groupset(18), THERMODYN_TRAV_FINT))
        self.to_ascii(DbNodalValueExtractor(groupset(18), Field1D(TY,RE), SortByDist0(0.,0.,0.), 1))
        