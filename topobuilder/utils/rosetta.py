"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""

# Standard Libraries
import textwrap
import math
from typing import Dict, Union, Optional, List

# External Libraries
from jinja2 import Template
from bs4 import BeautifulSoup

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore


__all__ = ['rosettascript', 'get_weight_patches', 'funfoldes', 'hybridize', 'constraint_design']


class ScriptPieces( dict ):
    def __add__( self, other ):
        data = ScriptPieces()
        for k in ['scorefxns', 'residueselectors', 'packerpalette', 'taskoperations', 'movemapfactory',
                  'simplemetrics', 'filters', 'movers', 'protocols', 'output']:
            if k in self or k in other:
                data.setdefault(k, [x for x in self.get(k, ['', ]) + other.get(k, ['', ]) if len(x) > 0])
        return data

def get_weight_patches():
    """Defines the weights sets used during folding.
    """
    wts0 = textwrap.dedent("""\
    # score0 from rosetta++, used in stage 1 of the
    # ClassicAbinitio protocol.
    # Score 0 has a vdw weight of 1, in R++, but then it divides
    # the vdw score by 10 and rounds down to nearest integer.
    # Mini does not round down to the nearest integer.
    env     0.0
    pair    0.0
    cbeta   0.0
    vdw     0.1
    rg      0.0
    cenpack 0.0
    hs_pair 0.0
    ss_pair 0.0
    rsigma  0.0 """)

    wts0_patch = textwrap.dedent("""\
    env     = 0.0
    pair    = 0.0
    cbeta   = 0.0
    vdw     = 0.1
    rg      = 0.0
    cenpack = 0.0
    hs_pair = 0.0
    ss_pair = 0.0
    rsigma  = 0.0
    hbond_sr_bb = 0.3
    hbond_lr_bb = 0.7
    rsigma  = 0.2
    sheet   = 0.2
    ss_pair = 0.2
    hs_pair = 0.2 """)

    wts1 = textwrap.dedent("""\
    # score1 from rosetta++, used in stage 2 of ClassicAbinitio
    # protocol from Rosetta++.
    env     1.0
    pair    1.0
    cbeta   0.0
    vdw     1.0
    rg      0.0
    cenpack 0.0
    hs_pair 1.0
    ss_pair 0.3
    rsigma  0.0
    sheet   1.0
    STRAND_STRAND_WEIGHTS 1 11 """)

    wts1_patch = textwrap.dedent("""\
    env     = 1.0
    pair    = 1.0
    cbeta   = 0.0
    vdw     = 1.0
    rg      = 0.0
    cenpack = 0.0
    rsigma  = 1.1
    sheet   = 1.0
    ss_pair = 1.0
    hs_pair = 1.0
    hbond_sr_bb = 1.17
    hbond_lr_bb = 2.0
    STRAND_STRAND_WEIGHTS 1 3 """)

    wts2 = textwrap.dedent("""\
    # score2 from rosetta++, used in stage 3 of ClassicAbinitio protocol.
    env     1.0
    pair    1.0
    cbeta   0.25
    cenpack 0.5
    vdw     1.0
    rg      0.0
    hs_pair 1.0
    ss_pair 1.0
    rsigma  0.0
    sheet   1.0
    STRAND_STRAND_WEIGHTS 1 6 """)

    wts2_patch = textwrap.dedent("""\
    env     = 1.0
    pair    = 1.0
    cbeta   = 0.25
    cenpack = 0.5
    vdw     = 1.0
    rg      = 0.0
    rsigma  = 0.7
    sheet   = 0.7
    ss_pair = 0.7
    hs_pair = 0.7
    hbond_sr_bb = 1.17
    hbond_lr_bb = 2.0
    angle_constraint = 0.3
    dihedral_constraint = 0.3
    STRAND_STRAND_WEIGHTS 1 3 """)

    wts3 = textwrap.dedent("""\
    # score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol.
    env     1.0
    pair    1.0
    cbeta   1.0
    vdw     1.0
    rg      3.0
    cenpack 1.0
    hs_pair 1.0
    ss_pair 1.0
    rsigma  1.0
    sheet   1.0
    STRAND_STRAND_WEIGHTS 1 6 """)

    wts3_patch = textwrap.dedent("""\
    env     = 1.0
    pair    = 1.0
    cbeta   = 1.0
    vdw     = 1.0
    rg      = 2.0
    cenpack = 1.0
    hbond_sr_bb = 1.17
    rsigma  = 1.17
    sheet   = 1.17
    ss_pair = 1.17
    hs_pair = 1.17
    hbond_lr_bb = 2.0
    angle_constraint = 0.3
    dihedral_constraint = 0.3
    STRAND_STRAND_WEIGHTS 1 2 """)

    wts5 = textwrap.dedent("""\
    # score5.wts, used in stage 3 of ClassicAbinitio protocol
    env     1.0
    pair    1.0
    cbeta   0.25
    cenpack 0.5
    hs_pair 1.0
    ss_pair 1.0
    rsigma  0.0
    sheet   1.17
    rg      0.0
    vdw     1.0
    STRAND_STRAND_WEIGHTS 1 6 """)

    wts5_patch = textwrap.dedent("""\
    env     = 1.0
    pair    = 1.0
    cbeta   = 0.25
    cenpack = 0.5
    rg      = 0.0
    vdw     = 1.0
    rsigma  = 1.17
    sheet   = 1.17
    ss_pair = 1.17
    hs_pair = 1.17
    hbond_sr_bb = 1.17
    hbond_lr_bb = 2.0
    angle_constraint = 0.3
    dihedral_constraint = 0.3
    STRAND_STRAND_WEIGHTS 1 3 """)

    wts_pieces = [wts0, wts0_patch,
                  wts1, wts1_patch,
                  wts2, wts2_patch,
                  wts3, wts3_patch,
                  wts5, wts5_patch,]
    return wts_pieces


def constraint_minimization(  case: Case, natbias: float ) -> ScriptPieces:
    """Creates the minimization script upon a :class:`.Case` with constraints.

    :param case: The case to use for constraints.
    :param natbias: The bias weight to favour correct SSE formation.
    """
    scorefxns = [textwrap.dedent("""\
    <ScoreFunction name="sfxn_cstmin" weights="ref2015">
        <Reweight scoretype="atom_pair_constraint" weight="1.0" />
        <Reweight scoretype="angle_constraint" weight="1.0" />
        <Reweight scoretype="dihedral_constraint" weight="1.0" />
        <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    """), ]

    residueselectors = [SELECTOR_SecondaryStructure('sse_cstmin', case), ]

    filters = [textwrap.dedent("""\
    <RmsdFromResidueSelectorFilter name="rmsd_cstmin" reference_selector="sse_cstmin"
            reference_name="eminPose_cstmin" query_selector="sse_cstmin" confidence="0." />
    """), ]

    movers = [textwrap.dedent("""\
    <SavePoseMover name="spose_cstmin" reference_name="eminPose_cstmin" restore_pose="0" />
    <AddConstraints name="cst_cstmin" >
        <SegmentedAtomPairConstraintGenerator name="cst_seg_cstmin" residue_selector="sse_cstmin" >
            <Outer sd="2.0" weight="1." ca_only="1"
             use_harmonic="1" unweighted="0" max_distance="40" />
        </SegmentedAtomPairConstraintGenerator>
        <AutomaticSheetConstraintGenerator name="cst_sheet_cstmin" sd="2.0" distance="6.1" />
    </AddConstraints>
    <MinMover name="fast_cstmin" scorefxn="sfxn_cstmin" chi="1" bb="1" />
    """), MOVER_SetSecStructEnergies( 'ssse_cstmin', 'sfxn_cstmin', natbias, case )]

    protocols = [textwrap.dedent("""\
    <Add mover="spose_cstmin" />
    <Add mover="cst_cstmin" />
    <Add mover="ssse_cstmin" />
    <Add mover="fast_cstmin" />
    <Add filter="rmsd_cstmin" />
    """), ]

    bf = PROTOCOL_BasicFilters(case, '_cstmin')
    return ScriptPieces({'scorefxns': scorefxns, 'movers': movers, 'filters': filters,
                         'residueselectors': residueselectors, 'protocols': protocols}) + bf


def constraint_design( case: Case, natbias: float, layer_design: bool = True,
                       motif: Optional[List] = None,
                       binder: Optional[List] = None,
                       hotspots: Optional[List] = None,
                       profile: Optional[bool] = False,
                      ) -> ScriptPieces:
    """Creates the design script upon a :class:`.Case` with constraints.

    :param case: The case to use for constraints.
    :param motif: The motif to be inserted and kept fixed.
    :param binder: The context structure to be added during design.
    :param hotspots: Motif hotspot residue to keep completely fixed.
    :param profile: Sequence profile from fragments to guide design (default: False).
    :param natbias: The bias weight to favour correct SSE formation.
    """
    scorefxns = [textwrap.dedent("""\
    <ScoreFunction name="sfxn_cstdes" weights="ref2015">
        <Reweight scoretype="atom_pair_constraint" weight="1.0" />
        <Reweight scoretype="angle_constraint" weight="1.0" />
        <Reweight scoretype="dihedral_constraint" weight="1.0" />
        <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    <ScoreFunction name="sfxn_cstdes_cart" weights="ref2015_cart">
        <Reweight scoretype="atom_pair_constraint" weight="1.0" />
        <Reweight scoretype="angle_constraint" weight="1.0" />
        <Reweight scoretype="dihedral_constraint" weight="1.0" />
        <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    """), ]
    if motif and binder and hotspots: # all
        coldspots = [s for s in motif if s not in hotspots]
        template  = list(set([s[-1] for s in motif]))
        residueselectors = [textwrap.dedent("""<Index name="piece_ffd" resnums="{}" />""").format(','.join(motif)),
                            textwrap.dedent("""<Chain name="template" chains="{}" />""").format(','.join(template)),
                            textwrap.dedent("""<Not name="binder" selector="template"/>"""),
                            textwrap.dedent("""<Index name="hotspots" resnums="{}" />""").format(','.join(hotspots)),
                            textwrap.dedent("""<Index name="coldspots" resnums="{}" />""").format(','.join(coldspots)),
                            SELECTOR_SecondaryStructure('sse_cstdes', case)]

        taskoperations = [textwrap.dedent("""\
        <OperateOnResidueSubset name="no_repacking_hotspots" selector="hotspots" >
            <PreventRepackingRLT/>
        </OperateOnResidueSubset>
        <OperateOnResidueSubset name="no_repacking_binder" selector="binder" >
            <PreventRepackingRLT/>
        </OperateOnResidueSubset>
        <OperateOnResidueSubset name="no_C_coldspots" selector="coldspots" >
            <DisallowIfNonnativeRLT disallow_aas="C" />
        </OperateOnResidueSubset>
        <OperateOnResidueSubset name="no_CM_template" selector="template" >
            <DisallowIfNonnativeRLT disallow_aas="CM" />
        </OperateOnResidueSubset>
        """)]
        tskop = ['no_repacking_hotspots', 'no_repacking_binder', 'no_C_coldspots', 'no_CM_template']

        movemapfactory = [textwrap.dedent("""\
        <MoveMapFactory name="mmap" bb="false" chi="false" nu="false" branches="false" jumps="false" >
            <Chi      enable="true" residue_selector="coldspots" />
            <Backbone enable="true" residue_selector="template" />
            <Chi      enable="true" residue_selector="template" />
        </MoveMapFactory>
        """)]
    elif motif and hotspots: # no binder
        coldspots = [s for s in attach if s not in hotspots]
        template  = list(set([s[-1] for s in motif]))
        residueselectors = [textwrap.dedent("""<Index name="piece_ffd" resnums="{}" />""").format(','.join(motif)),
                            textwrap.dedent("""<Chain name="template" chains="{}" />""").format(','.join(template)),
                            textwrap.dedent("""<Index name="hotspots" resnums="{}" />""").format(','.join(hotspots)),
                            textwrap.dedent("""<Index name="coldspots" resnums="{}" />""").format(','.join(coldspots)),
                            SELECTOR_SecondaryStructure('sse_cstdes', case)]

        taskoperations = [textwrap.dedent("""\
        <OperateOnResidueSubset name="no_repacking_hotspots" selector="hotspots" >
            <PreventRepackingRLT/>
        </OperateOnResidueSubset>
        <OperateOnResidueSubset name="no_C_coldspots" selector="coldspots" >
            <DisallowIfNonnativeRLT disallow_aas="C" />
        </OperateOnResidueSubset>
        <OperateOnResidueSubset name="no_CM_template" selector="template" >
            <DisallowIfNonnativeRLT disallow_aas="CM" />
        </OperateOnResidueSubset>
        """)]
        tskop = ['no_repacking_hotspots', 'no_C_coldspots', 'no_CM_template']

        movemapfactory = [textwrap.dedent("""\
        <MoveMapFactory name="mmap" bb="false" chi="false" nu="false" branches="false" jumps="false" >
            <Chi      enable="true" residue_selector="coldspots" />
            <Backbone enable="true" residue_selector="template" />
            <Chi      enable="true" residue_selector="template" />
        </MoveMapFactory>
        """)]
    elif motif: # no binder and no hotspots
        template  = list(set([s[-1] for s in motif]))
        residueselectors = [textwrap.dedent("""<Index name="piece_ffd" resnums="{}"/>""").format(','.join(motif)),
                            textwrap.dedent("""<Chain name="template" chains="{}"/>""").format(','.join(template)),
                            SELECTOR_SecondaryStructure('sse_cstdes', case)]

        taskoperations = [textwrap.dedent("""\
        <OperateOnResidueSubset name="no_CM_template" selector="template">
            <DisallowIfNonnativeRLT disallow_aas="CM"/>
        </OperateOnResidueSubset>
        """)]
        tskop = ['no_CM_template']

        movemapfactory = [textwrap.dedent("""\
        <MoveMapFactory name="mmap" bb="false" chi="false" nu="false" branches="false" jumps="false">
            <Backbone enable="true" residue_selector="template"/>
            <Chi      enable="true" residue_selector="template"/>
        </MoveMapFactory>
        """)]
    else: # no motif
        residueselectors = [SELECTOR_SecondaryStructure('sse_cstdes', case), ]
        taskoperations   = None

    filters = [textwrap.dedent("""\
    <RmsdFromResidueSelectorFilter name="rmsd_cstdes" reference_selector="sse_cstdes"
            reference_name="eminPose_cstdes" query_selector="sse_cstdes" confidence="0." />
    """), ]

    movers = [textwrap.dedent("""\
    <SavePoseMover name="spose_cstdes" reference_name="eminPose_cstdes" restore_pose="0" />
    <AddConstraints name="cst_cstdes" >
        <SegmentedAtomPairConstraintGenerator name="cst_seg_cstdes" residue_selector="sse_cstdes" >
            <Outer sd="2.0" weight="1.0" ca_only="1" use_harmonic="1" unweighted="0" max_distance="40" />
        </SegmentedAtomPairConstraintGenerator>
        <!--AutomaticSheetConstraintGenerator name="cst_sheet_cstdes" sd="2.0" distance="6.1" /-->
    </AddConstraints>
    <RemoveConstraints name="rm_cstdes" constraint_generators="cst_seg_cstdes" />
    <StructFragmentMover name="makeFrags_ffd" prefix="frags" small_frag_file="{}" large_frag_file="{}" />
    <NubInitioLoopClosureMover name="close_loops" fragments_id="frags"
      break_side_ramp="true" design="true" fullatom_scorefxn="sfxn_cstdes" max_trials="100"/>
    """).format(*case['metadata.fragments.files']),
    MOVER_SetSecStructEnergies( 'ssse_cstdes', 'sfxn_cstdes', natbias, case ),
    MOVER_SetSecStructEnergies( 'ssse_cstdes_cart', 'sfxn_cstdes_cart', natbias, case )]
    if profile is True:
        movers.append(textwrap.dedent("""\
        <FavorSequenceProfile name="set_profile" scaling="none" weight="1" pssm="{}" scorefxns="sfxn_cstdes,sfxn_cstdes_cart" />
        """).format(case['metadata.fragments.profile']))
    if not taskoperations:
        movers.append(textwrap.dedent("""\
        <FastDesign name="design_cstdes" scorefxn="sfxn_cstdes_cart" relaxscript="MonomerDesign2019"
                task_operations="layer_design" ramp_down_constraints="false" repeats="5" dualspace="true"/>
        <FastRelax name="cst_cstrel" scorefxn="sfxn_cstdes_cart" repeats="5" cartesian="true"/>
        """)) if layer_design else movers.append(textwrap.dedent("""\
        <FastDesign name="design_cstdes" scorefxn="sfxn_cstdes_cart" relaxscript="MonomerDesign2019"
                ramp_down_constraints="false" repeats="5" dualspace="true"/>
        <FastRelax name="cst_cstrel" scorefxn="sfxn_cstdes_cart" repeats="5" cartesian="true"/>"""))
    else:
        movers.append(textwrap.dedent("""\
    <FastDesign name="design_cstdes" scorefxn="sfxn_cstdes_cart" relaxscript="MonomerDesign2019"
          task_operations="layer_design,{}"  movemap_factory="mmap"
          ramp_down_constraints="false" repeats="5" dualspace="true"/>
    <FastRelax name="cst_cstrel" scorefxn="sfxn_cstdes_cart"  movemap_factory="mmap" repeats="5" cartesian="true"/>
    """).format(','.join(tskop)) ) if layer_design else movers.append(textwrap.dedent("""\
    <FastDesign name="design_cstdes" scorefxn="sfxn_cstdes_cart" relaxscript="MonomerDesign2019"
          task_operations="{}"  movemap_factory="mmap"
          ramp_down_constraints="false" repeats="5" dualspace="true"/>
    <FastRelax name="cst_cstrel" scorefxn="sfxn_cstdes_cart"  movemap_factory="mmap" repeats="5" cartesian="true"/>
          """).format(','.join(tskop)) )

    if binder:
        movers.append(textwrap.dedent("""\
        <DeleteRegionMover name="remove_chains" residue_selector="binder"/>"""))

    protocols = []
    if profile is True:
        protocols.append(textwrap.dedent("""\
        <Add mover="set_profile"/>"""))
    protocols.append(textwrap.dedent("""\
    <Add mover="ssse_cstdes" />
    <Add mover="ssse_cstdes_cart" />
    <Add mover="spose_cstdes" />
    <Add mover="cst_cstdes" />
    <Add mover="design_cstdes" />
    <Add mover="rm_cstdes" />
    <Add mover="cst_cstrel" />
    <Add mover="makeFrags_ffd" />
    <Add mover="close_loops" />
    <!--Add filter="rmsd_cstdes" /-->
    """))

    ld = PROTOCOL_LayerDesign(case) if layer_design else ScriptPieces()
    bf = PROTOCOL_BasicFilters(case, '_cstdes')
    if not taskoperations:
        return ScriptPieces({'scorefxns': scorefxns, 'movers': movers, 'filters': filters,
                             'residueselectors': residueselectors, 'protocols': protocols}) + ld + bf
    else:
        return ScriptPieces({'scorefxns': scorefxns, 'movers': movers, 'filters': filters,
                             'taskoperations': taskoperations, 'movemapfactory': movemapfactory,
                             'residueselectors': residueselectors, 'protocols': protocols}) + ld + bf


def funfoldes( case: Case,
               motif: Optional[List] = None,
               binder: Optional[List] = None,
               hotspots: Optional[List] = None
              ) -> str:
    """
    The general FunFoldDesMover script where input pose is the template, where the motif will be grafted into.
    Thus, if a binder is present, one needs to remove that one before running the NubInitioMover.
    The binder will be parsed together with the motif. Note, that after the NubInitioMover, the binder will
    be in the pose. We remove the binder as it is easier to work with files containing only the design afterwards.

    :param case: The case to use for constraints.
    :param motif: The motif to be inserted and kept fixed.
    :param binder: The context structure to be added during design.
    :param hotspots: Motif hotspot residue to keep completely fixed.
    """
    if motif and binder: # binder loaded after
        coldspots = [s[:-1] for s in motif if s not in hotspots]
        template = list(set([s[-1] for s in motif]))
        residueselectors = [textwrap.dedent("""<Index name="piece_ffd" resnums="{}"/>""").format(','.join(motif)),
                            textwrap.dedent("""<Chain name="template" chains="{}"/>""").format(','.join(template)),
                            textwrap.dedent("""<Not name="binder" selector="template"/>"""),
                            SELECTOR_SecondaryStructure('sse_ffd', case)]
    elif motif: # binder loaded after
        coldspots = [s[:-1] for s in motif if s not in hotspots]
        residueselectors = [textwrap.dedent("""<Index name="piece_ffd" resnums="{}"/>""").format(','.join(motif)),
                            SELECTOR_SecondaryStructure('sse_ffd', case)]
    else:
        mid = 2 # math.floor(len(case.secondary_structure) / 2)
        residueselectors = [textwrap.dedent("""<Index name="piece_ffd" resnums="{}-{}"/>""").format(mid - 1, mid + 1),
                            SELECTOR_SecondaryStructure('sse_ffd', case)]

    filters = [textwrap.dedent("""\
        <RmsdFromResidueSelectorFilter name="rmsd_ffd" reference_selector="sse_ffd"
            reference_name="sketchPose_ffd" query_selector="sse_ffd" confidence="0."/>""")]

    movers = [MOVER_PeptideStubMover('add_loops_ffd', case), textwrap.dedent("""\
        <SavePoseMover name="save_ffd" reference_name="sketchPose_ffd" restore_pose="0"/>
        <StructFragmentMover name="makeFrags_ffd" prefix="frags" small_frag_file="{}" large_frag_file="{}"/>
        <AddConstraints name="foldingCST_ffd">
            <SegmentedAtomPairConstraintGenerator name="foldCST" residue_selector="sse_ffd">
                <!--Inner sd="1.2" weight="1." ca_only="1"
                    use_harmonic="true" unweighted="false" min_seq_sep="4"/-->
                <Outer sd="2" weight="2." ca_only="1"
                    use_harmonic="true" unweighted="false"  max_distance="40"/>
            </SegmentedAtomPairConstraintGenerator>
            <!--AutomaticSheetConstraintGenerator name="sheetCST" sd="2.0" distance="6.1"/-->
        </AddConstraints>
        """).format(*case['metadata.fragments.files'])]
    if motif and template:
        movers.append(textwrap.dedent("""\
            <NubInitioMover name="FFL_ffd" fragments_id="frags" template_motif_selector="piece_ffd" rmsd_threshold="10"
                correction_weights="0" clear_motif_cst="0">
                <Nub reference_name="sketchPose_ffd" residue_selector="piece_ffd" binder_selector="binder">
                    <!--Segment order="1" n_term_flex="1" c_term_flex="1" editable="{}"/-->
                </Nub>
            """).format(','.join(coldspots)))
        movers.append(textwrap.dedent("""</NubInitioMover>"""))
        movers.append(textwrap.dedent("""<DeleteRegionMover name="remove_chains" residue_selector="binder"/>"""))
    elif motif:
        movers.append( textwrap.dedent("""\
            <NubInitioMover name="FFL_ffd" fragments_id="frags" template_motif_selector="piece_ffd" rmsd_threshold="10"
                correction_weights="0" clear_motif_cst="0">
                <Nub reference_name="sketchPose_ffd" residue_selector="piece_ffd">
                    <!--Segment order="1" n_term_flex="1" c_term_flex="1" editable="{}"/-->
                </Nub>
            """).format(','.join(coldspots)))
        movers.append(textwrap.dedent("""</NubInitioMover>"""))
    else:
        movers.append( textwrap.dedent("""\
            <NubInitioMover name="FFL_ffd" fragments_id="frags" template_motif_selector="piece_ffd" rmsd_threshold="10"
                correction_weights="0" clear_motif_cst="0">
                <Nub reference_name="sketchPose_ffd" residue_selector="piece_ffd">
                    <Segment order="1" n_term_flex="2" c_term_flex="1" editable="1,2,3"/>
                </Nub>
            """) )
        movers.append(textwrap.dedent("""</NubInitioMover>"""))
    if binder:
        protocols = [textwrap.dedent("""\
            <Add mover="add_loops_ffd"/>
            <Add mover="save_ffd"/>
            <Add mover="remove_chains"/>
            <Add mover="makeFrags_ffd"/>
            <Add mover="foldingCST_ffd"/>
            <Add mover="FFL_ffd"/>
            <!--Add mover="remove_chains"/-->
            <!--Add filter="rmsd_ffd"/-->""")]
    else:
        protocols = [textwrap.dedent("""\
            <Add mover="add_loops_ffd"/>
            <Add mover="save_ffd"/>
            <Add mover="makeFrags_ffd"/>
            <Add mover="foldingCST_ffd"/>
            <Add mover="FFL_ffd"/>
            <!--Add mover="remove_chains"/-->
            <!--Add filter="rmsd_ffd"/-->""")]

    with TBcore.on_option_value('psipred', 'script', None):
        bf = PROTOCOL_BasicFilters(case, '_ffd')
    return ScriptPieces({'movers': movers, 'filters': filters,
                         'residueselectors': residueselectors, 'protocols': protocols}) + bf


def hybridize( case: Case, template: str, natbias: float ) -> str:
    """
    Creates the hybridize script upon a :class:`.Case` with constraints.

    .. caution::
        Please be careful with the energy function setup.

    .. admonition:: To Developers

        This implementation has not been fully tested and optimized.

    :param case: The case to use for constraints.
    :param template: The template apply the hybridization onto.
    :param natbias: The bias weight to favour correct SSE formation.
    """
    #mid = math.floor(len(case.secondary_structure) / 2)
    mid = 2
    scorefxns = [textwrap.dedent("""\
    <ScoreFunction name="fa" weights="ref2015">
      <Reweight scoretype="coordinate_constraint" weight="1" />
      <Reweight scoretype="atom_pair_constraint" weight="1" />
      <Reweight scoretype="dihedral_constraint" weight="1" />
      <Reweight scoretype="angle_constraint" weight="1" />
      <Reweight scoretype="rsigma" weight="1.17" />
      <Reweight scoretype="sheet" weight="2.0" />
      <Reweight scoretype="ss_pair" weight="2.0" />
      <Reweight scoretype="hs_pair" weight="2.0" />
      <Reweight scoretype="hbond_sr_bb" weight="1.17" />
      <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    <ScoreFunction name="stage1" weights="score3">
      <Reweight scoretype="coordinate_constraint" weight="1" />
      <Reweight scoretype="atom_pair_constraint" weight="1" />
      <Reweight scoretype="dihedral_constraint" weight="1" />
      <Reweight scoretype="angle_constraint" weight="1" />
      <Reweight scoretype="rsigma" weight="1.17" />
      <Reweight scoretype="sheet" weight="2.0" />
      <Reweight scoretype="ss_pair" weight="2.0" />
      <Reweight scoretype="hs_pair" weight="2.0" />
      <Reweight scoretype="hbond_sr_bb" weight="1.17" />
      <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    <ScoreFunction name="stage2" weights="score4_smooth_cart">
      <Reweight scoretype="coordinate_constraint" weight="1" />
      <Reweight scoretype="atom_pair_constraint" weight="1" />
      <Reweight scoretype="dihedral_constraint" weight="1" />
      <Reweight scoretype="angle_constraint" weight="1" />
      <Reweight scoretype="rsigma" weight="1.17" />
      <Reweight scoretype="sheet" weight="2.0" />
      <Reweight scoretype="ss_pair" weight="2.0" />
      <Reweight scoretype="hs_pair" weight="2.0" />
      <Reweight scoretype="hbond_sr_bb" weight="1.17" />
      <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    """), ]

    residueselectors = [textwrap.dedent("""\
    <Index name="piece_ffd" resnums="{}-{}" />
    """).format(mid - 1, mid + 1), SELECTOR_SecondaryStructure('sse_ffd', case)]

    filters = [textwrap.dedent("""\
        <RmsdFromResidueSelectorFilter name="rmsd_ffd" reference_selector="sse_ffd"
            reference_name="sketchPose_ffd" query_selector="sse_ffd" confidence="0." />""")]

    movers = [MOVER_PeptideStubMover('add_loops_ffd', case, 'VAL'),
    #MOVER_SetSecStructEnergies('sse_energies', 'fa', natbias, case),
    textwrap.dedent("""\
        <SavePoseMover name="save_ffd" reference_name="sketchPose_ffd" restore_pose="0" />
        <AddConstraints name="foldingCST_ffd" >
            <SegmentedAtomPairConstraintGenerator name="foldCST" residue_selector="sse_ffd" >
                <!--Inner sd="1.2" weight="1." ca_only="1"
                    use_harmonic="true" unweighted="false" min_seq_sep="4" /-->
                <Outer sd="2" weight="2." ca_only="1"
                    use_harmonic="true" unweighted="false"  max_distance="40" />
            </SegmentedAtomPairConstraintGenerator>
            <AutomaticSheetConstraintGenerator name="sheetCST" sd="2.0" distance="6.1" />
        </AddConstraints>"""),
        MOVER_SetSecStructEnergies( 'ssse_stage1', 'stage1', natbias, case ),
        MOVER_SetSecStructEnergies( 'ssse_stage2', 'stage2', natbias, case ),
        MOVER_SetSecStructEnergies( 'ssse_fa', 'fa', natbias, case ),
    textwrap.dedent("""\
        <Hybridize name="hybridize" stage1_scorefxn="stage1" stage2_scorefxn="stage2" fa_scorefxn="fa"
                   batch="1" stage1_increase_cycles="1.0" stage2_increase_cycles="1.0" linmin_only="0"
                   realign_domains="0" keep_pose_constraint="1" csts_from_frags="0" max_registry_shift="1">
                   <Fragments three_mers="{}" nine_mers="{}"/>
                   <Template pdb="{}" cst_file="AUTO" weight="1.000" />
        </Hybridize>""").format(*case['metadata.fragments.files'], template)]

    protocols = [textwrap.dedent("""\
        <Add mover="add_loops_ffd" />
        <Add mover="save_ffd" />
        <Add mover="foldingCST_ffd" />
        <Add mover="ssse_stage1" />
        <Add mover="ssse_stage2" />
        <Add mover="ssse_fa" />
        <Add mover="hybridize" />
        <Add filter="rmsd_ffd" />""")]

    with TBcore.on_option_value('psipred', 'script', None):
        bf = PROTOCOL_BasicFilters(case, '_ffd')
    return ScriptPieces({'scorefxns': scorefxns, 'movers': movers, 'filters': filters,
                         'residueselectors': residueselectors,
                         'protocols': protocols}) + bf


def rosettascript( data: Dict ) -> str:
    """Global RosettaScripts XML template script.
    """
    if 'output' in data and isinstance(data['output'], list):
        data['output'] = data['output'][0]
    content = BeautifulSoup(Template(textwrap.dedent("""\
    <ROSETTASCRIPTS>
        {% if scorefxns %}<SCOREFXNS>{% for item in scorefxns %}{{ item }}{% endfor %}</SCOREFXNS>{% endif %}
        {% if residueselectors %}
            <RESIDUE_SELECTORS>{% for item in residueselectors %}{{ item }}{% endfor %}</RESIDUE_SELECTORS>
        {% endif %}
        {% if packerpalettes %}<PACKER_PALETTES>{% for item in packerpalettes %}{{item}}{% endfor %}</PACKER_PALETTES>{% endif %}
        {% if taskoperations %}<TASKOPERATIONS>{% for item in taskoperations %}{{ item }}{% endfor %}</TASKOPERATIONS>{% endif %}
        {% if movemapfactory %}
            <MOVE_MAP_FACTORIES>{% for item in movemapfactory %}{{ item }}{% endfor %}</MOVE_MAP_FACTORIES>
        {% endif %}
        {% if simplemetrics %}<SIMPLE_METRICS>{% for item in simplemetrics %}{{ item }}{% endfor %}</SIMPLE_METRICS>{% endif %}
        {% if filters %}<FILTERS>{% for item in filters %}{{ item }}{% endfor %}</FILTERS>{% endif %}
        {% if movers %}<MOVERS>{% for item in movers %}{{ item }}{% endfor %}</MOVERS>{% endif %}
        {% if protocols %}<PROTOCOLS>{% for item in protocols %}{{ item }}{% endfor %}</PROTOCOLS>{% endif %}
        {% if output %}<OUTPUT scorefxn="{{ output }}"/>{% endif %}
    </ROSETTASCRIPTS>
    """)).render(data), 'xml').contents
    return '\n'.join(x.prettify() for x in content)


def SELECTOR_SecondaryStructure( name: str,
                                 sse: Union[str, Case],
                                 sse_type: str = 'HE',
                                 terminal_loops: bool = False
                                 ) -> str:
    """Selects the secondary structure elements (SSEs).

    :param name: Name to give to the selector.
    :param sse: The SSEs to take into account.
    :param sse_type: The SSE types to take into account,
        either helix (H) or strands (E) (default: HE).
    :param terminal_loops: Should terminal loops be included (default: False).
    """
    if isinstance(sse, Case):
        sse = sse.secondary_structure

    return textwrap.dedent("""\
        <SecondaryStructure name="{}" overlap="0" minE="1"  minH="1" ss="{}" include_terminal_loops="{}"
        use_dssp="0" pose_secstruct="{}" />""").format(name, sse_type, int(terminal_loops), sse)


def MOVER_SetSecStructEnergies( name: str, score: str, natbias: float, case: Case ) -> str:
    """Sets the SSE energy bias term.

    :param name: Name to give to the mover.
    :param score: The scorefunction to be modified.
    :param natbias: The bias weight.
    :param case: Case to take into account for SSE selection.
    """
    data = dict(zip(['ss_pair', 'hh_pair', 'hss_triplets'], case.sse_pairing))
    data['sse'] = case.secondary_structure
    data['score'] = score
    data['name'] = name
    data['natbias'] = natbias
    return Template(textwrap.dedent("""\
        <SetSecStructEnergies name="{{name}}" scorefxn="{{score}}"
            secstruct="{{sse}}" use_dssp="0"
            {% if hh_pair|length > 0 %}hh_pair="{{hh_pair}}"{% endif %}
            {% if ss_pair|length > 0 %}ss_pair="{{ss_pair}}"{% endif %}
            {% if hss_triplets|length > 0 %}hss_triplets="{{hss_triplets}}"{% endif %}
            {% if ss_pair|length > 0 %}natbias_ss="{{natbias}}"{% endif %}
            {% if hh_pair|length > 0 %}natbias_hh="{{natbias}}"{% endif %}
            {% if hss_triplets|length > 0 %}natbias_hs="{{natbias}}"{% endif %}
        />""")).render(data)


def MOVER_PeptideStubMover( name: str, case: Case, residue: str = 'VAL' ) -> str:
    """Adds residues to an existing pose.

    :param name: Name to give to the mover.
    :param residue: The residue type that shall be added (default: 'VAL').
    :param case: Case to add residues to.
    """
    return Template(textwrap.dedent("""\
        <PeptideStubMover name="{{name}}" reset="false">
        {% for item in insert %}
        <Insert resname="{{residue}}" repeat="1" jump="false" anchor_rsd="{{item}}" anchor_atom="C" connecting_atom="N" />
        {% endfor %}
        </PeptideStubMover>""")).render({'insert': [i for i, ltr in enumerate(case.secondary_structure) if ltr == 'L'],
                                         'name': name, 'residue': residue})


def PROTOCOL_LayerDesign( case: Case ) -> ScriptPieces:
    """Creates the layer design protocol, where a protein is split into a core, boundary and
    surface layer where different AA are prefered. This happens on a per-residue level.

    :param case: Case to calculate the layers from.
    """
    residueselectors = [textwrap.dedent("""\
        <Layer name="surface" select_core="0" select_boundary="0" select_surface="1" use_sidechain_neighbors="1"/>
        <Layer name="boundary" select_core="0" select_boundary="1" select_surface="0" use_sidechain_neighbors="1"/>
        <Layer name="core" select_core="1" select_boundary="0" select_surface="0" use_sidechain_neighbors="1"/>"""),
                        SELECTOR_SecondaryStructure('sheet', case, 'E'),
                        SELECTOR_SecondaryStructure('entire_helix', case, 'H'),
                        SELECTOR_SecondaryStructure('entire_loop', case, 'L', True), textwrap.dedent("""\
        <And name="helix_cap" selectors="entire_loop">
            <PrimarySequenceNeighborhood lower="1" upper="0" selector="entire_helix"/></And>
        <And name="helix_start" selectors="entire_helix">
            <PrimarySequenceNeighborhood lower="0" upper="1" selector="helix_cap"/></And>
        <And name="helix" selectors="entire_helix"><Not selector="helix_start"/></And>
        <And name="loop" selectors="entire_loop"><Not selector="helix_cap"/></And>
        """)]

    taskoperations = [textwrap.dedent("""\
        <DesignRestrictions name="layer_design">
            <Action selector_logic="surface AND helix_start" aas="DEHKPQR"/>
            <Action selector_logic="surface AND helix" aas="EHKQR"/>
            <Action selector_logic="surface AND sheet" aas="EHKNQRST"/>
            <Action selector_logic="surface AND loop" aas="DEGHKNPQRST"/>
            <Action selector_logic="boundary AND helix_start" aas="ADEHIKLMNPQRSTVWY"/>
            <Action selector_logic="boundary AND helix" aas="ADEHIKLMNQRSTVWY"/>
            <Action selector_logic="boundary AND sheet" aas="DEFHIKLMNQRSTVWY"/>
            <Action selector_logic="boundary AND loop" aas="ADEFGHIKLMNPQRSTVWY"/>
            <Action selector_logic="core AND helix_start" aas="AFILMPVWY"/>
            <Action selector_logic="core AND helix" aas="AFILMVWY"/>
            <Action selector_logic="core AND sheet" aas="FILMVWY"/>
            <Action selector_logic="core AND loop" aas="AFGILMPVWY"/>
            <Action selector_logic="helix_cap" aas="DNST"/>
        </DesignRestrictions>"""), ]

    return ScriptPieces({'residueselectors': residueselectors, 'taskoperations': taskoperations})


def PROTOCOL_BasicFilters( case: Case, suffix: str = '' ) -> ScriptPieces:
    """A wrapper protocol around basic filters and score terms.

    :param case: Case to be used for filter and score calculations.
    :param suffix: The suffix to be appended to the output score names.
    """
    sse = case.secondary_structure

    scorefxns = textwrap.dedent("""\
        <ScoreFunction name="bb_only" weights="empty.wts" >
          <Reweight scoretype="fa_rep" weight="0.1" />
          <Reweight scoretype="fa_atr" weight="0.2" />
          <Reweight scoretype="hbond_sr_bb" weight="2.0" />
          <Reweight scoretype="hbond_lr_bb" weight="2.0" />
          <Reweight scoretype="rama_prepro" weight="0.45" />
          <Reweight scoretype="omega" weight="0.4" />
          <Reweight scoretype="p_aa_pp" weight="0.6" />
        </ScoreFunction>
    """)

    residueselectors = textwrap.dedent("""\
        <True name="full_pose" />
        <Layer name="surface{suffix}" select_core="0" select_boundary="0" select_surface="1" use_sidechain_neighbors="1"/>
        <Layer name="boundary{suffix}" select_core="0" select_boundary="1" select_surface="0" use_sidechain_neighbors="1"/>
        <Layer name="core{suffix}" select_core="1" select_boundary="0" select_surface="0" use_sidechain_neighbors="1"/>
        """).format(suffix=suffix)

    filters = [textwrap.dedent("""\
    <PackStat name="pack{suffix}" confidence="0." />
    <CavityVolume name="cav_vol{suffix}" confidence="0." />
    <ScorePoseSegmentFromResidueSelectorFilter name="bbscore{suffix}" confidence="0"
    residue_selector="full_pose" scorefxn="bb_only" />
    """).format(sse1=sse, suffix=suffix)]

    #if not case.data['metadata']['binder']:
    if not 'binder' in case.data['metadata']:
        filters.append(textwrap.dedent("""\
        <SecondaryStructure name="sse_match{suffix}" ss="{sse1}" compute_pose_secstruct_by_dssp="true" confidence="0." />
        """).format(sse1=sse, suffix=suffix))

    movers = [textwrap.dedent("""\
    <LabelPoseFromResidueSelectorMover name="labelcore{suffix}" property="CORE" residue_selector="core{suffix}" />
    <LabelPoseFromResidueSelectorMover name="labelboundary{suffix}" property="BOUNDARY" residue_selector="boundary{suffix}" />
    <LabelPoseFromResidueSelectorMover name="labelsurface{suffix}" property="SURFACE" residue_selector="surface{suffix}" />
    <DisplayPoseLabelsMover name="labeldump{suffix}" use_dssp="1" write="1" />
    """).format(suffix=suffix), ]
    if TBcore.get_option('psipred', 'script', in_path_none=True) is not None:
        movers.append(textwrap.dedent("""\
        <WriteSSEMover name="sse_report{suffix}" cmd="{psipred}" dssp="1" write_phipsi="1" />
        """).format(suffix=suffix, psipred=TBcore.get_option('psipred', 'script')))
    else:
        movers.append(textwrap.dedent("""\
        <WriteSSEMover name="sse_report{suffix}" dssp="1" write_phipsi="1" />
        """).format(suffix=suffix))

    protocols = [textwrap.dedent("""\
    <Add mover="labelcore{suffix}"/>
    <Add mover="labelboundary{suffix}"/>
    <Add mover="labelsurface{suffix}"/>
    <Add mover="labeldump{suffix}"/>
    <Add mover="sse_report{suffix}"/>
    <Add filter="pack{suffix}" />
    <Add filter="cav_vol{suffix}" />
    """).format(suffix=suffix)]

    #if not case.data['metadata']['binder']:
    if not 'binder' in case.data['metadata']:
        protocols.append(textwrap.dedent("""\
        <Add filter="sse_match{suffix}" />
        """).format(suffix=suffix))

    return ScriptPieces({'scorefxns': [scorefxns, ], 'residueselectors': [residueselectors, ], 'filters': filters,
                         'movers': movers, 'protocols': protocols})
