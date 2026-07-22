#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.9.1"    

from time import sleep
import numpy as np
import sys
import os
import shutil
import argparse
from stb.core import structure_io, kspace
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source

REPORT_FILE = "stb_inputfile_report.txt"

# --- Internal Relax Template ---
CALC_RELAX_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include structure.fdf

#%block Geometry.Constraints
  #stress 1  # Fixes XX 
  #stress 2  # Fixes YY
  #stress 3  # Fixes ZZ 
  #stress 4  # Fixes YZ 
  #stress 5  # Fixes XZ 
  #stress 6  # Fixes XY 
#%endblock Geometry.Constraints

## Optional Parameters
  # %block Supercell
  # %endblock Supercell
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       DZP
PAO.EnergyShift     0.01 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH 
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [6  6  6]
Mesh.CutOff           320  Ry  
FilterCutoff          150  Ry

## Optional Parameters
  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-5  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.

## Optional Parameters
  # SCF.MustConverge true 
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   false

## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## STRUCTURE RELAXATION (MOLECULAR DYNAMICS - MD)
## ===================================================================
##Essential Parameters
MD.TypeOfRun            CG
MD.VariableCell         true
MD.MaxForceTol          0.01 eV/Ang
MD.Steps                150

## Optional Parameters
  # MD.RelaxCellOnly true
  # Constant.Volume false
  # MD.MaxStressTol 1 GPa
  # MD.MaxDispl 0.2 Bohr
  # MD.PreconditionVariableCell 5 Ang
  # MD.Broyden.History.Steps 5
  # MD.Broyden.Cycle.On.Mait true
  # Target.Pressure 0 GPa
  # MD.RemoveIntramáscularPressure false

## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates 
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize 
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.  #(Use previous Density Matrix if available)

## Optional Parameters (Set to false or inactive)
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.  
  # ON.UseSaveLWF  .false.  

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteMDHistory          .true.
WriteCoorStep           .true.
WriteForces             .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.

## Optional Parameters (Interface flags)
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""
# --- End of Relax Template ---


# --- Internal Total Energy Template ---
CALC_TOTAL_ENERGY_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include struct.fdf

## Optional Parameters
  # %block Supercell
  # %endblock Supercell
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       DZP
PAO.EnergyShift     0.02 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [4  13  1]
Mesh.CutOff           320  Ry   
FilterCutoff          150  Ry

## Optional Parameters
  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-5  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.

## Optional Parameters
  # SCF.MustConverge true
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   .false.
## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.
## Optional Parameters (Set to false or inactive)
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.
  # ON.UseSaveLWF  .false.

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteMDHistory          .true.
WriteCoorStep           .true.
WriteForces             .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.
## Optional Parameters (Interface flags)
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""
# --- End of Total Energy Template ---


# --- Internal AIMD Template ---
CALC_AIMD_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include struct.fdf

%block Supercell
   2   0   0
   0   2   0
   0   0   2
%endblock Supercell

#%block Geometry.Constraints
  #stress 1  # Fixes XX 
  #stress 2  # Fixes YY
  #stress 3  # Fixes ZZ 
  #stress 4  # Fixes YZ 
  #stress 5  # Fixes XZ 
  #stress 6  # Fixes XY 
#%endblock Geometry.Constraints


## Optional Parameters
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       SZ
PAO.EnergyShift     0.02 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [1  1  1]
Mesh.CutOff           320  Ry    
FilterCutoff          150  Ry

## Optional Parameters
  
  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-4  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.

## Optional Parameters
  # SCF.MustConverge true
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   .false.
## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## STRUCTURE RELAXATION (MOLECULAR DYNAMICS - MD)
## ===================================================================
##Essential Parameters
MD.TypeOfRun  Nose
MD.InitialTimeStep  1
MD.FinalTimeStep  5000
MD.LengthTimeStep  1.0  fs
MD.InitialTemperature  300.0  K
MD.TargetTemperature   300.0  K
MD.NoseMass  30.0  Ry*fs**2

## Optional Parameters
  # MD.RelaxCellOnly true
  # Constant.Volume false
  # MD.MaxStressTol 1 GPa
  # MD.MaxDispl 0.2 Bohr
  # MD.PreconditionVariableCell 5 Ang
  # MD.Broyden.History.Steps 5
  # MD.Broyden.Cycle.On.Mait true
  # Target.Pressure 0 GPa
  # MD.RemoveIntramáscularPressure false

## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.
MD.UseSaveXV            .true.
## Optional Parameters (Set to false or inactive)
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.
  # ON.UseSaveLWF  .false.

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteMDHistory          .true.
WriteCoorStep           .true.
WriteForces             .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.
## Optional Parameters (Interface flags)
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""
# --- End of AIMD Template ---


# --- Internal Band Structure Template ---
CALC_BANDS_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include struct.fdf

## Optional Parameters
  #%block Supercell  %endblock Supercell
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       DZP
PAO.EnergyShift     0.02 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [1  1  1]
Mesh.CutOff           320  Ry 
FilterCutoff          150  Ry

## Optional Parameters

  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-5  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.

## Optional Parameters
  # SCF.MustConverge true
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   .false.
## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## BAND STRUCTURE AND DENSITY OF STATES
## ===================================================================
##Essential Parameters
#___
#IMPORTANT: k-path file created by STB-SUITE.
#For manual k-path use optional
#Parameters "%block BandLines" and comment the bellow line.
%include kpath_bs.fdf
#___

%PDOS.kgrid_Monkhorst_Pack
 20     0    0  0.0
  0    20    0  0.0
  0     0   20  0.0
%end PDOS.kgrid_Monkhorst_Pack

%block  ProjectedDensityOfStates
-20.0  20.0  0.01  8000  eV
%endblock  ProjectedDensityOfStates

## Optional Parameters
#BandLinesScale  ReciprocalLatticeVectors
# %block BandLines
#  Important: Include here the k-path
# %endblock BandLines


## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.
## Optional Parameters (Set to false or inactive)
  # MD.UseSaveXV  .true.
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.
  # ON.UseSaveLWF .false.

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.
## Optional Parameters (Interface flags)
  #WriteMDHistory          .true.
  #WriteCoorStep           .true.
  #WriteForces             .true.
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""
# --- End of Band Structure Template ---


# --- Calculation Logic Functions ---

VACUUM_GAP_ANG = 10.0  # min empty span (wrapped) along an axis to treat it as vacuum-padded

BIB_FILE = "references.bib"

# Method/tool references, keyed the same way the generated calc.fdf actually
# exercises them: SIESTA itself and the PBE (XC.Functional GGA / XC.Authors
# PBE, hardcoded in every template above) are cited unconditionally; DFT-D3
# and the Nose thermostat only when the generated input file actually turns
# them on (DFTD3 .true., or an AIMD template, whose MD.TypeOfRun is Nose).
_BIB_SIESTA = ("Soler2002", """@article{Soler2002,
  author  = {Soler, Jos\\'e M. and Artacho, Emilio and Gale, Julian D. and Garc\\'ia, Alberto and Junquera, Javier and Ordej\\'on, Pablo and S\\'anchez-Portal, Daniel},
  title   = {The {SIESTA} method for ab initio order-{N} materials simulation},
  journal = {Journal of Physics: Condensed Matter},
  year    = {2002},
  volume  = {14},
  number  = {11},
  pages   = {2745--2779},
  doi     = {10.1088/0953-8984/14/11/302}
}""")

_BIB_SIESTA_RECENT = ("Garcia2020", """@article{Garcia2020,
  author  = {Garc\\'ia, Alberto and Papior, Nick and Akhtar, Arsalan and Nordstr\\"om, Emilio and Ordej\\'on, Pablo and Blum, Volker and Sanchez-Portal, Daniel},
  title   = {Siesta: Recent developments and applications},
  journal = {The Journal of Chemical Physics},
  year    = {2020},
  volume  = {152},
  number  = {20},
  pages   = {204108},
  doi     = {10.1063/5.0005077}
}""")

_BIB_PBE = ("Perdew1996", """@article{Perdew1996,
  author  = {Perdew, John P. and Burke, Kieron and Ernzerhof, Matthias},
  title   = {Generalized Gradient Approximation Made Simple},
  journal = {Physical Review Letters},
  year    = {1996},
  volume  = {77},
  number  = {18},
  pages   = {3865--3868},
  doi     = {10.1103/PhysRevLett.77.3865}
}""")

_BIB_DFTD3 = ("Grimme2010", """@article{Grimme2010,
  author  = {Grimme, Stefan and Antony, Jens and Ehrlich, Stephan and Krieg, Helge},
  title   = {A consistent and accurate ab initio parametrization of density functional dispersion correction ({DFT-D}) for the 94 elements {H-Pu}},
  journal = {The Journal of Chemical Physics},
  year    = {2010},
  volume  = {132},
  number  = {15},
  pages   = {154104},
  doi     = {10.1063/1.3382344}
}""")

_BIB_NOSE = ("Nose1984", """@article{Nose1984,
  author  = {Nos\\'e, Shuichi},
  title   = {A unified formulation of the constant temperature molecular dynamics methods},
  journal = {The Journal of Chemical Physics},
  year    = {1984},
  volume  = {81},
  number  = {1},
  pages   = {511--519},
  doi     = {10.1063/1.447334}
}""")

# Pseudopotential-bank citations, one per stb.core.pseudopotentials.BANKS key
# -- kept here rather than read from BANKS[...]['citation'] directly since
# that field is a plain human-readable string (printed by resolve_pseudo_
# source()), not already-formatted BibTeX.
_BIB_PSEUDO_DOJO = ("vanSetten2018", """@article{vanSetten2018,
  author  = {van Setten, M. J. and Giantomassi, M. and Bousquet, E. and Verstraete, M. J. and Hamann, D. R. and Gonze, X. and Rignanese, G.-M.},
  title   = {The {PseudoDojo}: Training and grading a 85 element optimized norm-conserving pseudopotential table},
  journal = {Computer Physics Communications},
  year    = {2018},
  volume  = {226},
  pages   = {39--54},
  doi     = {10.1016/j.cpc.2018.01.012}
}""")

_BIB_VIRTUAL_VAULT = ("VirtualVault", """@misc{VirtualVault,
  author       = {{NNIN/C}},
  title        = {{SIESTA} Pseudopotentials Virtual Vault},
  institution  = {Cornell NanoScale Science and Technology Facility},
  url          = {https://nninc.cnf.cornell.edu/}
}""")

_BIB_PSEUDO_BANK = {
    "dojo": _BIB_PSEUDO_DOJO,
    "virtual_vault": _BIB_VIRTUAL_VAULT,
}


def detect_pseudo_bank(pp_path):
    """If `pp_path` resolves to one of the bundled stb.core.pseudopotentials.
    BANKS folders, returns its bank name (so the matching citation can be
    added); None for a blank path or a user-supplied folder -- there's no
    reliable citation for pseudopotentials this suite didn't curate itself.
    """
    if not pp_path:
        return None
    normalized = os.path.normpath(pp_path)
    base = os.path.basename(normalized)
    parent = os.path.basename(os.path.dirname(normalized))
    if parent == "pseudopotentials" and base in BANKS:
        return base
    return None


def build_bib_entries(chosen_mode, d3_enabled, pseudo_bank=None):
    """Returns the (key, bibtex-entry) pairs for the methods this specific
    run's calc.fdf actually exercises, so the user has everything needed to
    cite the calculation without hunting down references by hand."""
    entries = [_BIB_SIESTA, _BIB_SIESTA_RECENT, _BIB_PBE]
    if d3_enabled:
        entries.append(_BIB_DFTD3)
    if chosen_mode in ('aimd', 'aimd+d3'):
        entries.append(_BIB_NOSE)
    if pseudo_bank in _BIB_PSEUDO_BANK:
        entries.append(_BIB_PSEUDO_BANK[pseudo_bank])
    return entries


def write_bib_file(bib_path, entries):
    """Writes a plain .bib file with one blank-line-separated entry per
    reference -- a standard BibTeX file, importable as-is into any reference
    manager (Zotero, JabRef, LaTeX's own bibtex/biber)."""
    with open(bib_path, "w") as f_bib:
        f_bib.write("% References for the SIESTA calculation generated by stb-inputfile\n")
        f_bib.write("% Cite these alongside the SIESTA code itself.\n\n")
        for _key, entry in entries:
            f_bib.write(entry + "\n\n")


def parse_structure_fdf(filename):
    """
    Parses a .fdf (Siesta) file and returns the lattice vectors, the list of
    chemical species symbols, and the per-axis vacuum flags (see
    kspace.detect_vacuum_axes) used to compute a dimensionality-aware k-grid.
    """
    structure = structure_io.read_fdf(filename)
    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, VACUUM_GAP_ANG)
    return structure.lattice, structure_io.species_list(structure), vacuum_axes

def copy_pseudopotentials(species_list, pp_path, f_out=None):
    """
    Copies the required .psml or .psf files from the 'pp_path' folder to
    the current working directory. Prioritizes .psml if both exist.
    Returns a list of warnings/errors.
    """
    warnings = []
    if not os.path.isdir(pp_path):
        warnings.append(f"[WARNING] PP path '{pp_path}' is not a valid directory -- pseudopotentials not copied.")
        return warnings

    copied_files = 0
    for symbol in species_list:
        psml_filename = f"{symbol}.psml"
        psf_filename = f"{symbol}.psf"

        source_psml = os.path.join(pp_path, psml_filename)
        source_psf = os.path.join(pp_path, psf_filename)

        dest_psml = os.path.join(os.getcwd(), psml_filename)
        dest_psf = os.path.join(os.getcwd(), psf_filename)

        # Verifica se o .psml existe (prioridade 1)
        if os.path.exists(source_psml):
            try:
                shutil.copy2(source_psml, dest_psml)
                print_dual(f"[OK] Copied: {psml_filename}", f_out)
                copied_files += 1
            except Exception as e:
                warnings.append(f"[ERROR] Failed to copy {psml_filename}: {e}")

        # Se não tiver .psml, tenta o .psf (prioridade 2)
        elif os.path.exists(source_psf):
            try:
                shutil.copy2(source_psf, dest_psf)
                print_dual(f"[OK] Copied: {psf_filename}", f_out)
                copied_files += 1
            except Exception as e:
                warnings.append(f"[ERROR] Failed to copy {psf_filename}: {e}")

        # Nenhum dos dois foi encontrado
        else:
            warnings.append(f"[WARNING] Neither {psml_filename} nor {psf_filename} found in {pp_path}")

    if copied_files == len(species_list):
        print_dual(color_text("All pseudopotentials copied successfully.", 'green'), f_out)
    else:
        warnings.append("[WARNING] Some PP files were not found or could not be copied.")

    return warnings

# --- Main Generation Logic Function ---
def generate_calculation(struct_file, chosen_mode, pp_path, f_out=None):
    """
    Main function that executes the file generation logic.
    """
    try:
        # --- 1. Input Validation ---
        if not os.path.exists(struct_file):
            print_dual(color_text(f"[ERROR] Structure file '{struct_file}' not found.", 'red'), f_out)
            return False  # Abort the function

        print_section("[0] RUN METADATA", f_out)
        print_dual(f"Structure file : {struct_file}", f_out)
        print_dual(f"Mode selected  : {chosen_mode}", f_out)
        if not pp_path:
            print_dual("Pseudopotentials : (none -- skipping pseudopotential copy)", f_out)
        else:
            print_dual(f"Pseudopotentials : {pp_path}", f_out)

        output_file = "calc.fdf"  # Fixed output name

        print_section("[1] TEMPLATE SELECTION", f_out)
        # --- Select Template ---
        template_string = ""
        if chosen_mode in ['relax', 'relax+d3']:
            template_string = CALC_RELAX_TEMPLATE
            print_dual("Template : Relax", f_out)
        elif chosen_mode in ['total_energy', 'total_energy+d3']:
            template_string = CALC_TOTAL_ENERGY_TEMPLATE
            print_dual("Template : Total Energy", f_out)
        elif chosen_mode in ['aimd', 'aimd+d3']:
            template_string = CALC_AIMD_TEMPLATE
            print_dual("Template : AIMD", f_out)
        elif chosen_mode in ['bands', 'bands+d3']:
            template_string = CALC_BANDS_TEMPLATE
            print_dual("Template : Bands", f_out)

        # --- Set D3 flag ---
        d3_flag = ".false."
        if chosen_mode in ['relax+d3', 'total_energy+d3', 'aimd+d3', 'bands+d3']:
            d3_flag = ".true."
            print_dual("DFT-D3 (van der Waals) correction : ENABLED", f_out)
        else:
            # 'relax', 'total_energy', 'aimd', or 'bands'
            print_dual("DFT-D3 (van der Waals) correction : DISABLED", f_out)

        d3_line_new = f"DFTD3                   {d3_flag}"
        template_lines = template_string.splitlines(keepends=True)

        # Parse structure and get species
        lattice, species, vacuum_axes = parse_structure_fdf(struct_file)
        print_dual(f"Species found : {', '.join(species)}", f_out)

        print_section("[2] K-GRID", f_out)
        # --- K-Grid Conditional Logic ---
        replace_kgrid = False
        kgrid_line_new = ""  # Initialize

        if chosen_mode in ['aimd', 'aimd+d3']:
            # For AIMD, do nothing, keep the K-grid from the template
            replace_kgrid = False
            print_dual("AIMD mode selected -- K-grid will NOT be modified.", f_out)
        else:
            # For Total Energy, Relax, and Bands, calculate the K-grid
            replace_kgrid = True
            print_dual("Calculating K-grid (density = 0.2 1/Ang) ...", f_out)
            k_density = 0.2
            kgrid_divs = kspace.compute_monkhorts(lattice[0], lattice[1], lattice[2], k_density, vacuum_axes)
            kgrid_line_new = f"kgrid.MonkhorstPack   [{kgrid_divs[0]}  {kgrid_divs[1]}  {kgrid_divs[2]}]"
            print_dual(f"Suggested K-grid : {kgrid_divs[0]} {kgrid_divs[1]} {kgrid_divs[2]}", f_out)

        # Use only the base file name for the include
        include_line_new = f"%include {os.path.basename(struct_file)}"

        print_section("[3] WRITING INPUT FILE", f_out)
        # Process the template and replace lines
        output_lines = []
        for line in template_lines:
            line_stripped_lower = line.strip().lower()

            # Replace the structure include line
            if line_stripped_lower.startswith('%include') and 'struct' in line_stripped_lower:
                output_lines.append(include_line_new + '\n')

            # Replace the k-grid (if not AIMD)
            elif line_stripped_lower.startswith('kgrid.monkhorstpack'):
                if replace_kgrid:
                    output_lines.append(kgrid_line_new + '\n')
                else:
                    # Keep the original template line (for AIMD)
                    output_lines.append(line)

            # Replace the D3 flag
            elif line_stripped_lower.startswith('dftd3'):
                output_lines.append(d3_line_new + '\n')

            else:
                # Keep the original template line
                output_lines.append(line)

        # Write the output file
        with open(output_file, 'w') as f_output:
            f_output.writelines(output_lines)

        print_dual(color_text(f"[OK] File '{output_file}' generated successfully.", 'green'), f_out)

        print_section("[4] REFERENCES", f_out)
        pseudo_bank = detect_pseudo_bank(pp_path)
        bib_entries = build_bib_entries(chosen_mode, d3_flag == ".true.", pseudo_bank)
        write_bib_file(BIB_FILE, bib_entries)
        if pseudo_bank:
            print_dual(f"Pseudopotential source : bundled bank '{pseudo_bank}'", f_out)
        print_dual(color_text(
            f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
            f"({len(bib_entries)} entries).", 'green'), f_out)

        # --- 5. Copy Pseudopotentials ---
        pp_warnings = []
        if pp_path:
            print_section("[5] PSEUDOPOTENTIALS", f_out)
            pp_warnings = copy_pseudopotentials(species, pp_path, f_out)

        if pp_warnings:
            print_section("[5b] PSEUDOPOTENTIAL WARNINGS", f_out)
            for warning in pp_warnings:
                print_dual(color_text(warning, 'yellow'), f_out)

        print_section("[6] SUMMARY & FILES", f_out)
        print_dual("Status       : OK", f_out)
        print_dual(f"Input file   : {output_file}", f_out)
        print_dual(f"References   : {BIB_FILE}", f_out)
        return True

    except Exception as e:
        # Catch any errors (e.g., File not found, Zero cell volume)
        print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
        print_dual(color_text("Operation aborted.", 'red'), f_out)
        return False

# --- Main Function (Argument Parser) ---
def main():
    """
    Processes the command-line arguments and calls the generation logic.
    """
    # --- Definition of valid modes ---
    mode_list = [
        'total_energy', 'total_energy+d3',
        'relax', 'relax+d3',
        'aimd', 'aimd+d3',
        'bands', 'bands+d3'
    ]

    # --- Configure argparse ---
    parser = argparse.ArgumentParser(
        description="Script to prepare a Siesta input file (.fdf).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s struct.fdf -t relax+d3 -p /path/to/pps\n"
               "  %(prog)s system.fdf -t bands\n"
               "  %(prog)s system.fdf --type aimd\n"
    )
    
    # Argument 1: Structure File (Required, Positional)
    parser.add_argument(
        "structure_file",
        type=str,
        help="Path to the structure file (e.g., struct.fdf)"
    )
    
    # Argument 2: Calculation Mode (Required, with Flag)
    parser.add_argument(
        "-t", "--type",
        type=str,
        required=True, # <-- The -t/--type flag is now mandatory
        choices=mode_list,
        dest="calc_type", # The value will be stored in 'args.calc_type'
        help=f"Calculation mode. Valid choices: {', '.join(mode_list)}"
    )
    
    # Argument 3: PPs Path (Optional, with flag)
    parser.add_argument(
        "-p", "--pp-path",
        type=str,
        default="",
        dest="pp_path",
        help=f"Pseudopotentials source (optional): a bundled bank ({', '.join(BANKS)}) "
             "or a folder path."
    )

    
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-inputfile {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.pp_path in BANKS:
        # Only resolve a recognized bank name -- an invalid value is left
        # as-is so copy_pseudopotentials()'s own os.path.isdir() check
        # produces its existing (non-fatal) "not a valid directory" warning,
        # same lenient behavior as before bundled banks existed.
        args.pp_path = resolve_pseudo_source(args.pp_path)

    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-INPUTFILE REPORT =====", 'magenta'), f_out)

    # Call the main logic function with the processed arguments
    ok = generate_calculation(args.structure_file, args.calc_type, args.pp_path, f_out)

    if report_path:
        print_dual(f"Report       : {report_path}", f_out)

    if f_out:
        f_out.close()

    if not ok:
        sys.exit(1)


# --- Script Entry Point ---
if __name__ == "__main__":
    main()
