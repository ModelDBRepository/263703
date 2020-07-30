c 28 October 2015, further speedup of NMDA onset, from alphaX6
c 25 October 2015: start alpha266.f
c Start by trying to speed up onset kinetics of NMDA EPSCs between suppyrRS
! 23 Dec 2014: start with spikewaveS5 (delta).
! For attempts to simulate L4 alpha in V1 in vitro (Miles and Karen)
! TO DO THIS: deep layers likely silent.  WE PRETEND THAT SUPPYRRS cells ARE L4 PYRAMIDS
! THAT BEHAVE LIKE scortalphaY674 and cell2 fig.pptx (IN ALPHA FOLDER), i.e. LOTS of gCaT
! SUPERFICIAL INTERNEURONS: WE PRETEND THESE ARE L4 INTERNEURONS

c 26 Feb. 2013, start with spikewaveR70, which Miles liked, alter synaptic inputs
c to nontuftRS as per spikewaveR73; then put integration routine for nontuftRS
c into main program, so as to be able to alter intrinsic conductances.

c 21 Feb. 2013, start with spikewave69.f, but incorporate suppyrRS integration
c routine into main file, so as to be able to alter intrinsic properties.

c 7 Dec. 2011, from deltaAIX108.f: for use with delta and nicotinic spike-wave,
c involving all types of cortical cells - see notebook figs from Miles, this date.
! 4 June 2009, from Buster's delta107.f - rework to run on 16 nodes.
! 27 Oct. 2006, add ability to modify Ih in tuftIB.
! 22 Oct. 2006, further modification to all GABA-B (if large enough in local
! population of cells) to shut off ectopics in tuftIB.
! 19 Oct. 2006, from interactB.f, has GABA-B, for use in delta oscillations.
c Last interactB data: INTERACTB30....
! 17 Sept. 2006: start with /interact/interact.f, and add neurogliaform cells:
! supng and deepng.  These cells make GABA-B IPSCs in cortical principal cells.
! In addition, nRT cells can make GABA-B IPSCs in each other and in TCR cells.
! Integration subroutines testable, and tidied up - see red lab notebook.

! 7 Nov. 2005: make suppyrRS & FRB axons Colbert/Pan.
! 24 Oct. 2005, take isoldeepVFOKM2.f and restore synaptic connections as in groucho84.f.
! To inactivate deep or superficial layers, try hyperpolarizing axons of respective cell
! populations.

! 22 Oct. 2005, start with isoldeepVFOKM2.f; use different tuftIB integration
! subroutine, allowing passage of a variable for shifting gKM rate functions.
! 31 Aug. 2005: in integration subroutine for tuftIB, VK relatively depolarized
! (in groucho it was -95 mV); AND scaling of gCaL depends on the cell.

! 29 July 2005, begin revising isoldeepVFO.f so that cells of a given type can
! be distributed across many nodes.

! Alter distal axon structure: assumption in groucho.f is that at most 1000 cells
! of any one cell type (27 July 2005); all runs of isoldeepVFO and isoldeepcut
! up to this point, with more than 1000 tuftIB cells, are therefore suspect.

! 14 June 2005, taken from isoldeep/isoldeep99.f; allows passing to tuftIB
! integration subroutine the relative shift of axonal gNa kinetics.
! 26 Jan. 2005: take Diego42.f & separate layers 5 & 6 from thalamus, and layers 
!  4 and 2/3; perhaps also block all NMDA in deep layers.
! reduce drive to nontuftRS
! pass parameters to regulate gK(AHP) & gNaP in tuftIB cells, and gKM and gKA
! 11 October 2004, modify diego.f and also integrate_nontuftRS.f in /diego so that
! nontuftRS_nFRB of the num_spinstell cells are FRB cells, as in layVtup.f.3Feb04.

! 15 August 2004: copy gutnick.f; use modified stellate cell as in /spinst/bryan.f;
! For study of spiny stellate response to thalamic stimulation, in full network.

! 25 Aug. 2004: necessary to explore TCR -> deep interneuron connectivity, so
! that deep interneurons fire very soon after TCR volley, e.g. make connections
! proximal and strong....

c 14 June 2004, version of groucho.f (from groucho94.f), with which, in this
c directory (gutnick) integrate_spinstell.f is modified: NMDA conductance is
c neither voltage nor Mg dependent.  c.f. Fleidervish et al. 1998 NEURON
c In addition, in PROGRAM GROUCHO, tau NMDA for these connections likely to
c be shortened.


c Cell types: superficial pyramidal RS; superficial pyramidal FRB;
c superficial basket; superficial axoaxonic; superficial LTS;
c spiny stellate; tufted deep IB; tufted deep RS; nontufted deep RS;
c deep basket; deep axoaxonic; deep LTS; thalamocortical relay (TCR);
c nucleus reticularis thalami (nRT).

               PROGRAM GROUCHO

! Reduce numbers of some non-deep cells, to speed things up
        PARAMETER (num_suppyrRS =  240, num_suppyrFRB =  10,
     & num_supbask = 90, num_supaxax = 90, num_supLTS = 90,
     & num_spinstell = 240, num_tuftIB =1000, num_tuftRS = 200, ! Diego unsure how many tuftRS there are
     & num_nontuftRS = 500, nontuftRS_nFRB =   2,
! Last nontuftRS_nFRB nontuftRS cells are FRB 
     & num_deepbask = 100, 
     & num_deepaxax = 100, num_deepLTS = 100, num_TCR = 100,
     & num_supng = 90, num_deepng = 100,
     & num_nRT = 100, num_allcells = 4950, num_celltypes = 16)

        PARAMETER (ncellspernode_suppyrRS =  60, 
     &    nodesfor_suppyrRS  = 4,
     &             ncellspernode_suppyrFRB=  10,
     &    nodesfor_suppyrFRB = 1,
     &             ncellspernode_supbask  =  90,
     &    nodesfor_supbask   = 1,
     &             ncellspernode_supaxax  =  90,
     &    nodesfor_supaxax   = 1,
     &             ncellspernode_supLTS   =  90,
     &    nodesfor_supLTS    = 1,
     &             ncellspernode_spinstell= 240,  
     &    nodesfor_spinstell = 1,
     &             ncellspernode_tuftIB   = 200,
     &    nodesfor_tuftIB    = 5,
     &             ncellspernode_tuftRS   = 200,
     &    nodesfor_tuftRS    = 1,
     &             ncellspernode_nontuftRS= 250,  
     &    nodesfor_nontuftRS = 2,
     &             ncellspernode_deepbask = 100,
     &    nodesfor_deepbask  = 1,
     &             ncellspernode_deepaxax = 100,
     &    nodesfor_deepaxax  = 1,
     &             ncellspernode_deepLTS  = 100, 
     &    nodesfor_deepLTS   = 1,
     &             ncellspernode_TCR      = 100,
     &    nodesfor_TCR       = 1,
     &             ncellspernode_nRT      = 100,            
     &    nodesfor_nRT       = 1,
     &             ncellspernode_supng    =  90,
     &    nodesfor_supng     = 1,
     &             ncellspernode_deepng   = 100,
     &    nodesfor_deepng    = 1)

c       PARAMETER (numnodes = 27)  ! Check manually for consistency.
        PARAMETER (numnodes = 24)  ! Check manually for consistency.
        PARAMETER (maxcellspernode = 250)
c       PARAMETER (maxcellspernode = 400)
            

        PARAMETER (numcomp_suppyrRS   = 74,
     &             numcomp_suppyrFRB  = 74,
     &             numcomp_supbask    = 59,
     &             numcomp_supaxax    = 59,
     &             numcomp_supLTS     = 59,
     &             numcomp_spinstell  = 59,
     &             numcomp_tuftIB     = 61,
     &             numcomp_tuftRS     = 61,
     &             numcomp_nontuftRS  = 50,
     &             numcomp_deepbask   = 59,
     &             numcomp_deepaxax   = 59,
     &             numcomp_deepLTS    = 59,
     &             numcomp_TCR        =137,
     &             numcomp_nRT        = 59,
     &             numcomp_supng      = 59,
     &             numcomp_deepng     = 59)

        PARAMETER (num_suppyrRS_to_suppyrRS = 50,
     &   num_suppyrRS_to_suppyrFRB = 50,
     &   num_suppyrRS_to_supbask   = 90, ! note
     &   num_suppyrRS_to_supaxax   = 90, ! note
     &   num_suppyrRS_to_supLTS    = 90, ! note
     &   num_suppyrRS_to_supng     = 90, ! note
     &   num_suppyrRS_to_spinstell = 10, 
     &   num_suppyrRS_to_tuftIB    =  1, 
     &   num_suppyrRS_to_tuftRS    =  1,
     &   num_suppyrRS_to_deepbask  = 10, 
! connections from sup pyrs to deep interneurons not known, although exist: try smaller value
     &   num_suppyrRS_to_deepaxax  = 10, ! see above
     &   num_suppyrRS_to_deepLTS   = 10, ! see above
! No suppyrRS to deepng

     &   num_suppyrRS_to_nontuftRS =  3, 
     &   num_suppyrFRB_to_suppyrRS = 1,
     &   num_suppyrFRB_to_suppyrFRB= 1,
     &   num_suppyrFRB_to_supbask  = 1,
     &   num_suppyrFRB_to_supaxax  = 1,
     &   num_suppyrFRB_to_supLTS   = 1,
     &   num_suppyrFRB_to_supng    = 1,
     &   num_suppyrFRB_to_spinstell= 1, 
     &   num_suppyrFRB_to_tuftIB   = 1,
     &   num_suppyrFRB_to_tuftRS   = 1,
     &   num_suppyrFRB_to_deepbask = 1, ! see note above with suppyrRS_to_deepbask
     &   num_suppyrFRB_to_deepaxax = 1, ! see note above
     &   num_suppyrFRB_to_deepLTS  = 1, ! see note above
! No suppyrFRB to deepng
     &   num_suppyrFRB_to_nontuftRS= 1) ! small per Thomson & Bannister
        PARAMETER
     &  (num_supbask_to_suppyrRS   = 20,
     &   num_supbask_to_suppyrFRB  = 20,
     &   num_supbask_to_supbask    = 20,
     &   num_supbask_to_supaxax    = 20,
     &   num_supbask_to_supLTS     = 20,
     &   num_supbask_to_supng      = 20,
     &   num_supbask_to_spinstell  = 20,

     &   num_supaxax_to_suppyrRS   = 20, ! note
     &   num_supaxax_to_suppyrFRB  = 20,  ! note
     &   num_supaxax_to_spinstell  = 5,
     &   num_supaxax_to_tuftIB     = 5,
     &   num_supaxax_to_tuftRS     = 5,
     &   num_supaxax_to_nontuftRS  = 5,

     &   num_supLTS_to_suppyrRS    = 20,
     &   num_supLTS_to_suppyrFRB   = 20,
     &   num_supLTS_to_supbask     = 20,
     &   num_supLTS_to_supaxax     = 20,
     &   num_supLTS_to_supLTS      = 20,
     &   num_supLTS_to_spinstell   = 20,
     &   num_supLTS_to_tuftIB      = 20)

        PARAMETER
     &  (num_supng_to_suppyrRS     = 20,
     &   num_supng_to_suppyrFRB    = 20,
     &   num_supng_to_nontuftRS    = 20,
     &   num_supng_to_tuftIB       = 20,
     &   num_supng_to_tuftRS       = 20,
     &   num_supng_to_supng        = 10,
     &   num_supng_to_supbask      = 10)

        PARAMETER
     &  (num_supLTS_to_tuftRS      = 20,
     &   num_supLTS_to_deepbask    = 20,
     &   num_supLTS_to_deepaxax    = 20,
     &   num_supLTS_to_deepLTS     = 20,
     &   num_supLTS_to_nontuftRS   = 20,

     &   num_spinstell_to_suppyrRS = 20,
     &   num_spinstell_to_suppyrFRB= 20,
     &   num_spinstell_to_supbask  = 20,
     &   num_spinstell_to_supaxax  = 20,
     &   num_spinstell_to_supLTS   = 20,
     &   num_spinstell_to_spinstell= 30,
     &   num_spinstell_to_tuftIB   = 20,
     &   num_spinstell_to_tuftRS   = 20,
     &   num_spinstell_to_deepbask = 20,
     &   num_spinstell_to_deepaxax = 20,
c    &   num_spinstell_to_deepLTS  = 20,
     &   num_spinstell_to_deepLTS  = 40,
     &   num_spinstell_to_deepng   = 20,
     &   num_spinstell_to_nontuftRS= 20,

     &   num_tuftIB_to_suppyrRS    =  2,
     &   num_tuftIB_to_suppyrFRB   =  2, ! small per Thomson & Bannister
     &   num_tuftIB_to_supbask     = 10) ! 17 Nov. 2004, try smaller value

        PARAMETER
     &  (num_tuftIB_to_supaxax     = 10, ! 17 Nov. 2004, try smaller value
     &   num_tuftIB_to_supLTS      = 10, ! 17 Nov. 2004, try smaller value
     &   num_tuftIB_to_spinstell   = 20,
     &   num_tuftIB_to_tuftIB      = 50,
     &   num_tuftIB_to_tuftRS      = 20,
     &   num_tuftIB_to_deepbask    = 20,
     &   num_tuftIB_to_deepaxax    = 20,
     &   num_tuftIB_to_deepLTS     = 20,
     &   num_tuftIB_to_deepng      = 50,
     &   num_tuftIB_to_nontuftRS   = 20,

     &   num_tuftRS_to_suppyrRS    =  2, ! small per Thomson & Bannister
     &   num_tuftRS_to_suppyrFRB   =  2, ! small per Thomson & Bannister
     &   num_tuftRS_to_supbask     = 20,
     &   num_tuftRS_to_supaxax     = 20,
     &   num_tuftRS_to_supLTS      = 20,
     &   num_tuftRS_to_spinstell   = 20,
     &   num_tuftRS_to_tuftIB      = 20,
     &   num_tuftRS_to_tuftRS      = 10,
     &   num_tuftRS_to_deepbask    = 20,
     &   num_tuftRS_to_deepaxax    = 20,
     &   num_tuftRS_to_deepLTS     = 20,
     &   num_tuftRS_to_deepng      = 20,
     &   num_tuftRS_to_nontuftRS   = 20)

        PARAMETER
     &  (num_deepbask_to_spinstell = 20,
     &   num_deepbask_to_tuftIB    = 20,
     &   num_deepbask_to_tuftRS    = 20,
     &   num_deepbask_to_deepbask  = 20,
     &   num_deepbask_to_deepaxax  = 20,
     &   num_deepbask_to_deepLTS   = 20,
     &   num_deepbask_to_deepng    = 20,
     &   num_deepbask_to_nontuftRS = 20,

     &   num_deepaxax_to_suppyrRS  =  5,
     &   num_deepaxax_to_suppyrFRB =  5,
     &   num_deepaxax_to_spinstell =  5,
     &   num_deepaxax_to_tuftIB    =  5,
     &   num_deepaxax_to_tuftRS    =  5,
     &   num_deepaxax_to_nontuftRS =  5,

     &   num_deepLTS_to_suppyrRS   = 10)
        PARAMETER
     &  (num_deepLTS_to_suppyrFRB  = 10,
     &   num_deepLTS_to_supbask    = 10,
     &   num_deepLTS_to_supaxax    = 10,
     &   num_deepLTS_to_supLTS     = 10,
c    &   num_deepLTS_to_spinstell  = 20,
     &   num_deepLTS_to_spinstell  = 40,
     &   num_deepLTS_to_tuftIB     = 20,
     &   num_deepLTS_to_tuftRS     = 20,
     &   num_deepLTS_to_deepbask   = 20,
     &   num_deepLTS_to_deepaxax   = 20,
     &   num_deepLTS_to_deepLTS    = 20,
     &   num_deepLTS_to_nontuftRS  = 20)

        PARAMETER
     &  (num_deepng_to_tuftIB      = 20,
     &   num_deepng_to_tuftRS      = 20,
     &   num_deepng_to_nontuftRS   = 20,
     &   num_deepng_to_spinstell   = 20,
     &   num_deepng_to_deepng      = 10,
     &   num_deepng_to_deepbask    = 10,

     &   num_TCR_to_suppyrRS       = 10,
     &   num_TCR_to_suppyrFRB      = 10,
     &   num_TCR_to_supbask        = 10,
     &   num_TCR_to_supaxax        = 10,
     &   num_TCR_to_supng          = 10,
     &   num_TCR_to_spinstell      = 20,
     &   num_TCR_to_tuftIB         = 10,
     &   num_TCR_to_tuftRS         = 10,
     &   num_TCR_to_deepbask       = 20,
     &   num_TCR_to_deepaxax       = 10,
     &   num_TCR_to_deepng         = 10,
     &   num_TCR_to_nRT            = 25, ! note
     &   num_TCR_to_nontuftRS      = 10,

     &   num_nRT_to_TCR            = 15, ! note
     &   num_nRT_to_nRT            = 10)
        PARAMETER
     &  (num_nontuftRS_to_suppyrRS = 10,
     &   num_nontuftRS_to_suppyrFRB= 10,
     &   num_nontuftRS_to_supbask  = 10,
     &   num_nontuftRS_to_supaxax  = 10,
     &   num_nontuftRS_to_supLTS   = 10,
     &   num_nontuftRS_to_spinstell= 10,
     &   num_nontuftRS_to_tuftIB   = 10,
     &   num_nontuftRS_to_tuftRS   = 10,
     &   num_nontuftRS_to_deepbask = 10,
     &   num_nontuftRS_to_deepaxax = 10,
     &   num_nontuftRS_to_deepLTS  = 10,
     &   num_nontuftRS_to_deepng   = 10,
     &   num_nontuftRS_to_TCR      = 20,
     &   num_nontuftRS_to_nRT      = 20,
     &   num_nontuftRS_to_nontuftRS= 20)

c Begin definition of number of compartments that can be
c contacted for each type of synaptic connection.
c       PARAMETER (ncompallow_suppyrRS_to_suppyrRS = 36, ! original, for connections to basal and oblique
        PARAMETER (ncompallow_suppyrRS_to_suppyrRS = 28, ! for connections to distal apical
     &   ncompallow_suppyrRS_to_suppyrFRB = 36,
     &   ncompallow_suppyrRS_to_supbask   = 24,
     &   ncompallow_suppyrRS_to_supaxax   = 24,
     &   ncompallow_suppyrRS_to_supLTS    = 24,
     &   ncompallow_suppyrRS_to_supng     = 52,
     &   ncompallow_suppyrRS_to_spinstell = 24,
     &   ncompallow_suppyrRS_to_tuftIB    =  8,
     &   ncompallow_suppyrRS_to_tuftRS    =  8,
     &   ncompallow_suppyrRS_to_deepbask  = 24,
     &   ncompallow_suppyrRS_to_deepaxax  = 24,
     &   ncompallow_suppyrRS_to_deepLTS   = 24,
     &   ncompallow_suppyrRS_to_nontuftRS =  7)

        PARAMETER (ncompallow_suppyrFRB_to_suppyrRS = 36,
     &   ncompallow_suppyrFRB_to_suppyrFRB = 36,
     &   ncompallow_suppyrFRB_to_supbask   = 24,
     &   ncompallow_suppyrFRB_to_supaxax   = 24,
     &   ncompallow_suppyrFRB_to_supLTS    = 24,
     &   ncompallow_suppyrFRB_to_supng     = 52,
     &   ncompallow_suppyrFRB_to_spinstell = 24,
     &   ncompallow_suppyrFRB_to_tuftIB    =  8,
     &   ncompallow_suppyrFRB_to_tuftRS    =  8,
     &   ncompallow_suppyrFRB_to_deepbask  = 24,
     &   ncompallow_suppyrFRB_to_deepaxax  = 24,
     &   ncompallow_suppyrFRB_to_deepLTS   = 24,
     &   ncompallow_suppyrFRB_to_nontuftRS =  7)

        PARAMETER (ncompallow_supbask_to_suppyrRS   = 11,
     &   ncompallow_supbask_to_suppyrFRB   = 11,
     &   ncompallow_supbask_to_supbask     = 24,
     &   ncompallow_supbask_to_supng       = 4, 
     &   ncompallow_supbask_to_supaxax     = 24,
     &   ncompallow_supbask_to_supLTS      = 24,
     &   ncompallow_supbask_to_spinstell   =  5)

        PARAMETER (ncompallow_supLTS_to_suppyrRS    = 53,
     &   ncompallow_supLTS_to_suppyrFRB    = 53,
     &   ncompallow_supLTS_to_supbask      = 40,
     &   ncompallow_supLTS_to_supaxax      = 40,
     &   ncompallow_supLTS_to_supLTS       = 40,
     &   ncompallow_supLTS_to_spinstell    = 40,
     &   ncompallow_supLTS_to_tuftIB       = 40,
     &   ncompallow_supLTS_to_tuftRS       = 40,
     &   ncompallow_supLTS_to_deepbask     = 20,
     &   ncompallow_supLTS_to_deepaxax     = 20,
     &   ncompallow_supLTS_to_deepLTS      = 20,
     &   ncompallow_supLTS_to_nontuftRS    = 29)

        PARAMETER (ncompallow_supng_to_suppyrRS = 64,
     &   ncompallow_supng_to_suppyrFRB     = 64,
     &   ncompallow_supng_to_nontuftRS     =  5,
     &   ncompallow_supng_to_tuftIB        = 14,
     &   ncompallow_supng_to_tuftRS        = 14,
     &   ncompallow_supng_to_supng         =  4,
     &   ncompallow_supng_to_supbask       =  4)

        PARAMETER (ncompallow_spinstell_to_suppyrRS = 24,
     &   ncompallow_spinstell_to_suppyrFRB = 24,
     &   ncompallow_spinstell_to_supbask   = 24,
     &   ncompallow_spinstell_to_supaxax   = 24,
     &   ncompallow_spinstell_to_supLTS    = 24,
     &   ncompallow_spinstell_to_spinstell = 24,
     &   ncompallow_spinstell_to_tuftIB    = 12,
     &   ncompallow_spinstell_to_tuftRS    = 12,
     &   ncompallow_spinstell_to_deepbask  = 24,
     &   ncompallow_spinstell_to_deepaxax  = 24,
     &   ncompallow_spinstell_to_deepLTS   = 24,
     &   ncompallow_spinstell_to_deepng    = 52,
     &   ncompallow_spinstell_to_nontuftRS =  5)

        PARAMETER (ncompallow_tuftIB_to_suppyrRS   = 13,
     &   ncompallow_tuftIB_to_suppyrFRB    = 13,
     &   ncompallow_tuftIB_to_supbask      = 24,
     &   ncompallow_tuftIB_to_supaxax      = 24,
     &   ncompallow_tuftIB_to_supLTS       = 24,
     &   ncompallow_tuftIB_to_spinstell    = 24,
     &   ncompallow_tuftIB_to_tuftIB       = 46,
     &   ncompallow_tuftIB_to_tuftRS       = 46,
     &   ncompallow_tuftIB_to_deepbask     = 24,
     &   ncompallow_tuftIB_to_deepaxax     = 24,
     &   ncompallow_tuftIB_to_deepLTS      = 24,
     &   ncompallow_tuftIB_to_deepng       = 52,
     &   ncompallow_tuftIB_to_nontuftRS    = 43)

        PARAMETER (ncompallow_tuftRS_to_suppyrRS   = 13,
     &   ncompallow_tuftRS_to_suppyrFRB    = 13,
     &   ncompallow_tuftRS_to_supbask      = 24,
     &   ncompallow_tuftRS_to_supaxax      = 24,
     &   ncompallow_tuftRS_to_supLTS       = 24,
     &   ncompallow_tuftRS_to_spinstell    = 24,
     &   ncompallow_tuftRS_to_tuftIB       = 46,
     &   ncompallow_tuftRS_to_tuftRS       = 46,
     &   ncompallow_tuftRS_to_deepbask     = 24,
     &   ncompallow_tuftRS_to_deepaxax     = 24,
     &   ncompallow_tuftRS_to_deepLTS      = 24,
     &   ncompallow_tuftRS_to_deepng       = 52,
     &   ncompallow_tuftRS_to_nontuftRS    = 43)

        PARAMETER (ncompallow_deepbask_to_spinstell = 5, 
     &   ncompallow_deepbask_to_tuftIB     =  8,
     &   ncompallow_deepbask_to_tuftRS     =  8,
     &   ncompallow_deepbask_to_deepbask   = 24,
     &   ncompallow_deepbask_to_deepaxax   = 24,
     &   ncompallow_deepbask_to_deepLTS    = 24,
     &   ncompallow_deepbask_to_deepng     =  4,
     &   ncompallow_deepbask_to_nontuftRS  =  8)

        PARAMETER (ncompallow_deepLTS_to_suppyrRS = 53,
     &   ncompallow_deepLTS_to_suppyrFRB   = 53,
     &   ncompallow_deepLTS_to_supbask     = 20,
     &   ncompallow_deepLTS_to_supaxax     = 20,
     &   ncompallow_deepLTS_to_supLTS      = 20,
     &   ncompallow_deepLTS_to_spinstell   = 40,
     &   ncompallow_deepLTS_to_tuftIB      = 40,
     &   ncompallow_deepLTS_to_tuftRS      = 40,
     &   ncompallow_deepLTS_to_deepbask    = 40,
     &   ncompallow_deepLTS_to_deepaxax    = 40,
     &   ncompallow_deepLTS_to_deepLTS     = 40,
     &   ncompallow_deepLTS_to_nontuftRS   = 29)

        PARAMETER (ncompallow_deepng_to_tuftIB = 33,
     &   ncompallow_deepng_to_tuftRS    = 33,
     &   ncompallow_deepng_to_nontuftRS = 33,
     &   ncompallow_deepng_to_spinstell = 52,
     &   ncompallow_deepng_to_deepng    =  4,
     &   ncompallow_deepng_to_deepbask  = 52)

        PARAMETER (ncompallow_TCR_to_suppyrRS = 24,
     &   ncompallow_TCR_to_suppyrFRB    = 24,
     &   ncompallow_TCR_to_supbask      = 12,
     &   ncompallow_TCR_to_supaxax      = 12,
     &   ncompallow_TCR_to_supng        = 52,
     &   ncompallow_TCR_to_spinstell    = 52,
     &   ncompallow_TCR_to_tuftIB       =  9,
     &   ncompallow_TCR_to_tuftRS       =  9,
!    &   ncompallow_TCR_to_deepbask     = 12,
     &   ncompallow_TCR_to_deepbask     =  5, ! make them proximal
!    &   ncompallow_TCR_to_deepaxax     = 12,
     &   ncompallow_TCR_to_deepaxax     =  5, ! make them proximal
     &   ncompallow_TCR_to_deepng       = 52, 
     &   ncompallow_TCR_to_nRT          = 12,
     &   ncompallow_TCR_to_nontuftRS    =  5)

        PARAMETER (ncompallow_nRT_to_TCR  = 11,
     &   ncompallow_nRT_to_nRT = 53)

        PARAMETER (ncompallow_nontuftRS_to_suppyrRS = 4,
     &    ncompallow_nontuftRS_to_suppyrFRB =  4,
     &    ncompallow_nontuftRS_to_supbask   = 24,
     &    ncompallow_nontuftRS_to_supaxax   = 24,
     &    ncompallow_nontuftRS_to_supLTS    = 24,
     &    ncompallow_nontuftRS_to_spinstell = 24,
     &    ncompallow_nontuftRS_to_tuftIB    = 46,
     &    ncompallow_nontuftRS_to_tuftRS    = 46,
     &    ncompallow_nontuftRS_to_deepbask  = 24,
     &    ncompallow_nontuftRS_to_deepaxax  = 24,
     &    ncompallow_nontuftRS_to_deepLTS   = 24,
     &    ncompallow_nontuftRS_to_deepng    = 52,
     &    ncompallow_nontuftRS_to_TCR       = 90,
     &    ncompallow_nontuftRS_to_nRT       = 12,
     &    ncompallow_nontuftRS_to_nontuftRS = 43)
c End definition of number of allowed compartments that
c can be contacted for each sort of connection

c Note that gj form only between cells of a given type,
c  except suppyrRS/suppyrFRB & tuftIB/tuftRS
c gj/cell = 2 x total gj / # cells
c for proportions, see /home/traub/supergj/tests.f
       integer, parameter :: totaxgj_suppyrRS =  150
!      integer, parameter :: totaxgj_suppyrRS =  270
c      integer, parameter :: totaxgj_suppyrRS =  225
c      integer, parameter :: totaxgj_suppyrRS =  195
       integer, parameter :: totaxgj_suppyrFRB =   2 
       integer, parameter :: totaxgj_suppyr    =   1 
! totaxgj_suppyr = number of "mixed" gj between RS suppyr
! (1st col. of table) and FRB suppyr (3rd col. of table) cells
       integer, parameter :: totSDgj_supbask   = 200 
       integer, parameter :: totSDgj_supaxax   =   0 
       integer, parameter :: totSDgj_supLTS    = 200 
!      integer, parameter :: totaxgj_spinstell = 240 
       integer, parameter :: totaxgj_spinstell = 250 ! keep small: axons now very excitable 
       integer, parameter :: totaxgj_tuftIB    =1000 
       integer, parameter :: totaxgj_tuftRS    = 160 
       integer, parameter :: totaxgj_tuft      =   1 ! decr. antidr. bursts in IB 
! totaxgj_tuft for mixed gj between tuftIB (1st) and tuftRS (next)
       integer, parameter :: totaxgj_nontuftRS = 500 
       integer, parameter :: totSDgj_deepbask  = 250 
       integer, parameter :: totSDgj_deepaxax  =   0 
       integer, parameter :: totSDgj_deepLTS   = 250 
       integer, parameter :: totaxgj_TCR       = 100  
       integer, parameter :: totSDgj_nRT       = 250
       integer, parameter :: totSDgj_supng     = 250
       integer, parameter :: totSDgj_deepng    = 250
c Note: no gj between axoaxonic cells.

c Define number of compartments on a cell where a gj might form
       integer, parameter :: num_axgjcompallow_suppyrRS = 1
       integer, parameter :: num_axgjcompallow_suppyrFRB= 1
       integer, parameter :: num_SDgjcompallow_supbask  = 8
       integer, parameter :: num_SDgjcompallow_supng    = 8
       integer, parameter :: num_SDgjcompallow_supLTS   = 8
       integer, parameter :: num_axgjcompallow_spinstell= 1
       integer, parameter :: num_axgjcompallow_tuftIB   = 1
       integer, parameter :: num_axgjcompallow_tuftRS   = 1
       integer, parameter :: num_axgjcompallow_nontuftRS= 1
       integer, parameter :: num_SDgjcompallow_deepbask = 8
       integer, parameter :: num_SDgjcompallow_deepng   = 8
       integer, parameter :: num_SDgjcompallow_deepLTS  = 8
       integer, parameter :: num_axgjcompallow_TCR      = 1
       integer, parameter :: num_SDgjcompallow_nRT      = 8

c Define gap junction conductances.
c      double precision, parameter :: gapcon_suppyrRS  = 4.0d-3
       double precision, parameter :: gapcon_suppyrRS  = 0.d-3 ! to see if superf. lay. can follow 40 Hz
       double precision, parameter :: gapcon_suppyrFRB = 1.0d-3
!      double precision, parameter :: gapcon_suppyrFRB = 0.d-3 ! to see if superf. lay. can follow 40 Hz
       double precision, parameter :: gapcon_supbask   = 1.d-3
       double precision, parameter :: gapcon_supng     = 0.5d-3
       double precision, parameter :: gapcon_supaxax   = 0.d-3
       double precision, parameter :: gapcon_supLTS    = 1.d-3

c      double precision, parameter :: gapcon_spinstell = 3.d-3
       double precision, parameter :: gapcon_spinstell = 0.d-3 

       double precision, parameter :: gapcon_tuftIB    = 4.d-3
!      double precision, parameter :: gapcon_tuftIB    = 0.d-3 ! to decr. antidr. bursting
       double precision, parameter :: gapcon_tuftRS    = 1.d-3
!      double precision, parameter :: gapcon_tuftRS    = 0.d-3 ! now follow 40 Hz?
c      double precision, parameter :: gapcon_nontuftRS = 4.d-3
       double precision, parameter :: gapcon_nontuftRS = 0.d-3 ! to abolish VFO in lay. 6
       double precision, parameter :: gapcon_deepbask  = 1.d-3
!      double precision, parameter :: gapcon_deepbask  = 0.d-3
       double precision, parameter :: gapcon_deepng    = 0.5d-3
       double precision, parameter :: gapcon_deepaxax  = 0.d-3
       double precision, parameter :: gapcon_deepLTS   = 1.d-3
!      double precision, parameter :: gapcon_deepLTS   = 0.d-3
!      double precision, parameter :: gapcon_TCR       = 3.d-3
       double precision, parameter :: gapcon_TCR       = 0.d-3
c      double precision, parameter :: gapcon_nRT       = 1.d-3
       double precision, parameter :: gapcon_nRT       = 0.4d-3

! Parameters for scaling tuftIB conductances
       double precision, parameter :: scale_tuftIB_gAR  =1.0d0
c      double precision, parameter :: scale_tuftIB_gKAHP=0.50d0 ! this value in isoldeep99.f
c      double precision, parameter :: scale_tuftIB_gKAHP=1.00d0 ! this value isoldeep80.f
       double precision, parameter :: scale_tuftIB_gKAHP=0.10d0 ! used in delta32
c      double precision, parameter :: scale_tuftIB_gKAHP=0.00d0 ! used in delta32
c      double precision, parameter :: scale_tuftIB_gKAHP=2.00d0
c      double precision, parameter :: scale_tuftIB_gNaP =1.0d0 ! 25 Mar. 2005, define below
c      double precision, parameter :: scale_tuftIB_gKM  =1.0d0 ! 10 Feb. 2005, this one defined below
       double precision, parameter :: scale_tuftIB_gKA  =1.0d0
c      double precision, parameter :: scale_tuftIB_gKC  =0.1d0 ! this value in isoldeep99.f
       double precision, parameter :: scale_tuftIB_gKC  =1.0d0 ! this value in isoldeep80.f
c      double precision, parameter :: scale_tuftIB_gCaL =0.7d0 ! define below



c Assorted parameters
         double precision, parameter :: dt = 0.002d0
!        double precision, parameter :: Mg = 0.10  ! for relatively Mg-independent NMDA in L4
         double precision, parameter :: Mg = 0.50 ! for NMDA-dependent CCh delta, try lower Mg
! Castro-Alamancos J Physiol, disinhib. neocortex in vitro, uses
! Mg = 1.3
         double precision, parameter :: NMDA_saturation_fact
!    &                                   = 5.d0
     &                                   = 80.d0
c NMDA conductance developed on one postsynaptic compartment,
c from one type of presynaptic cell, can be at most this
c factor x unitary conductance
c UNFORTUNATELY, with this scheme,if one NMDA cond. set to 0
c on a cell type, all NMDA conductances will be forced to 0
c on that cell type...

       double precision, parameter :: thal_cort_delay = 1.d0
       double precision, parameter :: cort_thal_delay = 5.d0
       integer, parameter :: how_often = 50
! how_often defines how many time steps between synaptic conductance
! updates, and between broadcastings of axonal voltages.
       double precision, parameter :: axon_refrac_time = 1.5d0

c For these ectopic rate parameters, assume noisepe checked
c every 200 time steps = 0.4 ms = 1./2.5 ms
      double precision, parameter :: noisepe_suppyrRS   =
     &      0.d0 
!    &            1.d0 / (2.5d0 * 1000.d0) ! USUAL
c    &            1.d0 / (2.5d0 *   400.d0) 
      double precision, parameter :: noisepe_suppyrFRB  =
c    &            1.d0 / (2.5d0 * 1000.d0) ! USUAL
     &            1.d0 / (2.5d0 * 10000.d0)
c    &      0.d0 
      double precision, parameter :: noisepe_spinstell  =
c    &            1.d0 / (2.5d0 * 10000.d0)
     &            1.d0 / (2.5d0 *  1000.d0) ! is this actually used?
! Note that noisepe_tuftIB will be time-dependent, and sensitive
! to total average GABA-B conductance in cells on each node
      double precision, parameter :: noisepe_tuftIB_save     =
c    &            1.d0 / (2.5d0 * 5000.d0)
c    &            1.d0 / (2.5d0 * 1000.d0)
     &            1.d0 / (2.5d0 * 10000.d0)
c    &            1.d0 / (2.5d0 *  250.d0) ! 29 March 2005
      double precision, parameter :: noisepe_tuftRS_save =
c this one also will be time_dependent
c    &            1.d0 / (2.5d0 * 1000.d0)
     &            1.d0 / (2.5d0 * 10000.d0)
      double precision, parameter :: noisepe_nontuftRS  =
c    &            1.d0 / (2.5d0 * 1000.d0)
     &            0.d0 / (2.5d0 * 1000.d0)
      double precision, parameter :: noisepe_TCR        =
     &            1.d0 / (2.5d0 * 1000.d0)


c Synaptic conductance time constants. 
      real*8, parameter :: tauAMPA_suppyrRS_to_suppyrRS=2.d0 
c     real*8, parameter :: tauNMDA_suppyrRS_to_suppyrRS=130.5d0 
c     real*8, parameter :: tauNMDA_suppyrRS_to_suppyrRS= 60.0d0 
c     real*8, parameter :: tauNMDA_suppyrRS_to_suppyrRS= 45.0d0 
      real*8, parameter :: tauNMDA_suppyrRS_to_suppyrRS= 15.0d0 
      real*8, parameter :: tauAMPA_suppyrRS_to_suppyrFRB=2.d0 
      real*8, parameter :: tauNMDA_suppyrRS_to_suppyrFRB=130.d0  
      real*8, parameter :: tauAMPA_suppyrRS_to_supbask  =.8d0   
      real*8, parameter :: tauNMDA_suppyrRS_to_supbask  =100.d0 
      real*8, parameter :: tauAMPA_suppyrRS_to_supng    =.8d0   
      real*8, parameter :: tauNMDA_suppyrRS_to_supng    =100.d0 
      real*8, parameter :: tauAMPA_suppyrRS_to_supaxax  =.8d0  
      real*8, parameter :: tauNMDA_suppyrRS_to_supaxax  =100.d0 
      real*8, parameter :: tauAMPA_suppyrRS_to_supLTS   =1.d0  
      real*8, parameter :: tauNMDA_suppyrRS_to_supLTS   =100.d0 
      real*8, parameter :: tauAMPA_suppyrRS_to_spinstell=2.d0   
      real*8, parameter :: tauNMDA_suppyrRS_to_spinstell=130.d0 
      real*8, parameter :: tauAMPA_suppyrRS_to_tuftIB   =2.d0   
      real*8, parameter :: tauNMDA_suppyrRS_to_tuftIB   =130.d0 
      real*8, parameter :: tauAMPA_suppyrRS_to_tuftRS   =2.d0   
      real*8, parameter :: tauNMDA_suppyrRS_to_tuftRS   =130.d0 
      real*8, parameter :: tauAMPA_suppyrRS_to_deepbask =.8d0   
      real*8, parameter :: tauNMDA_suppyrRS_to_deepbask =100.d0 
      real*8, parameter :: tauAMPA_suppyrRS_to_deepaxax =.8d0   
      real*8, parameter :: tauNMDA_suppyrRS_to_deepaxax =100.d0 
      real*8, parameter :: tauAMPA_suppyrRS_to_deepLTS  =1.d0   
      real*8, parameter :: tauNMDA_suppyrRS_to_deepLTS  =100.d0 
      real*8, parameter :: tauAMPA_suppyrRS_to_nontuftRS=2.d0   
      real*8, parameter :: tauNMDA_suppyrRS_to_nontuftRS=130.d0 

      real*8, parameter :: tauAMPA_suppyrFRB_to_suppyrRS=2.d0   
      real*8, parameter :: tauNMDA_suppyrFRB_to_suppyrRS=130.d0  
      real*8, parameter :: tauAMPA_suppyrFRB_to_suppyrFRB=2.d0   
      real*8, parameter :: tauNMDA_suppyrFRB_to_suppyrFRB=130.d0 
      real*8, parameter :: tauAMPA_suppyrFRB_to_supbask  =.8d0   
      real*8, parameter :: tauNMDA_suppyrFRB_to_supbask  =100.d0 
      real*8, parameter :: tauAMPA_suppyrFRB_to_supng    =.8d0   
      real*8, parameter :: tauNMDA_suppyrFRB_to_supng    =100.d0 
      real*8, parameter :: tauAMPA_suppyrFRB_to_supaxax  =.8d0  
      real*8, parameter :: tauNMDA_suppyrFRB_to_supaxax  =100.d0 
      real*8, parameter :: tauAMPA_suppyrFRB_to_supLTS   =1.d0  
      real*8, parameter :: tauNMDA_suppyrFRB_to_supLTS   =100.d0 
      real*8, parameter :: tauAMPA_suppyrFRB_to_spinstell=2.d0   
      real*8, parameter :: tauNMDA_suppyrFRB_to_spinstell=130.d0 
      real*8, parameter :: tauAMPA_suppyrFRB_to_tuftIB   =2.d0   
      real*8, parameter :: tauNMDA_suppyrFRB_to_tuftIB   =130.d0 
      real*8, parameter :: tauAMPA_suppyrFRB_to_tuftRS   =2.d0   
      real*8, parameter :: tauNMDA_suppyrFRB_to_tuftRS   =130.d0
      real*8, parameter :: tauAMPA_suppyrFRB_to_deepbask =.8d0   
      real*8, parameter :: tauNMDA_suppyrFRB_to_deepbask =100.d0 
      real*8, parameter :: tauAMPA_suppyrFRB_to_deepaxax =.8d0   
      real*8, parameter :: tauNMDA_suppyrFRB_to_deepaxax =100.d0 
      real*8, parameter :: tauAMPA_suppyrFRB_to_deepLTS  =1.d0   
      real*8, parameter :: tauNMDA_suppyrFRB_to_deepLTS  =100.d0 
      real*8, parameter :: tauAMPA_suppyrFRB_to_nontuftRS=2.d0   
      real*8, parameter :: tauNMDA_suppyrFRB_to_nontuftRS=130.d0 

      real*8,  parameter :: tauGABA_supbask_to_suppyrRS   =6.d0 ! 29 Nov. 2005, reduce from 6 to 5, also below. 
      real*8,  parameter :: tauGABA_supbask_to_suppyrFRB  =6.d0
      real*8,  parameter :: tauGABA_supbask_to_supbask    =3.d0  
      real*8,  parameter :: tauGABA_supbask_to_supng      =3.d0  
      real*8,  parameter :: tauGABA_supbask_to_supaxax    =3.d0  
      real*8,  parameter :: tauGABA_supbask_to_supLTS     =3.d0  
      real*8,  parameter :: tauGABA_supbask_to_spinstell  =6.d0 

      real*8,  parameter :: tauGABA_supaxax_to_suppyrRS   =6.d0 
      real*8,  parameter :: tauGABA_supaxax_to_suppyrFRB  =6.d0 
      real*8,  parameter :: tauGABA_supaxax_to_spinstell  =6.d0 
      real*8,  parameter :: tauGABA_supaxax_to_tuftIB     =6.d0 
      real*8,  parameter :: tauGABA_supaxax_to_tuftRS     =6.d0 
      real*8,  parameter :: tauGABA_supaxax_to_nontuftRS  =6.d0 

      real*8, parameter :: tauGABA_supLTS_to_suppyrRS    =20.d0 
      real*8, parameter :: tauGABA_supLTS_to_suppyrFRB   =20.d0 
      real*8, parameter :: tauGABA_supLTS_to_supbask     =20.d0 
      real*8, parameter :: tauGABA_supLTS_to_supaxax     =20.d0 
      real*8, parameter :: tauGABA_supLTS_to_supLTS      =20.d0 
      real*8, parameter :: tauGABA_supLTS_to_spinstell   =20.d0 
      real*8, parameter :: tauGABA_supLTS_to_tuftIB      =20.d0 
      real*8, parameter :: tauGABA_supLTS_to_tuftRS      =20.d0 
      real*8, parameter :: tauGABA_supLTS_to_deepbask    =20.d0 
      real*8, parameter :: tauGABA_supLTS_to_deepaxax    =20.d0 
      real*8, parameter :: tauGABA_supLTS_to_deepLTS     =20.d0 
      real*8, parameter :: tauGABA_supLTS_to_nontuftRS   =20.d0  

      real*8, parameter:: tauGABA_supng_to_suppyrRS      =6.d0
      real*8, parameter:: tauGABA_supng_to_suppyrFRB     =6.d0
      real*8, parameter:: tauGABA_supng_to_nontuftRS     =6.d0
      real*8, parameter:: tauGABA_supng_to_tuftIB        =6.d0
      real*8, parameter:: tauGABA_supng_to_tuftRS        =6.d0
      real*8, parameter:: tauGABA_supng_to_supng         =3.d0
      real*8, parameter:: tauGABA_supng_to_supbask       =3.d0

      real*8, parameter :: tauAMPA_spinstell_to_suppyrRS =2.d0  
      real*8, parameter :: tauNMDA_spinstell_to_suppyrRS =130.d0 
      real*8, parameter :: tauAMPA_spinstell_to_suppyrFRB=2.d0  
      real*8, parameter :: tauNMDA_spinstell_to_suppyrFRB=130.d0 
      real*8, parameter :: tauAMPA_spinstell_to_supbask  =.8d0  
      real*8, parameter :: tauNMDA_spinstell_to_supbask  =100.d0
      real*8, parameter :: tauAMPA_spinstell_to_supaxax  =.8d0  
      real*8, parameter :: tauNMDA_spinstell_to_supaxax  =100.d0
      real*8, parameter :: tauAMPA_spinstell_to_supLTS   =1.d0  
      real*8, parameter :: tauNMDA_spinstell_to_supLTS   =100.d0
      real*8, parameter :: tauAMPA_spinstell_to_spinstell=2.d0  
!     real*8, parameter :: tauNMDA_spinstell_to_spinstell=130.d0 
      real*8, parameter :: tauNMDA_spinstell_to_spinstell= 15.d0 ! small tau per Fleidervish et al., NEURON 
      real*8, parameter :: tauAMPA_spinstell_to_tuftIB   =2.d0  
      real*8, parameter :: tauNMDA_spinstell_to_tuftIB   =130.d0 
      real*8, parameter :: tauAMPA_spinstell_to_tuftRS   =2.d0  
      real*8, parameter :: tauNMDA_spinstell_to_tuftRS   =130.d0
      real*8, parameter :: tauAMPA_spinstell_to_deepbask =.8d0  
      real*8, parameter :: tauNMDA_spinstell_to_deepbask =100.d0
      real*8, parameter :: tauAMPA_spinstell_to_deepng   =.8d0  
      real*8, parameter :: tauNMDA_spinstell_to_deepng   =100.d0
      real*8, parameter :: tauAMPA_spinstell_to_deepaxax =.8d0  
      real*8, parameter :: tauNMDA_spinstell_to_deepaxax =100.d0
      real*8, parameter :: tauAMPA_spinstell_to_deepLTS  =1.d0  
      real*8, parameter :: tauNMDA_spinstell_to_deepLTS  =100.d0
      real*8, parameter :: tauAMPA_spinstell_to_nontuftRS=2.d0  
      real*8, parameter :: tauNMDA_spinstell_to_nontuftRS=130.d0

      real*8, parameter :: tauAMPA_tuftIB_to_suppyrRS    =2.d0 
      real*8, parameter :: tauNMDA_tuftIB_to_suppyrRS    =130.d0
      real*8, parameter :: tauAMPA_tuftIB_to_suppyrFRB   =2.d0 
      real*8, parameter :: tauNMDA_tuftIB_to_suppyrFRB   =130.d0
      real*8, parameter :: tauAMPA_tuftIB_to_supbask     =.8d0  
      real*8, parameter :: tauNMDA_tuftIB_to_supbask     =100.d0 
      real*8, parameter :: tauAMPA_tuftIB_to_supaxax     =.8d0  
      real*8, parameter :: tauNMDA_tuftIB_to_supaxax     =100.d0 
      real*8, parameter :: tauAMPA_tuftIB_to_supLTS      =1.d0  
      real*8, parameter :: tauNMDA_tuftIB_to_supLTS      =100.d0 
      real*8, parameter :: tauAMPA_tuftIB_to_spinstell   =2.d0   
      real*8, parameter :: tauNMDA_tuftIB_to_spinstell   =130.d0 
      real*8, parameter :: tauAMPA_tuftIB_to_tuftIB      =2.d0  
      real*8, parameter :: tauNMDA_tuftIB_to_tuftIB      =130.d0 
      real*8, parameter :: tauAMPA_tuftIB_to_tuftRS      =2.0d0 
      real*8, parameter :: tauNMDA_tuftIB_to_tuftRS      =130.d0 
      real*8, parameter :: tauAMPA_tuftIB_to_deepbask    =.8d0  
      real*8, parameter :: tauNMDA_tuftIB_to_deepbask    =100.d0 
      real*8, parameter :: tauAMPA_tuftIB_to_deepng      =.8d0  
      real*8, parameter :: tauNMDA_tuftIB_to_deepng      =100.d0 
      real*8, parameter :: tauAMPA_tuftIB_to_deepaxax    =.8d0  
      real*8, parameter :: tauNMDA_tuftIB_to_deepaxax    =100.d0 
      real*8, parameter :: tauAMPA_tuftIB_to_deepLTS     =1.d0  
      real*8, parameter :: tauNMDA_tuftIB_to_deepLTS     =100.d0 
      real*8, parameter :: tauAMPA_tuftIB_to_nontuftRS   =2.0d0 
      real*8, parameter :: tauNMDA_tuftIB_to_nontuftRS   =130.d0 

      real*8, parameter :: tauAMPA_tuftRS_to_suppyrRS    =2.d0 
      real*8, parameter :: tauNMDA_tuftRS_to_suppyrRS    =130.d0
      real*8, parameter :: tauAMPA_tuftRS_to_suppyrFRB   =2.d0 
      real*8, parameter :: tauNMDA_tuftRS_to_suppyrFRB   =130.d0
      real*8, parameter :: tauAMPA_tuftRS_to_supbask     =.8d0  
      real*8, parameter :: tauNMDA_tuftRS_to_supbask     =100.d0 
      real*8, parameter :: tauAMPA_tuftRS_to_supaxax     =.8d0  
      real*8, parameter :: tauNMDA_tuftRS_to_supaxax     =100.d0 
      real*8, parameter :: tauAMPA_tuftRS_to_supLTS      =1.d0  
      real*8, parameter :: tauNMDA_tuftRS_to_supLTS      =100.d0 
      real*8, parameter :: tauAMPA_tuftRS_to_spinstell   =2.d0  
      real*8, parameter :: tauNMDA_tuftRS_to_spinstell   =130.d0 
      real*8, parameter :: tauAMPA_tuftRS_to_tuftIB      =2.d0  
      real*8, parameter :: tauNMDA_tuftRS_to_tuftIB      =130.d0 
      real*8, parameter :: tauAMPA_tuftRS_to_tuftRS      =2.d0  
      real*8, parameter :: tauNMDA_tuftRS_to_tuftRS      =130.d0 
      real*8, parameter :: tauAMPA_tuftRS_to_deepbask    =.8d0  
      real*8, parameter :: tauNMDA_tuftRS_to_deepbask    =100.d0 
      real*8, parameter :: tauAMPA_tuftRS_to_deepng      =.8d0  
      real*8, parameter :: tauNMDA_tuftRS_to_deepng      =100.d0 
      real*8, parameter :: tauAMPA_tuftRS_to_deepaxax    =.8d0  
      real*8, parameter :: tauNMDA_tuftRS_to_deepaxax    =100.d0 
      real*8, parameter :: tauAMPA_tuftRS_to_deepLTS     =1.d0   
      real*8, parameter :: tauNMDA_tuftRS_to_deepLTS     =100.d0 
      real*8, parameter :: tauAMPA_tuftRS_to_nontuftRS   =2.d0  
      real*8, parameter :: tauNMDA_tuftRS_to_nontuftRS   =130.d0 

      real*8, parameter :: tauGABA_deepbask_to_spinstell =6.d0 
      real*8, parameter :: tauGABA_deepbask_to_tuftIB    =6.d0  
      real*8, parameter :: tauGABA_deepbask_to_tuftRS    =6.d0  
      real*8, parameter :: tauGABA_deepbask_to_deepbask  =3.d0  
      real*8, parameter :: tauGABA_deepbask_to_deepaxax  =3.d0  
      real*8, parameter :: tauGABA_deepbask_to_deepLTS   =3.d0  
      real*8, parameter :: tauGABA_deepbask_to_deepng    =3.d0  
      real*8, parameter :: tauGABA_deepbask_to_nontuftRS =6.d0  

      real*8, parameter :: tauGABA_deepaxax_to_suppyrRS   =6.d0 
      real*8, parameter :: tauGABA_deepaxax_to_suppyrFRB  =6.d0 
      real*8, parameter :: tauGABA_deepaxax_to_spinstell  =6.d0 
      real*8, parameter :: tauGABA_deepaxax_to_tuftIB     =6.d0 
      real*8, parameter :: tauGABA_deepaxax_to_tuftRS     =6.d0 
      real*8, parameter :: tauGABA_deepaxax_to_nontuftRS  =6.d0 

      real*8, parameter :: tauGABA_deepLTS_to_suppyrRS    =20.d0 
      real*8, parameter :: tauGABA_deepLTS_to_suppyrFRB   =20.d0 
      real*8, parameter :: tauGABA_deepLTS_to_supbask     =20.d0 
      real*8, parameter :: tauGABA_deepLTS_to_supaxax     =20.d0 
      real*8, parameter :: tauGABA_deepLTS_to_supLTS      =20.d0 
      real*8, parameter :: tauGABA_deepLTS_to_spinstell   =20.d0 
      real*8, parameter :: tauGABA_deepLTS_to_tuftIB      =20.d0 
      real*8, parameter :: tauGABA_deepLTS_to_tuftRS      =20.d0 
      real*8, parameter :: tauGABA_deepLTS_to_deepbask    =20.d0 
      real*8, parameter :: tauGABA_deepLTS_to_deepaxax    =20.d0 
      real*8, parameter :: tauGABA_deepLTS_to_deepLTS     =20.d0 
      real*8, parameter :: tauGABA_deepLTS_to_nontuftRS   =20.d0 

      real*8, parameter :: tauGABA_deepng_to_tuftIB       =6.d0
      real*8, parameter :: tauGABA_deepng_to_tuftRS       =6.d0
      real*8, parameter :: tauGABA_deepng_to_nontuftRS    =6.d0
      real*8, parameter :: tauGABA_deepng_to_spinstell    =6.d0
      real*8, parameter :: tauGABA_deepng_to_deepng       =3.d0
      real*8, parameter :: tauGABA_deepng_to_deepbask     =3.d0

      real*8, parameter :: tauAMPA_TCR_to_suppyrRS        =2.d0  
      real*8, parameter :: tauNMDA_TCR_to_suppyrRS        =130.d0
      real*8, parameter :: tauAMPA_TCR_to_suppyrFRB       =2.d0  
      real*8, parameter :: tauNMDA_TCR_to_suppyrFRB       =130.d0
c     real*8, parameter :: tauAMPA_TCR_to_supbask         =1.d0  
      real*8, parameter :: tauAMPA_TCR_to_supbask         =0.75d0  
      real*8, parameter :: tauNMDA_TCR_to_supbask         =100.d0
      real*8, parameter :: tauAMPA_TCR_to_supng           =0.75d0  
      real*8, parameter :: tauNMDA_TCR_to_supng           =100.d0
      real*8, parameter :: tauAMPA_TCR_to_supaxax         =1.d0  
      real*8, parameter :: tauNMDA_TCR_to_supaxax         =100.d0 
      real*8, parameter :: tauAMPA_TCR_to_spinstell       =2.0d0 
      real*8, parameter :: tauNMDA_TCR_to_spinstell       =130.d0
      real*8, parameter :: tauAMPA_TCR_to_tuftIB          =2.d0  
      real*8, parameter :: tauNMDA_TCR_to_tuftIB          =130.d0
      real*8, parameter :: tauAMPA_TCR_to_tuftRS          =2.d0  
      real*8, parameter :: tauNMDA_TCR_to_tuftRS          =130.d0
!     real*8, parameter :: tauAMPA_TCR_to_deepbask        =1.d0  
      real*8, parameter :: tauAMPA_TCR_to_deepbask        =0.75d0  
      real*8, parameter :: tauNMDA_TCR_to_deepbask        =100.d0
      real*8, parameter :: tauAMPA_TCR_to_deepng          =0.75d0  
      real*8, parameter :: tauNMDA_TCR_to_deepng          =100.d0
!     real*8, parameter :: tauAMPA_TCR_to_deepaxax        =1.d0  
      real*8, parameter :: tauAMPA_TCR_to_deepaxax        =0.75d0  
      real*8, parameter :: tauNMDA_TCR_to_deepaxax        =100.d0
      real*8, parameter :: tauAMPA_TCR_to_nRT             =2.0d0      
      real*8, parameter :: tauNMDA_TCR_to_nRT             =150.d0
      real*8, parameter :: tauAMPA_TCR_to_nontuftRS       =2.0d0     
      real*8, parameter :: tauNMDA_TCR_to_nontuftRS       =130.d0

!     real*8, parameter :: tauGABA1_nRT_to_TCR             =10.d0 
!     real*8, parameter :: tauGABA2_nRT_to_TCR             =30.d0 
!     real*8, parameter :: tauGABA1_nRT_to_nRT             =18.d0 
!     real*8, parameter :: tauGABA2_nRT_to_nRT             =89.d0 
! See notebook entry of 17 Feb. 2004.
! Speed these up per Huntsman & Huguenard (2000)
      real*8, parameter :: tauGABA1_nRT_to_TCR             =3.30d0 
      real*8, parameter :: tauGABA2_nRT_to_TCR             =10.d0 
      real*8, parameter :: tauGABA1_nRT_to_nRT             = 9.d0 
      real*8, parameter :: tauGABA2_nRT_to_nRT             =44.5d0 

      real*8, parameter :: tauAMPA_nontuftRS_to_suppyrRS  =2.d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_suppyrRS  =130.d0
      real*8, parameter :: tauAMPA_nontuftRS_to_suppyrFRB =2.d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_suppyrFRB =130.d0
      real*8, parameter :: tauAMPA_nontuftRS_to_supbask   =.8d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_supbask   =100.d0
      real*8, parameter :: tauAMPA_nontuftRS_to_supaxax   =.8d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_supaxax   =100.d0 
      real*8, parameter :: tauAMPA_nontuftRS_to_supLTS    =1.0d0 
      real*8, parameter :: tauNMDA_nontuftRS_to_supLTS    =100.d0
      real*8, parameter :: tauAMPA_nontuftRS_to_spinstell =2.d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_spinstell =130.d0
      real*8, parameter :: tauAMPA_nontuftRS_to_tuftIB    =2.d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_tuftIB    =130.d0
      real*8, parameter :: tauAMPA_nontuftRS_to_tuftRS    =2.d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_tuftRS    =130.d0
      real*8, parameter :: tauAMPA_nontuftRS_to_deepbask  =.8d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_deepbask  =100.d0
      real*8, parameter :: tauAMPA_nontuftRS_to_deepng    =.8d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_deepng    =100.d0
      real*8, parameter :: tauAMPA_nontuftRS_to_deepaxax  =.8d0   
      real*8, parameter :: tauNMDA_nontuftRS_to_deepaxax  =100.d0
      real*8, parameter :: tauAMPA_nontuftRS_to_deepLTS   =1.d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_deepLTS   =100.d0
      real*8, parameter :: tauAMPA_nontuftRS_to_TCR       =2.d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_TCR       =130.d0 
      real*8, parameter :: tauAMPA_nontuftRS_to_nRT       =2.0d0 
      real*8, parameter :: tauNMDA_nontuftRS_to_nRT       =100.d0 
      real*8, parameter :: tauAMPA_nontuftRS_to_nontuftRS =2.d0  
      real*8, parameter :: tauNMDA_nontuftRS_to_nontuftRS =130.d0 
c End definition of synaptic time constants

c Synaptic conductance scaling factors.
      real*8 gAMPA_suppyrRS_to_suppyrRS/0.05d-3/
c     real*8 gNMDA_suppyrRS_to_suppyrRS/0.010d-3/
c     real*8 gNMDA_suppyrRS_to_suppyrRS/1.500d-3/
      real*8 gNMDA_suppyrRS_to_suppyrRS/1.100d-3/
      real*8 gAMPA_suppyrRS_to_suppyrFRB/0.25d-3/
      real*8 gNMDA_suppyrRS_to_suppyrFRB/0.025d-3/
      real*8 gAMPA_suppyrRS_to_supbask  /1.00d-3/
      real*8 gNMDA_suppyrRS_to_supbask  /0.15d-3/
c     real*8 gAMPA_suppyrRS_to_supng    /2.00d-3/
      real*8 gAMPA_suppyrRS_to_supng    /0.20d-3/
      real*8 gNMDA_suppyrRS_to_supng    /0.01d-3/
      real*8 gAMPA_suppyrRS_to_supaxax  /1.0d-3/
      real*8 gNMDA_suppyrRS_to_supaxax  /0.15d-3/
      real*8 gAMPA_suppyrRS_to_supLTS   /2.0d-3/
      real*8 gNMDA_suppyrRS_to_supLTS   /0.15d-3/
      real*8 gAMPA_suppyrRS_to_spinstell/0.05d-3/
      real*8 gNMDA_suppyrRS_to_spinstell/0.05d-3/
      real*8 gAMPA_suppyrRS_to_tuftIB   /0.10d-3/
      real*8 gNMDA_suppyrRS_to_tuftIB   /0.01d-3/
      real*8 gAMPA_suppyrRS_to_tuftRS   /0.10d-3/
      real*8 gNMDA_suppyrRS_to_tuftRS   /0.01d-3/
      real*8 gAMPA_suppyrRS_to_deepbask /1.00d-3/
      real*8 gNMDA_suppyrRS_to_deepbask /0.10d-3/
      real*8 gAMPA_suppyrRS_to_deepaxax /1.00d-3/
      real*8 gNMDA_suppyrRS_to_deepaxax /0.10d-3/
      real*8 gAMPA_suppyrRS_to_deepLTS  /1.00d-3/
      real*8 gNMDA_suppyrRS_to_deepLTS  /0.05d-3/
      real*8 gAMPA_suppyrRS_to_nontuftRS/0.10d-3/
      real*8 gNMDA_suppyrRS_to_nontuftRS/0.05d-3/

      real*8 gAMPA_suppyrFRB_to_suppyrRS/0.25d-3/
      real*8 gNMDA_suppyrFRB_to_suppyrRS/0.025d-3/
      real*8 gAMPA_suppyrFRB_to_suppyrFRB/0.25d-3/
      real*8 gNMDA_suppyrFRB_to_suppyrFRB/.025d-3/
      real*8 gAMPA_suppyrFRB_to_supbask  /3.00d-3/
      real*8 gNMDA_suppyrFRB_to_supbask  /0.10d-3/
c     real*8 gAMPA_suppyrFRB_to_supng    /2.00d-3/
      real*8 gAMPA_suppyrFRB_to_supng    /0.20d-3/
      real*8 gNMDA_suppyrFRB_to_supng    /0.01d-3/
      real*8 gAMPA_suppyrFRB_to_supaxax  /3.0d-3/
      real*8 gNMDA_suppyrFRB_to_supaxax  /0.10d-3/
      real*8 gAMPA_suppyrFRB_to_supLTS   /2.0d-3/
      real*8 gNMDA_suppyrFRB_to_supLTS   /0.10d-3/
      real*8 gAMPA_suppyrFRB_to_spinstell/0.10d-3/
      real*8 gNMDA_suppyrFRB_to_spinstell/0.01d-3/
      real*8 gAMPA_suppyrFRB_to_tuftIB   /0.10d-3/
      real*8 gNMDA_suppyrFRB_to_tuftIB   /0.01d-3/
      real*8 gAMPA_suppyrFRB_to_tuftRS   /0.10d-3/
      real*8 gNMDA_suppyrFRB_to_tuftRS   /0.01d-3/
      real*8 gAMPA_suppyrFRB_to_deepbask /1.00d-3/
      real*8 gNMDA_suppyrFRB_to_deepbask /0.10d-3/
      real*8 gAMPA_suppyrFRB_to_deepaxax /1.00d-3/
      real*8 gNMDA_suppyrFRB_to_deepaxax /0.10d-3/
      real*8 gAMPA_suppyrFRB_to_deepLTS  /1.00d-3/
      real*8 gNMDA_suppyrFRB_to_deepLTS  /0.05d-3/
      real*8 gAMPA_suppyrFRB_to_nontuftRS/0.05d-3/
      real*8 gNMDA_suppyrFRB_to_nontuftRS/0.05d-3/

      real*8 gGABA_supbask_to_suppyrRS   /0.0d-3/
!     real*8 gGABA_supbask_to_suppyrRS   /3.0d-3/
c     real*8 gGABA_supbask_to_suppyrFRB  /1.2d-3/
      real*8 gGABA_supbask_to_suppyrFRB  /0.0d-3/ ! try to make gamma robust
      real*8 gGABA_supbask_to_supbask    /0.2d-3/
      real*8 gGABA_supbask_to_supng      /0.2d-3/
      real*8 gGABA_supbask_to_supaxax    /0.2d-3/
      real*8 gGABA_supbask_to_supLTS     /0.5d-3/
      real*8 gGABA_supbask_to_spinstell  /0.2d-3/ 

      real*8 gGABA_supng_to_suppyrRS     /0.0d-3/
!     real*8 gGABA_supng_to_suppyrRS     /0.1d-3/
      real*8 gGABA_supng_to_suppyrFRB    /0.0d-3/
      real*8 gGABA_supng_to_nontuftRS    /0.8d-3/
      real*8 gGABA_supng_to_tuftIB       /0.8d-3/
      real*8 gGABA_supng_to_tuftRS       /0.8d-3/
      real*8 gGABA_supng_to_supng        /0.2d-3/
      real*8 gGABA_supng_to_supbask      /0.2d-3/

c     real*8 gGABAB_supng_to_suppyrRS    /0.001d-3/
      real*8 gGABAB_supng_to_suppyrRS    /0.000d-3/
c     real*8 gGABAB_supng_to_suppyrFRB   /0.025d-3/
      real*8 gGABAB_supng_to_suppyrFRB   /0.005d-3/
      real*8 gGABAB_supng_to_nontuftRS   /0.0003d-3/
      real*8 gGABAB_supng_to_tuftIB      /0.035d-3/
      real*8 gGABAB_supng_to_tuftRS      /0.025d-3/

      real*8 gGABA_supaxax_to_suppyrRS   /0.25d-3/
c     real*8 gGABA_supaxax_to_suppyrRS   /0.00d-3/
      real*8 gGABA_supaxax_to_suppyrFRB  /0.0d-3/
      real*8 gGABA_supaxax_to_spinstell  /0.0d-3/ ! if main inhib. to spinstell from deep int.
      real*8 gGABA_supaxax_to_tuftIB     /1.0d-3/
      real*8 gGABA_supaxax_to_tuftRS     /1.0d-3/
      real*8 gGABA_supaxax_to_nontuftRS  /0.5d-3/


      real*8 gGABA_supLTS_to_suppyrRS    /.00d-3/
c     real*8 gGABA_supLTS_to_suppyrRS    /.01d-3/
!     real*8 gGABA_supLTS_to_suppyrRS    /1.500d-3/
      real*8 gGABA_supLTS_to_suppyrFRB   /.00d-3/
      real*8 gGABA_supLTS_to_supbask     /.01d-3/
      real*8 gGABA_supLTS_to_supaxax     /.01d-3/
      real*8 gGABA_supLTS_to_supLTS      /.05d-3/
c     real*8 gGABA_supLTS_to_spinstell   /.01d-3/
      real*8 gGABA_supLTS_to_spinstell   /1.50d-3/
      real*8 gGABA_supLTS_to_tuftIB      /.02d-3/
      real*8 gGABA_supLTS_to_tuftRS      /.02d-3/
      real*8 gGABA_supLTS_to_deepbask    /.01d-3/
      real*8 gGABA_supLTS_to_deepaxax    /.00d-3/
      real*8 gGABA_supLTS_to_deepLTS     /.05d-3/
      real*8 gGABA_supLTS_to_nontuftRS   /.02d-3/

      real*8 gAMPA_spinstell_to_suppyrRS /0.1d-3/
      real*8 gNMDA_spinstell_to_suppyrRS /0.1d-3/
      real*8 gAMPA_spinstell_to_suppyrFRB/1.0d-3/
      real*8 gNMDA_spinstell_to_suppyrFRB/0.1d-3/
      real*8 gAMPA_spinstell_to_supbask  /1.0d-3/
      real*8 gNMDA_spinstell_to_supbask  /.15d-3/
      real*8 gAMPA_spinstell_to_supaxax  /1.0d-3/
      real*8 gNMDA_spinstell_to_supaxax  /.15d-3/
c     real*8 gAMPA_spinstell_to_supLTS   /1.0d-3/
      real*8 gAMPA_spinstell_to_supLTS   /2.0d-3/
      real*8 gNMDA_spinstell_to_supLTS   /.15d-3/
      real*8 gAMPA_spinstell_to_spinstell/0.2d-3/
c     real*8 gAMPA_spinstell_to_spinstell/0.6d-3/
      real*8 gNMDA_spinstell_to_spinstell/0.1d-3/
      real*8 gAMPA_spinstell_to_tuftIB   /1.0d-3/
      real*8 gNMDA_spinstell_to_tuftIB   /0.1d-3/
      real*8 gAMPA_spinstell_to_tuftRS   /1.0d-3/
      real*8 gNMDA_spinstell_to_tuftRS   /0.1d-3/
      real*8 gAMPA_spinstell_to_deepbask /1.0d-3/
      real*8 gNMDA_spinstell_to_deepbask /.15d-3/
c     real*8 gAMPA_spinstell_to_deepng   /1.0d-3/
      real*8 gAMPA_spinstell_to_deepng   /0.4d-3/
      real*8 gNMDA_spinstell_to_deepng   /.15d-3/
      real*8 gAMPA_spinstell_to_deepaxax /1.0d-3/
      real*8 gNMDA_spinstell_to_deepaxax /.15d-3/
      real*8 gAMPA_spinstell_to_deepLTS  /0.5d-3/
      real*8 gNMDA_spinstell_to_deepLTS  /.01d-3/
      real*8 gAMPA_spinstell_to_nontuftRS/0.1d-3/
      real*8 gNMDA_spinstell_to_nontuftRS/0.1d-3/

      real*8 gAMPA_tuftIB_to_suppyrRS    /0.0d-3/
      real*8 gNMDA_tuftIB_to_suppyrRS    /0.00d-3/
      real*8 gAMPA_tuftIB_to_suppyrFRB   /0.0d-3/
      real*8 gNMDA_tuftIB_to_suppyrFRB   /0.00d-3/
      real*8 gAMPA_tuftIB_to_supbask     /0.0d-3/
      real*8 gNMDA_tuftIB_to_supbask     /0.00d-3/
      real*8 gAMPA_tuftIB_to_supaxax     /0.0d-3/
      real*8 gNMDA_tuftIB_to_supaxax     /0.00d-3/
      real*8 gAMPA_tuftIB_to_supLTS      /0.0d-3/
      real*8 gNMDA_tuftIB_to_supLTS      /0.00d-3/
      real*8 gAMPA_tuftIB_to_spinstell   /0.0d-3/
      real*8 gNMDA_tuftIB_to_spinstell   /0.00d-3/
c     real*8 gAMPA_tuftIB_to_tuftIB      /2.0d-3/
      real*8 gAMPA_tuftIB_to_tuftIB      /3.0d-3/
c     real*8 gNMDA_tuftIB_to_tuftIB      /0.20d-3/
      real*8 gNMDA_tuftIB_to_tuftIB      /0.07d-3/
      real*8 gAMPA_tuftIB_to_tuftRS      /2.0d-3/
      real*8 gNMDA_tuftIB_to_tuftRS      /0.20d-3/
      real*8 gAMPA_tuftIB_to_deepbask    /3.0d-3/
      real*8 gNMDA_tuftIB_to_deepbask    /0.15d-3/
c     real*8 gAMPA_tuftIB_to_deepng      /2.0d-3/
      real*8 gAMPA_tuftIB_to_deepng      /1.2d-3/
c     real*8 gAMPA_tuftIB_to_deepng      /3.0d-3/
      real*8 gNMDA_tuftIB_to_deepng      /0.10d-3/
c     real*8 gNMDA_tuftIB_to_deepng      /0.50d-3/
      real*8 gAMPA_tuftIB_to_deepaxax    /3.0d-3/
      real*8 gNMDA_tuftIB_to_deepaxax    /0.15d-3/
      real*8 gAMPA_tuftIB_to_deepLTS     /1.0d-3/
c     real*8 gNMDA_tuftIB_to_deepLTS     /0.05d-3/
      real*8 gNMDA_tuftIB_to_deepLTS     /0.01d-3/
      real*8 gAMPA_tuftIB_to_nontuftRS   /0.40d-3/
      real*8 gNMDA_tuftIB_to_nontuftRS   /0.01d-3/

      real*8 gAMPA_tuftRS_to_suppyrRS    /0.0d-3/
      real*8 gNMDA_tuftRS_to_suppyrRS    /0.00d-3/
      real*8 gAMPA_tuftRS_to_suppyrFRB   /0.0d-3/
      real*8 gNMDA_tuftRS_to_suppyrFRB   /0.00d-3/
      real*8 gAMPA_tuftRS_to_supbask     /0.0d-3/
      real*8 gNMDA_tuftRS_to_supbask     /0.00d-3/
      real*8 gAMPA_tuftRS_to_supaxax     /0.0d-3/
      real*8 gNMDA_tuftRS_to_supaxax     /0.00d-3/
      real*8 gAMPA_tuftRS_to_supLTS      /0.0d-3/
      real*8 gNMDA_tuftRS_to_supLTS      /0.00d-3/
      real*8 gAMPA_tuftRS_to_spinstell   /0.0d-3/
      real*8 gNMDA_tuftRS_to_spinstell   /0.00d-3/
      real*8 gAMPA_tuftRS_to_tuftIB      /1.0d-3/
      real*8 gNMDA_tuftRS_to_tuftIB      /0.10d-3/
      real*8 gAMPA_tuftRS_to_tuftRS      /1.0d-3/
      real*8 gNMDA_tuftRS_to_tuftRS      /0.10d-3/
      real*8 gAMPA_tuftRS_to_deepbask    /3.0d-3/
      real*8 gNMDA_tuftRS_to_deepbask    /0.10d-3/
c     real*8 gAMPA_tuftRS_to_deepng      /2.0d-3/
      real*8 gAMPA_tuftRS_to_deepng      /0.8d-3/
      real*8 gNMDA_tuftRS_to_deepng      /0.10d-3/
      real*8 gAMPA_tuftRS_to_deepaxax    /3.0d-3/
      real*8 gNMDA_tuftRS_to_deepaxax    /0.10d-3/
      real*8 gAMPA_tuftRS_to_deepLTS     /1.0d-3/
      real*8 gNMDA_tuftRS_to_deepLTS     /0.01d-3/
      real*8 gAMPA_tuftRS_to_nontuftRS   /0.1d-3/
      real*8 gNMDA_tuftRS_to_nontuftRS   /0.01d-3/

!     real*8 gGABA_deepbask_to_spinstell /1.0d-3/
      real*8 gGABA_deepbask_to_spinstell /0.00d-3/ ! ? suppress spiny stellate bursts ?
c     real*8 gGABA_deepbask_to_tuftIB    /0.7d-3/
      real*8 gGABA_deepbask_to_tuftIB    /0.3d-3/
      real*8 gGABA_deepbask_to_tuftRS    /0.7d-3/
      real*8 gGABA_deepbask_to_deepbask  /0.2d-3/
c     real*8 gGABA_deepbask_to_deepng    /0.2d-3/
      real*8 gGABA_deepbask_to_deepng    /0.1d-3/
      real*8 gGABA_deepbask_to_deepaxax  /0.2d-3/
      real*8 gGABA_deepbask_to_deepLTS   /0.7d-3/
      real*8 gGABA_deepbask_to_nontuftRS /1.2d-3/

c     real*8 gGABA_deepng_to_tuftIB      /0.8d-3/
      real*8 gGABA_deepng_to_tuftIB      /0.3d-3/
      real*8 gGABA_deepng_to_tuftRS      /0.8d-3/
      real*8 gGABA_deepng_to_nontuftRS   /0.8d-3/
      real*8 gGABA_deepng_to_spinstell   /0.0d-3/
c     real*8 gGABA_deepng_to_deepng      /0.2d-3/
      real*8 gGABA_deepng_to_deepng      /0.1d-3/
      real*8 gGABA_deepng_to_deepbask    /0.2d-3/

c     real*8 gGABAB_deepng_to_tuftIB      /0.015d-3/
      real*8 gGABAB_deepng_to_tuftIB      /0.030d-3/
      real*8 gGABAB_deepng_to_tuftRS      /0.025d-3/
      real*8 gGABAB_deepng_to_nontuftRS   /0.0003d-3/
      real*8 gGABAB_deepng_to_spinstell   /0.000d-3/

      real*8 gGABA_deepaxax_to_suppyrRS   /0.0d-3/
      real*8 gGABA_deepaxax_to_suppyrFRB  /0.0d-3/
!     real*8 gGABA_deepaxax_to_spinstell  /1.0d-3/
      real*8 gGABA_deepaxax_to_spinstell  /0.0d-3/ ! ? suppress spiny stellate bursts ?
      real*8 gGABA_deepaxax_to_tuftIB     /1.0d-3/
      real*8 gGABA_deepaxax_to_tuftRS     /1.0d-3/
      real*8 gGABA_deepaxax_to_nontuftRS  /1.00d-3/

      real*8 gGABA_deepLTS_to_suppyrRS    /.00d-3/
      real*8 gGABA_deepLTS_to_suppyrFRB   /.00d-3/
      real*8 gGABA_deepLTS_to_supbask     /.00d-3/
      real*8 gGABA_deepLTS_to_supaxax     /.00d-3/
      real*8 gGABA_deepLTS_to_supLTS      /.00d-3/
c     real*8 gGABA_deepLTS_to_spinstell   /.00d-3/
      real*8 gGABA_deepLTS_to_spinstell   /1.0d-3/
!     real*8 gGABA_deepLTS_to_tuftIB      /.02d-3/
      real*8 gGABA_deepLTS_to_tuftIB      /.05d-3/ ! will this help suppress bursting?
      real*8 gGABA_deepLTS_to_tuftRS      /.02d-3/
      real*8 gGABA_deepLTS_to_deepbask    /.01d-3/
      real*8 gGABA_deepLTS_to_deepaxax    /.01d-3/
c     real*8 gGABA_deepLTS_to_deepLTS     /.05d-3/
      real*8 gGABA_deepLTS_to_deepLTS     /.20d-3/
      real*8 gGABA_deepLTS_to_nontuftRS   /.20d-3/

      real*8 gAMPA_TCR_to_suppyrRS        /0.0d-3/
      real*8 gNMDA_TCR_to_suppyrRS        /0.00d-3/
      real*8 gAMPA_TCR_to_suppyrFRB       /0.0d-3/
      real*8 gNMDA_TCR_to_suppyrFRB       /0.00d-3/
!     real*8 gAMPA_TCR_to_supbask         /1.0d-3/
      real*8 gAMPA_TCR_to_supbask         /0.0d-3/
! try a variation in which main feedforward inhibtion from thalamus
! is via deep interneurons.  May be necessary later to include special
! layer 4 interneurons
!     real*8 gNMDA_TCR_to_supbask         /.10d-3/
      real*8 gNMDA_TCR_to_supbask         /.00d-3/
      real*8 gAMPA_TCR_to_supng           /0.0d-3/
      real*8 gNMDA_TCR_to_supng           /0.0d-3/
!     real*8 gAMPA_TCR_to_supaxax         /1.0d-3/
      real*8 gAMPA_TCR_to_supaxax         /0.0d-3/
!     real*8 gNMDA_TCR_to_supaxax         /.10d-3/
      real*8 gNMDA_TCR_to_supaxax         /.00d-3/
      real*8 gAMPA_TCR_to_spinstell       /0.0d-3/
      real*8 gNMDA_TCR_to_spinstell       /.00d-3/
      real*8 gAMPA_TCR_to_tuftIB          /1.5d-3/
      real*8 gNMDA_TCR_to_tuftIB          /.15d-3/
      real*8 gAMPA_TCR_to_tuftRS          /1.5d-3/
      real*8 gNMDA_TCR_to_tuftRS          /.15d-3/
!     real*8 gAMPA_TCR_to_deepbask        /1.0d-3/
      real*8 gAMPA_TCR_to_deepbask        /1.5d-3/
!     real*8 gAMPA_TCR_to_deepbask        /0.0d-3/ ! try for very fast FF excit.
      real*8 gNMDA_TCR_to_deepbask        /.10d-3/
      real*8 gAMPA_TCR_to_deepng          /1.5d-3/
      real*8 gNMDA_TCR_to_deepng          /0.1d-3/
      real*8 gAMPA_TCR_to_deepaxax        /1.0d-3/
!     real*8 gAMPA_TCR_to_deepaxax        /0.0d-3/ ! try for very fast FF excit.
      real*8 gNMDA_TCR_to_deepaxax        /.10d-3/
      real*8 gAMPA_TCR_to_nRT             /0.75d-3/   
      real*8 gNMDA_TCR_to_nRT             /.15d-3/
      real*8 gAMPA_TCR_to_nontuftRS       /0.0d-3/    
      real*8 gNMDA_TCR_to_nontuftRS       /.00d-3/

      real*8 gGABAB_nRT_to_TCR            /0.02d-3/
!     real*8 gGABA_nRT_to_TCR             /1.0d-3/
      real*8 gGABA_nRT_to_TCR(num_nRT)
! Values here need to be set below  
      real*8 gGABA_nRT_to_nRT             /0.30d-3/
      real*8 gGABAB_nRT_to_nRT            /0.020d-3/

      real*8 gAMPA_nontuftRS_to_suppyrRS  /0.0d-3/
      real*8 gNMDA_nontuftRS_to_suppyrRS  /0.00d-3/
      real*8 gAMPA_nontuftRS_to_suppyrFRB /0.0d-3/
      real*8 gNMDA_nontuftRS_to_suppyrFRB /0.00d-3/
      real*8 gAMPA_nontuftRS_to_supbask   /0.0d-3/
      real*8 gNMDA_nontuftRS_to_supbask   /0.0d-3/
      real*8 gAMPA_nontuftRS_to_supaxax   /0.0d-3/
      real*8 gNMDA_nontuftRS_to_supaxax   /0.0d-3/
      real*8 gAMPA_nontuftRS_to_supLTS    /0.0d-3/
      real*8 gNMDA_nontuftRS_to_supLTS    /0.0d-3/
      real*8 gAMPA_nontuftRS_to_spinstell /0.0d-3/
      real*8 gNMDA_nontuftRS_to_spinstell /0.00d-3/
      real*8 gAMPA_nontuftRS_to_tuftIB    /1.0d-3/
      real*8 gNMDA_nontuftRS_to_tuftIB    /0.1d-3/
      real*8 gAMPA_nontuftRS_to_tuftRS    /1.0d-3/
      real*8 gNMDA_nontuftRS_to_tuftRS    /0.1d-3/
      real*8 gAMPA_nontuftRS_to_deepbask  /3.0d-3/
      real*8 gNMDA_nontuftRS_to_deepbask  /.10d-3/
c     real*8 gAMPA_nontuftRS_to_deepng    /2.0d-3/
      real*8 gAMPA_nontuftRS_to_deepng    /0.8d-3/
      real*8 gNMDA_nontuftRS_to_deepng    /.10d-3/
      real*8 gAMPA_nontuftRS_to_deepaxax  /3.0d-3/
      real*8 gNMDA_nontuftRS_to_deepaxax  /.00d-3/
      real*8 gAMPA_nontuftRS_to_deepLTS   /1.0d-3/
      real*8 gNMDA_nontuftRS_to_deepLTS   /.01d-3/
      real*8 gAMPA_nontuftRS_to_TCR       /.15d-3/ ! make this small
      real*8 gNMDA_nontuftRS_to_TCR       /.015d-3/
      real*8 gAMPA_nontuftRS_to_nRT       /0.5d-3/
      real*8 gNMDA_nontuftRS_to_nRT       /0.05d-3/
      real*8 gAMPA_nontuftRS_to_nontuftRS /1.5d-3/
      real*8 gNMDA_nontuftRS_to_nontuftRS /0.01d-3/
c End defining synaptic conductance scaling factors

c Begin definition of compartments where synaptic connections
c can form.
       INTEGER compallow_suppyrRS_to_suppyrRS 
     &  (ncompallow_suppyrRS_to_suppyrRS)
     & /41,42,43,44,45,46,47,48,49,50,51,52,53,
     &  54,55,56,57,58,59,60,61,62,63,64,65,
     &  66,67,68/
c    &  /2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26,
c    & 27,28,29,30,31,32,33,10,11,12,13,22,23,24,25,
c    & 34,35,36,37/
       INTEGER compallow_suppyrRS_to_suppyrFRB
     &  (ncompallow_suppyrRS_to_suppyrFRB)
     &  /2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26,
     & 27,28,29,30,31,32,33,10,11,12,13,22,23,24,25,
     & 34,35,36,37/
       INTEGER compallow_suppyrRS_to_supbask  
     &  (ncompallow_suppyrRS_to_supbask  )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrRS_to_supng    
     &  (ncompallow_suppyrRS_to_supng    )
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_suppyrRS_to_supaxax  
     &  (ncompallow_suppyrRS_to_supaxax  )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrRS_to_supLTS   
     &  (ncompallow_suppyrRS_to_supLTS   )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrRS_to_spinstell
     &  (ncompallow_suppyrRS_to_spinstell)
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrRS_to_tuftIB   
     &  (ncompallow_suppyrRS_to_tuftIB   )
     &  /39,40,41,42,43,44,45,46/
       INTEGER compallow_suppyrRS_to_tuftRS   
     &  (ncompallow_suppyrRS_to_tuftRS   )
     &  /39,40,41,42,43,44,45,46/
       INTEGER compallow_suppyrRS_to_deepbask 
     &  (ncompallow_suppyrRS_to_deepbask )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrRS_to_deepaxax 
     &  (ncompallow_suppyrRS_to_deepaxax )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrRS_to_deepLTS  
     &  (ncompallow_suppyrRS_to_deepLTS  )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrRS_to_nontuftRS
     &  (ncompallow_suppyrRS_to_nontuftRS)
     &  /38,39,40,41,42,43,44/

       INTEGER compallow_suppyrFRB_to_suppyrRS 
     &  (ncompallow_suppyrFRB_to_suppyrRS)
     &  /2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26,
     & 27,28,29,30,31,32,33,10,11,12,13,22,23,24,25,
     & 34,35,36,37/
       INTEGER compallow_suppyrFRB_to_suppyrFRB
     &  (ncompallow_suppyrFRB_to_suppyrFRB)
     &  /2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,26,
     & 27,28,29,30,31,32,33,10,11,12,13,22,23,24,25,
     & 34,35,36,37/
       INTEGER compallow_suppyrFRB_to_supbask  
     &  (ncompallow_suppyrFRB_to_supbask  )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrFRB_to_supng    
     &  (ncompallow_suppyrFRB_to_supng    )
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_suppyrFRB_to_supaxax  
     &  (ncompallow_suppyrFRB_to_supaxax  )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrFRB_to_supLTS   
     &  (ncompallow_suppyrFRB_to_supLTS   )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrFRB_to_spinstell
     &  (ncompallow_suppyrFRB_to_spinstell)
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrFRB_to_tuftIB   
     &  (ncompallow_suppyrFRB_to_tuftIB   )
     &  /39,40,41,42,43,44,45,46/
       INTEGER compallow_suppyrFRB_to_tuftRS   
     &  (ncompallow_suppyrFRB_to_tuftRS   )
     &  /39,40,41,42,43,44,45,46/
       INTEGER compallow_suppyrFRB_to_deepbask 
     &  (ncompallow_suppyrFRB_to_deepbask )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrFRB_to_deepaxax 
     &  (ncompallow_suppyrFRB_to_deepaxax )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrFRB_to_deepLTS  
     &  (ncompallow_suppyrFRB_to_deepLTS  )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_suppyrFRB_to_nontuftRS
     &  (ncompallow_suppyrFRB_to_nontuftRS)
     &  /38,39,40,41,42,43,44/

       INTEGER compallow_supbask_to_suppyrRS
     &  (ncompallow_supbask_to_suppyrRS)
     & /1,2,3,4,5,6,7,8,9,38,39/
       INTEGER compallow_supbask_to_suppyrFRB
     &  (ncompallow_supbask_to_suppyrFRB)
     & /1,2,3,4,5,6,7,8,9,38,39/
       INTEGER compallow_supbask_to_supbask  
     &  (ncompallow_supbask_to_supbask  )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_supbask_to_supng    
     &  (ncompallow_supbask_to_supng    )
     & /2,15,28,41/
       INTEGER compallow_supbask_to_supaxax  
     &  (ncompallow_supbask_to_supaxax  )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_supbask_to_supLTS   
     &  (ncompallow_supbask_to_supLTS   )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_supbask_to_spinstell
     &  (ncompallow_supbask_to_spinstell)
     &  /1,2,15,28,41/

       INTEGER compallow_supng_to_suppyrRS 
     &  (ncompallow_supng_to_suppyrRS )
     & /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     & 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
     & 41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,
     & 58,59,60,61,62,63,64,65,66,67,68/
       INTEGER compallow_supng_to_suppyrFRB
     &  (ncompallow_supng_to_suppyrFRB)
     & /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     & 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
     & 41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,
     & 58,59,60,61,62,63,64,65,66,67,68/
       INTEGER compallow_supng_to_nontuftRS
     &  (ncompallow_supng_to_nontuftRS)
     & /40,41,42,43,44/
       INTEGER compallow_supng_to_tuftIB   
     &  (ncompallow_supng_to_tuftIB   )
     & /42,43,44,45,46,47,48,49,50,51,52,53,54,55/
       INTEGER compallow_supng_to_tuftRS   
     &  (ncompallow_supng_to_tuftRS   )
     & /42,43,44,45,46,47,48,49,50,51,52,53,54,55/
       INTEGER compallow_supng_to_supng    
     &  (ncompallow_supng_to_supng    )
     & /2,1,28,41/
       INTEGER compallow_supng_to_supbask  
     &  (ncompallow_supng_to_supbask  )
     & /2,1,28,41/

       INTEGER compallow_supLTS_to_suppyrRS
     &  (ncompallow_supLTS_to_suppyrRS)
     & /14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
     &  31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,
     &  50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
     &  67,68/
       INTEGER compallow_supLTS_to_suppyrFRB
     &  (ncompallow_supLTS_to_suppyrFRB)
     & /14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
     &  31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,
     &  50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
     &  67,68/
       INTEGER compallow_supLTS_to_supbask  
     &  (ncompallow_supLTS_to_supbask)  
     & /5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
     &  26,27,31,32,33,34,35,36,37,38,39,40,
     &  44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_supLTS_to_supaxax  
     &  (ncompallow_supLTS_to_supaxax)  
     & /5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
     &  26,27,31,32,33,34,35,36,37,38,39,40,
     &  44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_supLTS_to_supLTS   
     &  (ncompallow_supLTS_to_supLTS )  
     & /5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
     &  26,27,31,32,33,34,35,36,37,38,39,40,
     &  44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_supLTS_to_spinstell
     &  (ncompallow_supLTS_to_spinstell)
     & /5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
     &  26,27,31,32,33,34,35,36,37,38,39,40,
     &  44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_supLTS_to_tuftIB   
     &  (ncompallow_supLTS_to_tuftIB   )
     & / 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
     &   29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,
     &   48,49,50,51,52,53,54,55/
       INTEGER compallow_supLTS_to_tuftRS   
     &  (ncompallow_supLTS_to_tuftRS   )
     & / 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
     &   29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,
     &   48,49,50,51,52,53,54,55/
       INTEGER compallow_supLTS_to_deepbask 
     &  (ncompallow_supLTS_to_deepbask )
     & / 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
     &   47,48,49,50,51/ 
       INTEGER compallow_supLTS_to_deepaxax 
     &  (ncompallow_supLTS_to_deepaxax )
     & / 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
     &   47,48,49,50,51/ 
       INTEGER compallow_supLTS_to_deepLTS  
     &  (ncompallow_supLTS_to_deepLTS  )
     & / 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
     &   47,48,49,50,51/ 
       INTEGER compallow_supLTS_to_nontuftRS
     &  (ncompallow_supLTS_to_nontuftRS)
     & / 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
     &   29,30,31,32,33,34,38,39,40,41,42,43,44/ 

c      INTEGER compallow_spinstell_to_suppyrRS
c    &   (ncompallow_spinstell_to_suppyrRS)
c    & / 40,41,42,43,44,45,46,47,48,49,50,51,52/
c      INTEGER compallow_spinstell_to_suppyrFRB
c    &   (ncompallow_spinstell_to_suppyrFRB)
c    & / 40,41,42,43,44,45,46,47,48,49,50,51,52/
! 3 Mar. 2004: Feldmeyer, ..., Sakmann, J Physiol 2002 assert
! that in barrel ctx, spiny stellates go to basal dendrites of
! layer 2/3 pyramids
       INTEGER compallow_spinstell_to_suppyrRS
     &   (ncompallow_spinstell_to_suppyrRS)
     & /  2, 3, 4, 5, 6, 7, 8, 9,14,15,16,17,18,19,20,21,
     &   26,27,28,29,30,31,32,33/
       INTEGER compallow_spinstell_to_suppyrFRB
     &   (ncompallow_spinstell_to_suppyrFRB)
     & /  2, 3, 4, 5, 6, 7, 8, 9,14,15,16,17,18,19,20,21,
     &   26,27,28,29,30,31,32,33/
       INTEGER compallow_spinstell_to_supbask  
     &   (ncompallow_spinstell_to_supbask  )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_spinstell_to_supaxax  
     &   (ncompallow_spinstell_to_supaxax  )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_spinstell_to_supLTS   
     &   (ncompallow_spinstell_to_supLTS   )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_spinstell_to_spinstell
     &   (ncompallow_spinstell_to_spinstell)
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_spinstell_to_tuftIB   
     &   (ncompallow_spinstell_to_tuftIB   )
     &  / 7,8,9,10,11,12,36,37,38,39,40,41/
       INTEGER compallow_spinstell_to_tuftRS   
     &   (ncompallow_spinstell_to_tuftRS   )
     &  / 7,8,9,10,11,12,36,37,38,39,40,41/
       INTEGER compallow_spinstell_to_deepbask 
     &   (ncompallow_spinstell_to_deepbask )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_spinstell_to_deepng   
     &   (ncompallow_spinstell_to_deepng   )
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_spinstell_to_deepaxax 
     &   (ncompallow_spinstell_to_deepaxax )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_spinstell_to_deepLTS  
     &   (ncompallow_spinstell_to_deepLTS  )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_spinstell_to_nontuftRS
     &   (ncompallow_spinstell_to_nontuftRS)
     &  / 37,38,39,40,41/

       INTEGER compallow_tuftIB_to_suppyrRS
     &   (ncompallow_tuftIB_to_suppyrRS)
     & / 40,41,42,43,44,45,46,47,48,49,50,51,52/
       INTEGER compallow_tuftIB_to_suppyrFRB
     &   (ncompallow_tuftIB_to_suppyrFRB)
     & / 40,41,42,43,44,45,46,47,48,49,50,51,52/
       INTEGER compallow_tuftIB_to_supbask  
     &   (ncompallow_tuftIB_to_supbask)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftIB_to_supaxax  
     &   (ncompallow_tuftIB_to_supaxax)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftIB_to_supLTS   
     &   (ncompallow_tuftIB_to_supLTS )  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftIB_to_spinstell
     &   (ncompallow_tuftIB_to_spinstell) 
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftIB_to_tuftIB   
     &   (ncompallow_tuftIB_to_tuftIB)    
     &  / 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
     &   38,39,40,41,42,43,44,45,46,47/
       INTEGER compallow_tuftIB_to_tuftRS   
     &   (ncompallow_tuftIB_to_tuftRS)    
     &  / 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
     &   38,39,40,41,42,43,44,45,46,47/
       INTEGER compallow_tuftIB_to_deepbask 
     &   (ncompallow_tuftIB_to_deepbask)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftIB_to_deepng   
     &   (ncompallow_tuftIB_to_deepng  )  
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_tuftIB_to_deepaxax 
     &   (ncompallow_tuftIB_to_deepaxax)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftIB_to_deepLTS  
     &   (ncompallow_tuftIB_to_deepLTS )  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftIB_to_nontuftRS
     &   (ncompallow_tuftIB_to_nontuftRS) 
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &   37,38,39,40,41,42,43,44/

       INTEGER compallow_tuftRS_to_suppyrRS
     &   (ncompallow_tuftRS_to_suppyrRS)
     & / 40,41,42,43,44,45,46,47,48,49,50,51,52/
       INTEGER compallow_tuftRS_to_suppyrFRB
     &   (ncompallow_tuftRS_to_suppyrFRB)
     & / 40,41,42,43,44,45,46,47,48,49,50,51,52/
       INTEGER compallow_tuftRS_to_supbask  
     &   (ncompallow_tuftRS_to_supbask)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftRS_to_supaxax  
     &   (ncompallow_tuftRS_to_supaxax)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftRS_to_supLTS   
     &   (ncompallow_tuftRS_to_supLTS )  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftRS_to_spinstell
     &   (ncompallow_tuftRS_to_spinstell) 
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftRS_to_tuftIB   
     &   (ncompallow_tuftRS_to_tuftIB)    
     &  / 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
     &   38,39,40,41,42,43,44,45,46,47/
       INTEGER compallow_tuftRS_to_tuftRS   
     &   (ncompallow_tuftRS_to_tuftRS)    
     &  / 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
     &   38,39,40,41,42,43,44,45,46,47/
       INTEGER compallow_tuftRS_to_deepbask 
     &   (ncompallow_tuftRS_to_deepbask)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftRS_to_deepng   
     &   (ncompallow_tuftRS_to_deepng  )  
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_tuftRS_to_deepaxax 
     &   (ncompallow_tuftRS_to_deepaxax)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftRS_to_deepLTS  
     &   (ncompallow_tuftRS_to_deepLTS )  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_tuftRS_to_nontuftRS
     &   (ncompallow_tuftRS_to_nontuftRS) 
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &   37,38,39,40,41,42,43,44/

       INTEGER compallow_deepbask_to_spinstell
     &   (ncompallow_deepbask_to_spinstell)
     &  /1,2,15,28,41/
       INTEGER compallow_deepbask_to_tuftIB   
     &   (ncompallow_deepbask_to_tuftIB)   
     & / 1,2,3,4,5,6,35,36/
       INTEGER compallow_deepbask_to_tuftRS   
     &   (ncompallow_deepbask_to_tuftRS)   
     & / 1,2,3,4,5,6,35,36/
       INTEGER compallow_deepbask_to_deepbask 
     &   (ncompallow_deepbask_to_deepbask) 
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_deepbask_to_deepng   
     &   (ncompallow_deepbask_to_deepng  ) 
     &  /2,15,28,41/
       INTEGER compallow_deepbask_to_deepaxax 
     &   (ncompallow_deepbask_to_deepaxax) 
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_deepbask_to_deepLTS  
     &   (ncompallow_deepbask_to_deepLTS)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_deepbask_to_nontuftRS
     &   (ncompallow_deepbask_to_nontuftRS)
     &  /1,2,3,4,5,6,35,36/

       INTEGER compallow_deepng_to_tuftIB    
     &  (ncompallow_deepng_to_tuftIB   )
     & /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34/
       INTEGER compallow_deepng_to_tuftRS    
     &  (ncompallow_deepng_to_tuftRS   )
     & /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34/
       INTEGER compallow_deepng_to_nontuftRS 
     &  (ncompallow_deepng_to_nontuftRS)
     & /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34/
       INTEGER compallow_deepng_to_spinstell 
     &  (ncompallow_deepng_to_spinstell)
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_deepng_to_deepng    
     &  (ncompallow_deepng_to_deepng   )
     &  /2,15,28,41/
       INTEGER compallow_deepng_to_deepbask  
     &  (ncompallow_deepng_to_deepbask )
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/

       INTEGER compallow_deepLTS_to_suppyrRS
     &   (ncompallow_deepLTS_to_suppyrRS)
     & /14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
     &  31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,
     &  50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
     &  67,68/
       INTEGER compallow_deepLTS_to_suppyrFRB
     &   (ncompallow_deepLTS_to_suppyrFRB)
     & /14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
     &  31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,48,49,
     &  50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
     &  67,68/
       INTEGER compallow_deepLTS_to_supbask  
     &   (ncompallow_deepLTS_to_supbask)  
     & / 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
     &   47,48,49,50,51/ 
       INTEGER compallow_deepLTS_to_supaxax  
     &   (ncompallow_deepLTS_to_supaxax)  
     & / 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
     &   47,48,49,50,51/ 
       INTEGER compallow_deepLTS_to_supLTS   
     &   (ncompallow_deepLTS_to_supLTS)   
     & / 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
     &   47,48,49,50,51/ 
       INTEGER compallow_deepLTS_to_spinstell
     &   (ncompallow_deepLTS_to_spinstell)
     & /5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
     &  26,27,31,32,33,34,35,36,37,38,39,40,
     &  44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_deepLTS_to_tuftIB   
     &   (ncompallow_deepLTS_to_tuftIB)    
     & / 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
     &   29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,
     &   48,49,50,51,52,53,54,55/
       INTEGER compallow_deepLTS_to_tuftRS   
     &   (ncompallow_deepLTS_to_tuftRS)    
     & / 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
     &   29,30,31,32,33,34,38,39,40,41,42,43,44,45,46,47,
     &   48,49,50,51,52,53,54,55/
       INTEGER compallow_deepLTS_to_deepbask 
     &   (ncompallow_deepLTS_to_deepbask)  
     & /5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
     &  26,27,31,32,33,34,35,36,37,38,39,40,
     &  44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_deepLTS_to_deepaxax 
     &   (ncompallow_deepLTS_to_deepaxax)  
     & /5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
     &  26,27,31,32,33,34,35,36,37,38,39,40,
     &  44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_deepLTS_to_deepLTS  
     &   (ncompallow_deepLTS_to_deepLTS)   
     & /5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,24,25,
     &  26,27,31,32,33,34,35,36,37,38,39,40,
     &  44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_deepLTS_to_nontuftRS
     &   (ncompallow_deepLTS_to_nontuftRS) 
     & / 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
     &   29,30,31,32,33,34,38,39,40,41,42,43,44/ 

       INTEGER compallow_TCR_to_suppyrRS
     &   (ncompallow_TCR_to_suppyrRS)
     &  /45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
     &   61,62,63,64,65,66,67,68/
       INTEGER compallow_TCR_to_suppyrFRB
     &   (ncompallow_TCR_to_suppyrFRB)
     &  /45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
     &   61,62,63,64,65,66,67,68/
       INTEGER compallow_TCR_to_supbask  
     &   (ncompallow_TCR_to_supbask)  
     &  /2,3,4,15,16,17,28,29,30,41,42,43/
       INTEGER compallow_TCR_to_supng    
     &   (ncompallow_TCR_to_supng  )  
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_TCR_to_supaxax  
     &   (ncompallow_TCR_to_supaxax)  
     &  /2,3,4,15,16,17,28,29,30,41,42,43/
       INTEGER compallow_TCR_to_spinstell
     &   (ncompallow_TCR_to_spinstell)
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &   37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_TCR_to_tuftIB   
     &   (ncompallow_TCR_to_tuftIB)   
     &  / 47,48,49,50,51,52,53,54,55/
       INTEGER compallow_TCR_to_tuftRS   
     &   (ncompallow_TCR_to_tuftRS)   
     &  / 47,48,49,50,51,52,53,54,55/
       INTEGER compallow_TCR_to_deepbask 
     &   (ncompallow_TCR_to_deepbask) 
!    &  /2,3,4,15,16,17,28,29,30,41,42,43/
     &  /1,2,15,28,41/  ! soma & proximal dendrites
       INTEGER compallow_TCR_to_deepng   
     &   (ncompallow_TCR_to_deepng  ) 
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_TCR_to_deepaxax 
     &   (ncompallow_TCR_to_deepaxax) 
!    &  /2,3,4,15,16,17,28,29,30,41,42,43/
     &  /1,2,15,28,41/  ! soma & proximal dendrites
       INTEGER compallow_TCR_to_nRT      
     &   (ncompallow_TCR_to_nRT)      
     &  /2,3,4,15,16,17,28,29,30,41,42,43/
       INTEGER compallow_TCR_to_nontuftRS
     &   (ncompallow_TCR_to_nontuftRS)
     &  /40,41,42,43,44/

       INTEGER compallow_nRT_to_TCR
     &   (ncompallow_nRT_to_TCR)
     &  /1,2,15,28,41,54,67,80,93,106,119/
       INTEGER compallow_nRT_to_nRT
     &   (ncompallow_nRT_to_nRT)
     &  /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
     &   20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,
     &   36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,
     &   52,53/

        INTEGER compallow_nontuftRS_to_suppyrRS
     &    (ncompallow_nontuftRS_to_suppyrRS)
     &   / 41,42,43,44 /
        INTEGER compallow_nontuftRS_to_suppyrFRB
     &    (ncompallow_nontuftRS_to_suppyrFRB)
     &   / 41,42,43,44 /
        INTEGER compallow_nontuftRS_to_supbask  
     &    (ncompallow_nontuftRS_to_supbask)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
        INTEGER compallow_nontuftRS_to_supaxax  
     &    (ncompallow_nontuftRS_to_supaxax)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
        INTEGER compallow_nontuftRS_to_supLTS   
     &    (ncompallow_nontuftRS_to_supLTS)   
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
        INTEGER compallow_nontuftRS_to_spinstell
     &    (ncompallow_nontuftRS_to_spinstell)
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
        INTEGER compallow_nontuftRS_to_tuftIB   
     &    (ncompallow_nontuftRS_to_tuftIB)   
     &  / 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
     &   38,39,40,41,42,43,44,45,46,47/
        INTEGER compallow_nontuftRS_to_tuftRS   
     &    (ncompallow_nontuftRS_to_tuftRS)   
     &  / 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,
     &   38,39,40,41,42,43,44,45,46,47/
        INTEGER compallow_nontuftRS_to_deepbask 
     &    (ncompallow_nontuftRS_to_deepbask) 
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
        INTEGER compallow_nontuftRS_to_deepng   
     &    (ncompallow_nontuftRS_to_deepng  ) 
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
        INTEGER compallow_nontuftRS_to_deepaxax 
     &    (ncompallow_nontuftRS_to_deepaxax) 
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
        INTEGER compallow_nontuftRS_to_deepLTS  
     &    (ncompallow_nontuftRS_to_deepLTS)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
        INTEGER compallow_nontuftRS_to_TCR      
     &    (ncompallow_nontuftRS_to_TCR)      
     &  /  6,  7,  8,  9, 10, 11, 12, 13, 14,
     &    19, 20, 21, 22, 23, 24, 25, 26, 27,
     &    32, 33, 34, 35, 36, 37, 38, 39, 40,
     &    45, 46, 47, 48, 49, 50, 51, 52, 53,
     &    58, 59, 60, 61, 62, 63, 64, 65, 66,
     &    71, 72, 73, 74, 75, 76, 77, 78, 79,
     &    84, 85, 86, 87, 88, 89, 90, 91, 92,
     &    97, 98, 99,100,101,102,103,104,105,
     &   110,111,112,113,114,115,116,117,118,
     &   123,124,125,126,127,128,129,130,131/
        INTEGER compallow_nontuftRS_to_nRT      
     &    (ncompallow_nontuftRS_to_nRT)      
     & / 2,3,4,15,16,17,28,29,30,41,42,43/
        INTEGER compallow_nontuftRS_to_nontuftRS
     &    (ncompallow_nontuftRS_to_nontuftRS)
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &   37,38,39,40,41,42,43,44/


c Maps of synaptic connectivity.  For simplicity, all contacts
c only made to one compartment.  Axoaxonic cells forced to contact 
c initial axon segments; other compartments will be listed in arrays.
        INTEGER 
     & map_suppyrRS_to_suppyrRS(num_suppyrRS_to_suppyrRS,
     &                           num_suppyrRS),
     & map_suppyrRS_to_suppyrFRB(num_suppyrRS_to_suppyrFRB,
     &                           num_suppyrFRB),
     & map_suppyrRS_to_supbask(num_suppyrRS_to_supbask,  
     &                           num_supbask), 
     & map_suppyrRS_to_supng  (num_suppyrRS_to_supng  ,  
     &                           num_supng  ), 
     & map_suppyrRS_to_supaxax(num_suppyrRS_to_supaxax, 
     &                           num_supaxax),
     & map_suppyrRS_to_supLTS(num_suppyrRS_to_supLTS,   
     &                           num_supLTS),
     & map_suppyrRS_to_spinstell(num_suppyrRS_to_spinstell,
     &                           num_spinstell),
     & map_suppyrRS_to_tuftIB(num_suppyrRS_to_tuftIB,
     &                           num_tuftIB),  
     & map_suppyrRS_to_tuftRS(num_suppyrRS_to_tuftRS,
     &                           num_tuftRS), 
     & map_suppyrRS_to_deepbask(num_suppyrRS_to_deepbask,
     &                           num_deepbask), 
     & map_suppyrRS_to_deepaxax(num_suppyrRS_to_deepaxax,
     &                           num_deepaxax), 
     & map_suppyrRS_to_deepLTS(num_suppyrRS_to_deepLTS,
     &                           num_deepLTS), 
     & map_suppyrRS_to_nontuftrS(num_suppyrRS_to_nontuftRS,
     &                           num_nontuftRS), 
     & map_suppyrFRB_to_suppyrRS(num_suppyrFRB_to_suppyrRS,
     &                           num_suppyrRS) 
              INTEGER
     & map_suppyrFRB_to_suppyrFRB(num_suppyrFRB_to_suppyrFRB,
     &                           num_suppyrFRB),
     & map_suppyrFRB_to_supbask(num_suppyrFRB_to_supbask,
     &                           num_supbask), 
     & map_suppyrFRB_to_supng  (num_suppyrFRB_to_supng  ,
     &                           num_supng  ), 
     & map_suppyrFRB_to_supaxax(num_suppyrFRB_to_supaxax,
     &                           num_supaxax),
     & map_suppyrFRB_to_supLTS(num_suppyrFRB_to_supLTS,
     &                           num_supLTS),
     & map_suppyrFRB_to_spinstell(num_suppyrFRB_to_spinstell,
     &                           num_spinstell),
     & map_suppyrFRB_to_tuftIB(num_suppyrFRB_to_tuftIB,
     &                           num_tuftIB),  
     & map_suppyrFRB_to_tuftRS(num_suppyrFRB_to_tuftRS,
     &                           num_tuftRS),  
     & map_suppyrFRB_to_deepbask(num_suppyrFRB_to_deepbask,
     &                           num_deepbask),
     & map_suppyrFRB_to_deepaxax(num_suppyrFRB_to_deepaxax,
     &                           num_deepaxax),
     & map_suppyrFRB_to_deepLTS(num_suppyrFRB_to_deepLTS,
     &                           num_deepLTS), 
     & map_suppyrFRB_to_nontuftRS(num_suppyrFRB_to_nontuftRS,
     &                           num_nontuftRS),

     & map_supbask_to_suppyrRS(num_supbask_to_suppyrRS,
     &                           num_suppyrRS),  
     & map_supbask_to_suppyrFRB(num_supbask_to_suppyrFRB,
     &                           num_suppyrFRB),
     & map_supbask_to_supbask(num_supbask_to_supbask,
     &                           num_supbask), 
     & map_supbask_to_supng  (num_supbask_to_supng  ,
     &                           num_supng  ),  
     & map_supbask_to_supaxax(num_supbask_to_supaxax,
     &                           num_supaxax),
     & map_supbask_to_supLTS(num_supbask_to_supLTS,
     &                           num_supLTS),  
     & map_supbask_to_spinstell(num_supbask_to_spinstell,
     &                           num_spinstell)  
              INTEGER
     & map_supng_to_suppyrRS  (num_supng_to_suppyrRS ,
     &                           num_suppyrRS),
     & map_supng_to_suppyrFRB (num_supng_to_suppyrFRB,
     &                           num_suppyrFRB),
     & map_supng_to_nontuftRS (num_supng_to_nontuftRS,
     &                           num_nontuftRS),
     & map_supng_to_tuftIB    (num_supng_to_tuftIB   ,
     &                           num_tuftIB   ),
     & map_supng_to_tuftRS    (num_supng_to_tuftRS   ,
     &                           num_tuftRS   ),
     & map_supng_to_supng     (num_supng_to_supng    ,
     &                           num_supng    ),
     & map_supng_to_supbask   (num_supng_to_supbask  ,
     &                           num_supbask  ), 

     & map_supaxax_to_suppyrRS(num_supaxax_to_suppyrRS,
     &                           num_suppyrRS), 
     & map_supaxax_to_suppyrFRB(num_supaxax_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & map_supaxax_to_spinstell(num_supaxax_to_spinstell,
     &                           num_spinstell),
     & map_supaxax_to_tuftIB(num_supaxax_to_tuftIB,
     &                           num_tuftIB),  
     & map_supaxax_to_tuftRS(num_supaxax_to_tuftRS,
     &                           num_tuftRS), 
     & map_supaxax_to_nontuftRS(num_supaxax_to_nontuftRS,
     &                           num_nontuftRS), 
     & map_supLTS_to_suppyrRS(num_supLTS_to_suppyrRS,
     &                           num_suppyrRS),  
     & map_supLTS_to_suppyrFRB(num_supLTS_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & map_supLTS_to_supbask(num_supLTS_to_supbask,
     &                           num_supbask),  
     & map_supLTS_to_supaxax(num_supLTS_to_supaxax,
     &                           num_supaxax), 
     & map_supLTS_to_supLTS(num_supLTS_to_supLTS,
     &                           num_supLTS), 
     & map_supLTS_to_spinstell(num_supLTS_to_spinstell,
     &                           num_spinstell), 
     & map_supLTS_to_tuftIB(num_supLTS_to_tuftIB,
     &                           num_tuftIB),   
     & map_supLTS_to_tuftRS(num_supLTS_to_tuftRS,
     &                           num_tuftRS),  
     & map_supLTS_to_deepbask(num_supLTS_to_deepbask,
     &                           num_deepbask), 
     & map_supLTS_to_deepaxax(num_supLTS_to_deepaxax,
     &                           num_deepaxax), 
     & map_supLTS_to_deepLTS(num_supLTS_to_deepLTS,
     &                           num_deepLTS), 
     & map_supLTS_to_nontuftRS(num_supLTS_to_nontuftRS,
     &                           num_nontuftRS), 
     & map_spinstell_to_suppyrRS(num_spinstell_to_suppyrRS,
     &                           num_suppyrRS),
     & map_spinstell_to_suppyrFRB(num_spinstell_to_suppyrFRB,
     &                           num_suppyrFRB),
     & map_spinstell_to_supbask(num_spinstell_to_supbask,
     &                           num_supbask) 
               INTEGER
     & map_spinstell_to_supaxax(num_spinstell_to_supaxax,
     &                           num_supaxax),
     & map_spinstell_to_supLTS(num_spinstell_to_supLTS,
     &                           num_supLTS), 
     & map_spinstell_to_spinstell(num_spinstell_to_spinstell,
     &                           num_spinstell),
     & map_spinstell_to_tuftIB(num_spinstell_to_tuftIB,
     &                           num_tuftIB),  
     & map_spinstell_to_tuftRS(num_spinstell_to_tuftRS,
     &                           num_tuftRS), 
     & map_spinstell_to_deepbask(num_spinstell_to_deepbask,
     &                           num_deepbask), 
     & map_spinstell_to_deepng  (num_spinstell_to_deepng  ,
     &                           num_deepng  ), 
     & map_spinstell_to_deepaxax(num_spinstell_to_deepaxax,
     &                           num_deepaxax),
     & map_spinstell_to_deepLTS(num_spinstell_to_deepLTS,
     &                           num_deepLTS),
     & map_spinstell_to_nontuftRS(num_spinstell_to_nontuftRS,
     &                           num_nontuftRS),

     & map_tuftIB_to_suppyrRS(num_tuftIB_to_suppyrRS,
     &                           num_suppyrRS),   
     & map_tuftIB_to_suppyrFRB(num_tuftIB_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & map_tuftIB_to_supbask(num_tuftIB_to_supbask,
     &                           num_supbask),  
     & map_tuftIB_to_supaxax(num_tuftIB_to_supaxax,
     &                           num_supaxax), 
     & map_tuftIB_to_supLTS(num_tuftIB_to_supLTS,
     &                           num_supLTS), 
     & map_tuftIB_to_spinstell(num_tuftIB_to_spinstell,
     &                           num_spinstell), 
     & map_tuftIB_to_tuftIB(num_tuftIB_to_tuftIB,
     &                           num_tuftIB),   
     & map_tuftIB_to_tuftRS(num_tuftIB_to_tuftRS,
     &                           num_tuftRS),  
     & map_tuftIB_to_deepbask(num_tuftIB_to_deepbask,
     &                           num_deepbask), 
     & map_tuftIB_to_deepng  (num_tuftIB_to_deepng  ,
     &                           num_deepng  ), 
     & map_tuftIB_to_deepaxax(num_tuftIB_to_deepaxax,
     &                           num_deepaxax),  
     & map_tuftIB_to_deepLTS(num_tuftIB_to_deepLTS,
     &                           num_deepLTS),  
     & map_tuftIB_to_nontuftRS(num_tuftIB_to_nontuftRS,
     &                           num_nontuftRS), 
     & map_tuftRS_to_suppyrRS(num_tuftRS_to_suppyrRS,
     &                           num_suppyrRS), 
     & map_tuftRS_to_suppyrFRB(num_tuftRS_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & map_tuftRS_to_supbask(num_tuftRS_to_supbask,
     &                           num_supbask),  
     & map_tuftRS_to_supaxax(num_tuftRS_to_supaxax,
     &                           num_supaxax),   
     & map_tuftRS_to_supLTS(num_tuftRS_to_supLTS,
     &                           num_supLTS)     
            INTEGER
     & map_tuftRS_to_spinstell(num_tuftRS_to_spinstell,
     &                           num_spinstell), 
     & map_tuftRS_to_tuftIB(num_tuftRS_to_tuftIB,
     &                           num_tuftIB),   
     & map_tuftRS_to_tuftRS(num_tuftRS_to_tuftRS,
     &                           num_tuftRS),     
     & map_tuftRS_to_deepbask(num_tuftRS_to_deepbask,
     &                           num_deepbask),  
     & map_tuftRS_to_deepng  (num_tuftRS_to_deepng  ,
     &                           num_deepng  ),  
     & map_tuftRS_to_deepaxax(num_tuftRS_to_deepaxax,
     &                           num_deepaxax),   
     & map_tuftRS_to_deepLTS(num_tuftRS_to_deepLTS,
     &                           num_deepLTS),   
     & map_tuftRS_to_nontuftRS(num_tuftRS_to_nontuftRS,
     &                           num_nontuftRS),  
     & map_deepbask_to_spinstell(num_deepbask_to_spinstell,
     &                           num_spinstell), 
     & map_deepbask_to_tuftIB(num_deepbask_to_tuftIB,
     &                           num_tuftIB),   
     & map_deepbask_to_tuftRS(num_deepbask_to_tuftRS,
     &                           num_tuftRS),  
     & map_deepbask_to_deepbask(num_deepbask_to_deepbask,
     &                           num_deepbask), 
     & map_deepbask_to_deepng  (num_deepbask_to_deepng  ,
     &                           num_deepng  ), 
     & map_deepbask_to_deepaxax(num_deepbask_to_deepaxax,
     &                           num_deepaxax),  
     & map_deepbask_to_deepLTS(num_deepbask_to_deepLTS,
     &                           num_deepLTS)  
                INTEGER
     & map_deepbask_to_nontuftRS(num_deepbask_to_nontuftRS,
     &                           num_nontuftRS), 
     & map_deepng_to_tuftIB     (num_deepng_to_tuftIB     ,
     &                           num_tuftIB      ),
     & map_deepng_to_tuftRS     (num_deepng_to_tuftRS     ,
     &                           num_tuftRS      ),
     & map_deepng_to_nontuftRS  (num_deepng_to_nontuftRS  ,
     &                           num_nontuftRS   ),
     & map_deepng_to_spinstell  (num_deepng_to_spinstell  ,
     &                           num_spinstell   ),
     & map_deepng_to_deepng     (num_deepng_to_deepng     ,
     &                           num_deepng      ),
     & map_deepng_to_deepbask   (num_deepng_to_deepbask   ,
     &                           num_deepbask    ) 

                INTEGER
     & map_deepaxax_to_suppyrRS(num_deepaxax_to_suppyrRS,
     &                           num_suppyrRS), 
     & map_deepaxax_to_suppyrFRB(num_deepaxax_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & map_deepaxax_to_spinstell(num_deepaxax_to_spinstell,
     &                           num_spinstell),
     & map_deepaxax_to_tuftIB(num_deepaxax_to_tuftIB,
     &                           num_tuftIB), 
     & map_deepaxax_to_tuftRS(num_deepaxax_to_tuftRS,
     &                           num_tuftRS),    
     & map_deepaxax_to_nontuftRS(num_deepaxax_to_nontuftRS,
     &                           num_nontuftRS),
     & map_deepLTS_to_suppyrRS(num_deepLTS_to_suppyrRS,
     &                           num_suppyrRS), 
     & map_deepLTS_to_suppyrFRB(num_deepLTS_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & map_deepLTS_to_supbask(num_deepLTS_to_supbask,
     &                           num_supbask),  
     & map_deepLTS_to_supaxax(num_deepLTS_to_supaxax,
     &                           num_supaxax), 
     & map_deepLTS_to_supLTS(num_deepLTS_to_supLTS,
     &                           num_supLTS), 
     & map_deepLTS_to_spinstell(num_deepLTS_to_spinstell,
     &                           num_spinstell),
     & map_deepLTS_to_tuftIB(num_deepLTS_to_tuftIB,
     &                           num_tuftIB),  
     & map_deepLTS_to_tuftRS(num_deepLTS_to_tuftRS,
     &                            num_tuftRS), 
     & map_deepLTS_to_deepbask(num_deepLTS_to_deepbask,
     &                            num_deepbask), 
     & map_deepLTS_to_deepaxax(num_deepLTS_to_deepaxax,
     &                            num_deepaxax),  
     & map_deepLTS_to_deepLTS(num_deepLTS_to_deepLTS,
     &                            num_deepLTS),  
     & map_deepLTS_to_nontuftRS(num_deepLTS_to_nontuftRS,
     &                            num_nontuftRS), 
     & map_TCR_to_suppyrRS(num_TCR_to_suppyrRS,
     &                            num_suppyrRS),     
     & map_TCR_to_suppyrFRB(num_TCR_to_suppyrFRB,
     &                            num_suppyrFRB),   
     & map_TCR_to_supbask(num_TCR_to_supbask,
     &                            num_supbask)    
               INTEGER
     & map_TCR_to_supng  (num_TCR_to_supng  ,
     &                            num_supng  ),   
     & map_TCR_to_supaxax(num_TCR_to_supaxax,num_supaxax),   
     & map_TCR_to_spinstell(num_TCR_to_spinstell,num_spinstell),
     & map_TCR_to_tuftIB(num_TCR_to_tuftIB,num_tuftIB),  
     & map_TCR_to_tuftRS(num_TCR_to_tuftRS,num_tuftRS),    
     & map_TCR_to_deepbask(num_TCR_to_deepbask,num_deepbask), 
     & map_TCR_to_deepng  (num_TCR_to_deepng  ,num_deepng  ), 
     & map_TCR_to_deepaxax(num_TCR_to_deepaxax,num_deepaxax),
     & map_TCR_to_nRT(num_TCR_to_nRT,num_nRT),    
     & map_TCR_to_nontuftRS(num_TCR_to_nontuftRS,num_nontuftRS), 
     & map_nRT_to_TCR(num_nRT_to_TCR,num_TCR),      
     & map_nRT_to_nRT(num_nRT_to_nRT,num_nRT),     
     & map_nontuftRS_to_suppyrRS(num_nontuftRS_to_suppyrRS,
     &                             num_suppyrRS), 
     & map_nontuftRS_to_suppyrFRB(num_nontuftRS_to_suppyrFRB,
     &                             num_suppyrFRB),
     & map_nontuftRS_to_supbask(num_nontuftRS_to_supbask,
     &                             num_supbask), 
     & map_nontuftRS_to_supaxax(num_nontuftRS_to_supaxax,
     &                             num_supaxax),
     & map_nontuftRS_to_supLTS(num_nontuftRS_to_supLTS,
     &                             num_supLTS),  
     & map_nontuftRS_to_spinstell(num_nontuftRS_to_spinstell,
     &                             num_spinstell),
     & map_nontuftRS_to_tuftIB(num_nontuftRS_to_tuftIB,
     &                             num_tuftIB),  
     & map_nontuftRS_to_tuftRS(num_nontuftRS_to_tuftRS,
     &                             num_tuftRS),  
     & map_nontuftRS_to_deepbask(num_nontuftRS_to_deepbask,
     &                             num_deepbask), 
     & map_nontuftRS_to_deepng  (num_nontuftRS_to_deepng  ,
     &                             num_deepng  ), 
     & map_nontuftRS_to_deepaxax(num_nontuftRS_to_deepaxax,
     &                             num_deepaxax),
     & map_nontuftRS_to_deepLTS(num_nontuftRS_to_deepLTS,
     &                             num_deepLTS),
     & map_nontuftRS_to_TCR(num_nontuftRS_to_TCR,num_TCR),
     & map_nontuftRS_to_nRT(num_nontuftRS_to_nRT,num_nRT),  
     & map_nontuftRS_to_nontuftRS(num_nontuftRS_to_nontuftRS,
     &                             num_nontuftRS)

c Maps of synaptic compartments.  For simplicity, all contacts
c only made to one compartment.  Axoaxonic cells forced to contact 
c initial axon segments; other compartments will be listed in arrays.
        INTEGER 
     & com_suppyrRS_to_suppyrRS(num_suppyrRS_to_suppyrRS,
     &                           num_suppyrRS),
     & com_suppyrRS_to_suppyrFRB(num_suppyrRS_to_suppyrFRB,
     &                           num_suppyrFRB),
     & com_suppyrRS_to_supbask(num_suppyrRS_to_supbask,  
     &                           num_supbask), 
     & com_suppyrRS_to_supng  (num_suppyrRS_to_supng  ,  
     &                           num_supng  ), 
     & com_suppyrRS_to_supaxax(num_suppyrRS_to_supaxax, 
     &                           num_supaxax),
     & com_suppyrRS_to_supLTS(num_suppyrRS_to_supLTS,   
     &                           num_supLTS),
     & com_suppyrRS_to_spinstell(num_suppyrRS_to_spinstell,
     &                           num_spinstell),
     & com_suppyrRS_to_tuftIB(num_suppyrRS_to_tuftIB,
     &                           num_tuftIB),  
     & com_suppyrRS_to_tuftRS(num_suppyrRS_to_tuftRS,
     &                           num_tuftRS), 
     & com_suppyrRS_to_deepbask(num_suppyrRS_to_deepbask,
     &                           num_deepbask), 
     & com_suppyrRS_to_deepaxax(num_suppyrRS_to_deepaxax,
     &                           num_deepaxax), 
     & com_suppyrRS_to_deepLTS(num_suppyrRS_to_deepLTS,
     &                           num_deepLTS), 
     & com_suppyrRS_to_nontuftrS(num_suppyrRS_to_nontuftRS,
     &                           num_nontuftRS), 
     & com_suppyrFRB_to_suppyrRS(num_suppyrFRB_to_suppyrRS,
     &                           num_suppyrRS) 
              INTEGER
     & com_suppyrFRB_to_suppyrFRB(num_suppyrFRB_to_suppyrFRB,
     &                           num_suppyrFRB),
     & com_suppyrFRB_to_supbask(num_suppyrFRB_to_supbask,
     &                           num_supbask), 
     & com_suppyrFRB_to_supng  (num_suppyrFRB_to_supng  ,
     &                           num_supng  ), 
     & com_suppyrFRB_to_supaxax(num_suppyrFRB_to_supaxax,
     &                           num_supaxax),
     & com_suppyrFRB_to_supLTS(num_suppyrFRB_to_supLTS,
     &                           num_supLTS),
     & com_suppyrFRB_to_spinstell(num_suppyrFRB_to_spinstell,
     &                           num_spinstell),
     & com_suppyrFRB_to_tuftIB(num_suppyrFRB_to_tuftIB,
     &                           num_tuftIB),  
     & com_suppyrFRB_to_tuftRS(num_suppyrFRB_to_tuftRS,
     &                           num_tuftRS),  
     & com_suppyrFRB_to_deepbask(num_suppyrFRB_to_deepbask,
     &                           num_deepbask),
     & com_suppyrFRB_to_deepaxax(num_suppyrFRB_to_deepaxax,
     &                           num_deepaxax),
     & com_suppyrFRB_to_deepLTS(num_suppyrFRB_to_deepLTS,
     &                           num_deepLTS), 
     & com_suppyrFRB_to_nontuftRS(num_suppyrFRB_to_nontuftRS,
     &                           num_nontuftRS),
     & com_supbask_to_suppyrRS(num_supbask_to_suppyrRS,
     &                           num_suppyrRS),  
     & com_supbask_to_suppyrFRB(num_supbask_to_suppyrFRB,
     &                           num_suppyrFRB),
     & com_supbask_to_supbask(num_supbask_to_supbask,
     &                           num_supbask), 
     & com_supbask_to_supng  (num_supbask_to_supng  ,
     &                           num_supng  ), 
     & com_supbask_to_supaxax(num_supbask_to_supaxax,
     &                           num_supaxax),
     & com_supbask_to_supLTS(num_supbask_to_supLTS,
     &                           num_supLTS),  
     & com_supbask_to_spinstell(num_supbask_to_spinstell,
     &                           num_spinstell)  

          INTEGER
     & com_supng_to_suppyrRS  (num_supng_to_suppyrRS,
     &                         num_suppyrRS),
     & com_supng_to_suppyrFRB (num_supng_to_suppyrFRB,
     &                         num_suppyrFRB),
     & com_supng_to_nontuftRS (num_supng_to_nontuftRS,
     &                         num_nontuftRS),
     & com_supng_to_tuftIB    (num_supng_to_tuftIB   ,
     &                         num_tuftIB   ),
     & com_supng_to_tuftRS    (num_supng_to_tuftRS   ,
     &                         num_tuftRS   ),
     & com_supng_to_supng     (num_supng_to_supng    ,
     &                         num_supng    ),
     & com_supng_to_supbask   (num_supng_to_supbask  ,
     &                         num_supbask  ) 

          INTEGER
     & com_supaxax_to_suppyrRS(num_supaxax_to_suppyrRS,
     &                           num_suppyrRS), 
     & com_supaxax_to_suppyrFRB(num_supaxax_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & com_supaxax_to_spinstell(num_supaxax_to_spinstell,
     &                           num_spinstell)
           INTEGER
     & com_supaxax_to_tuftIB(num_supaxax_to_tuftIB,
     &                           num_tuftIB),  
     & com_supaxax_to_tuftRS(num_supaxax_to_tuftRS,
     &                           num_tuftRS), 
     & com_supaxax_to_nontuftRS(num_supaxax_to_nontuftRS,
     &                           num_nontuftRS), 
     & com_supLTS_to_suppyrRS(num_supLTS_to_suppyrRS,
     &                           num_suppyrRS),  
     & com_supLTS_to_suppyrFRB(num_supLTS_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & com_supLTS_to_supbask(num_supLTS_to_supbask,
     &                           num_supbask),  
     & com_supLTS_to_supaxax(num_supLTS_to_supaxax,
     &                           num_supaxax), 
     & com_supLTS_to_supLTS(num_supLTS_to_supLTS,
     &                           num_supLTS), 
     & com_supLTS_to_spinstell(num_supLTS_to_spinstell,
     &                           num_spinstell), 
     & com_supLTS_to_tuftIB(num_supLTS_to_tuftIB,
     &                           num_tuftIB),   
     & com_supLTS_to_tuftRS(num_supLTS_to_tuftRS,
     &                           num_tuftRS),  
     & com_supLTS_to_deepbask(num_supLTS_to_deepbask,
     &                           num_deepbask), 
     & com_supLTS_to_deepaxax(num_supLTS_to_deepaxax,
     &                           num_deepaxax), 
     & com_supLTS_to_deepLTS(num_supLTS_to_deepLTS,
     &                           num_deepLTS), 
     & com_supLTS_to_nontuftRS(num_supLTS_to_nontuftRS,
     &                           num_nontuftRS), 
     & com_spinstell_to_suppyrRS(num_spinstell_to_suppyrRS,
     &                           num_suppyrRS),
     & com_spinstell_to_suppyrFRB(num_spinstell_to_suppyrFRB,
     &                           num_suppyrFRB),
     & com_spinstell_to_supbask(num_spinstell_to_supbask,
     &                           num_supbask), 
     & com_spinstell_to_supaxax(num_spinstell_to_supaxax,
     &                           num_supaxax)
                INTEGER
     & com_spinstell_to_supLTS(num_spinstell_to_supLTS,
     &                           num_supLTS), 
     & com_spinstell_to_spinstell(num_spinstell_to_spinstell,
     &                           num_spinstell),
     & com_spinstell_to_tuftIB(num_spinstell_to_tuftIB,
     &                           num_tuftIB),  
     & com_spinstell_to_tuftRS(num_spinstell_to_tuftRS,
     &                           num_tuftRS), 
     & com_spinstell_to_deepbask(num_spinstell_to_deepbask,
     &                           num_deepbask), 
     & com_spinstell_to_deepng  (num_spinstell_to_deepng  ,
     &                           num_deepng  ), 
     & com_spinstell_to_deepaxax(num_spinstell_to_deepaxax,
     &                           num_deepaxax),
     & com_spinstell_to_deepLTS(num_spinstell_to_deepLTS,
     &                           num_deepLTS),
     & com_spinstell_to_nontuftRS(num_spinstell_to_nontuftRS,
     &                           num_nontuftRS),
     & com_tuftIB_to_suppyrRS(num_tuftIB_to_suppyrRS,
     &                           num_suppyrRS),   
     & com_tuftIB_to_suppyrFRB(num_tuftIB_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & com_tuftIB_to_supbask(num_tuftIB_to_supbask,
     &                           num_supbask),  
     & com_tuftIB_to_supaxax(num_tuftIB_to_supaxax,
     &                           num_supaxax), 
     & com_tuftIB_to_supLTS(num_tuftIB_to_supLTS,
     &                           num_supLTS), 
     & com_tuftIB_to_spinstell(num_tuftIB_to_spinstell,
     &                           num_spinstell), 
     & com_tuftIB_to_tuftIB(num_tuftIB_to_tuftIB,
     &                           num_tuftIB),   
     & com_tuftIB_to_tuftRS(num_tuftIB_to_tuftRS,
     &                           num_tuftRS),  
     & com_tuftIB_to_deepbask(num_tuftIB_to_deepbask,
     &                           num_deepbask), 
     & com_tuftIB_to_deepng  (num_tuftIB_to_deepng  ,
     &                           num_deepng  ), 
     & com_tuftIB_to_deepaxax(num_tuftIB_to_deepaxax,
     &                           num_deepaxax),  
     & com_tuftIB_to_deepLTS(num_tuftIB_to_deepLTS,
     &                           num_deepLTS),  
     & com_tuftIB_to_nontuftRS(num_tuftIB_to_nontuftRS,
     &                           num_nontuftRS) 
              INTEGER
     & com_tuftRS_to_suppyrRS(num_tuftRS_to_suppyrRS,
     &                           num_suppyrRS), 
     & com_tuftRS_to_suppyrFRB(num_tuftRS_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & com_tuftRS_to_supbask(num_tuftRS_to_supbask,
     &                           num_supbask),  
     & com_tuftRS_to_supaxax(num_tuftRS_to_supaxax,
     &                           num_supaxax),   
     & com_tuftRS_to_supLTS(num_tuftRS_to_supLTS,
     &                           num_supLTS),     
     & com_tuftRS_to_spinstell(num_tuftRS_to_spinstell,
     &                           num_spinstell), 
     & com_tuftRS_to_tuftIB(num_tuftRS_to_tuftIB,
     &                           num_tuftIB),   
     & com_tuftRS_to_tuftRS(num_tuftRS_to_tuftRS,
     &                           num_tuftRS),     
     & com_tuftRS_to_deepbask(num_tuftRS_to_deepbask,
     &                           num_deepbask),  
     & com_tuftRS_to_deepng  (num_tuftRS_to_deepng  ,
     &                           num_deepng  ),  
     & com_tuftRS_to_deepaxax(num_tuftRS_to_deepaxax,
     &                           num_deepaxax),   
     & com_tuftRS_to_deepLTS(num_tuftRS_to_deepLTS,
     &                           num_deepLTS),   
     & com_tuftRS_to_nontuftRS(num_tuftRS_to_nontuftRS,
     &                           num_nontuftRS),  
     & com_deepbask_to_spinstell(num_deepbask_to_spinstell,
     &                           num_spinstell), 
     & com_deepbask_to_tuftIB(num_deepbask_to_tuftIB,
     &                           num_tuftIB),   
     & com_deepbask_to_tuftRS(num_deepbask_to_tuftRS,
     &                           num_tuftRS),  
     & com_deepbask_to_deepbask(num_deepbask_to_deepbask,
     &                           num_deepbask), 
     & com_deepbask_to_deepng  (num_deepbask_to_deepng  ,
     &                           num_deepng  ), 
     & com_deepbask_to_deepaxax(num_deepbask_to_deepaxax,
     &                           num_deepaxax),  
     & com_deepbask_to_deepLTS(num_deepbask_to_deepLTS,
     &                           num_deepLTS),  
     & com_deepbask_to_nontuftRS(num_deepbask_to_nontuftRS,
     &                           num_nontuftRS) 
            INTEGER
     & com_deepng_to_tuftIB     (num_deepng_to_tuftIB    ,
     &                           num_tuftIB      ),
     & com_deepng_to_tuftRS     (num_deepng_to_tuftRS    ,
     &                           num_tuftRS      ),
     & com_deepng_to_nontuftRS  (num_deepng_to_nontuftRS ,
     &                           num_nontuftRS   ),
     & com_deepng_to_spinstell  (num_deepng_to_spinstell ,
     &                           num_spinstell   ),
     & com_deepng_to_deepng     (num_deepng_to_deepng    ,
     &                           num_deepng      ),
     & com_deepng_to_deepbask   (num_deepng_to_deepbask  ,
     &                           num_deepbask    ) 
            INTEGER
     & com_deepaxax_to_suppyrRS(num_deepaxax_to_suppyrRS,
     &                           num_suppyrRS), 
     & com_deepaxax_to_suppyrFRB(num_deepaxax_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & com_deepaxax_to_spinstell(num_deepaxax_to_spinstell,
     &                           num_spinstell),
     & com_deepaxax_to_tuftIB(num_deepaxax_to_tuftIB,
     &                           num_tuftIB), 
     & com_deepaxax_to_tuftRS(num_deepaxax_to_tuftRS,
     &                           num_tuftRS),    
     & com_deepaxax_to_nontuftRS(num_deepaxax_to_nontuftRS,
     &                           num_nontuftRS),
     & com_deepLTS_to_suppyrRS(num_deepLTS_to_suppyrRS,
     &                           num_suppyrRS), 
     & com_deepLTS_to_suppyrFRB(num_deepLTS_to_suppyrFRB,
     &                           num_suppyrFRB), 
     & com_deepLTS_to_supbask(num_deepLTS_to_supbask,
     &                           num_supbask),  
     & com_deepLTS_to_supaxax(num_deepLTS_to_supaxax,
     &                           num_supaxax), 
     & com_deepLTS_to_supLTS(num_deepLTS_to_supLTS,
     &                           num_supLTS), 
     & com_deepLTS_to_spinstell(num_deepLTS_to_spinstell,
     &                           num_spinstell),
     & com_deepLTS_to_tuftIB(num_deepLTS_to_tuftIB,
     &                           num_tuftIB),  
     & com_deepLTS_to_tuftRS(num_deepLTS_to_tuftRS,
     &                            num_tuftRS), 
     & com_deepLTS_to_deepbask(num_deepLTS_to_deepbask,
     &                            num_deepbask), 
     & com_deepLTS_to_deepaxax(num_deepLTS_to_deepaxax,
     &                            num_deepaxax),  
     & com_deepLTS_to_deepLTS(num_deepLTS_to_deepLTS,
     &                            num_deepLTS)  
           INTEGER
     & com_deepLTS_to_nontuftRS(num_deepLTS_to_nontuftRS,
     &                            num_nontuftRS), 
     & com_TCR_to_suppyrRS(num_TCR_to_suppyrRS,
     &                            num_suppyrRS),     
     & com_TCR_to_suppyrFRB(num_TCR_to_suppyrFRB,
     &                            num_suppyrFRB),   
     & com_TCR_to_supbask(num_TCR_to_supbask,
     &                            num_supbask),    
     & com_TCR_to_supng  (num_TCR_to_supng  ,
     &                            num_supng  ),    
     & com_TCR_to_supaxax(num_TCR_to_supaxax,num_supaxax),   
     & com_TCR_to_spinstell(num_TCR_to_spinstell,num_spinstell),
     & com_TCR_to_tuftIB(num_TCR_to_tuftIB,num_tuftIB),  
     & com_TCR_to_tuftRS(num_TCR_to_tuftRS,num_tuftRS),    
     & com_TCR_to_deepbask(num_TCR_to_deepbask,num_deepbask), 
     & com_TCR_to_deepng  (num_TCR_to_deepng  ,num_deepng  ), 
     & com_TCR_to_deepaxax(num_TCR_to_deepaxax,num_deepaxax),
     & com_TCR_to_nRT(num_TCR_to_nRT,num_nRT),    
     & com_TCR_to_nontuftRS(num_TCR_to_nontuftRS,num_nontuftRS), 
     & com_nRT_to_TCR(num_nRT_to_TCR,num_TCR),      
     & com_nRT_to_nRT(num_nRT_to_nRT,num_nRT),     
     & com_nontuftRS_to_suppyrRS(num_nontuftRS_to_suppyrRS,
     &                             num_suppyrRS), 
     & com_nontuftRS_to_suppyrFRB(num_nontuftRS_to_suppyrFRB,
     &                             num_suppyrFRB),
     & com_nontuftRS_to_supbask(num_nontuftRS_to_supbask,
     &                             num_supbask), 
     & com_nontuftRS_to_supaxax(num_nontuftRS_to_supaxax,
     &                             num_supaxax),
     & com_nontuftRS_to_supLTS(num_nontuftRS_to_supLTS,
     &                             num_supLTS),  
     & com_nontuftRS_to_spinstell(num_nontuftRS_to_spinstell,
     &                             num_spinstell),
     & com_nontuftRS_to_tuftIB(num_nontuftRS_to_tuftIB,
     &                             num_tuftIB)  
              INTEGER
     & com_nontuftRS_to_tuftRS(num_nontuftRS_to_tuftRS,
     &                             num_tuftRS),  
     & com_nontuftRS_to_deepbask(num_nontuftRS_to_deepbask,
     &                             num_deepbask), 
     & com_nontuftRS_to_deepng  (num_nontuftRS_to_deepng  ,
     &                             num_deepng  ), 
     & com_nontuftRS_to_deepaxax(num_nontuftRS_to_deepaxax,
     &                             num_deepaxax),
     & com_nontuftRS_to_deepLTS(num_nontuftRS_to_deepLTS,
     &                             num_deepLTS),
     & com_nontuftRS_to_TCR(num_nontuftRS_to_TCR,num_TCR),
     & com_nontuftRS_to_nRT(num_nontuftRS_to_nRT,num_nRT),  
     & com_nontuftRS_to_nontuftRS(num_nontuftRS_to_nontuftRS,
     &                             num_nontuftRS)

c Entries in gjtable are cell a, compart. of cell a with gj,
c  cell b, compart. of cell b with gj; entries not repeated,
c which means that, for given cell being integrated, table
c must be searched through cols. 1 and 3.
       integer gjtable_suppyrRS(totaxgj_suppyrRS,4),
     &   gjtable_suppyrFRB(totaxgj_suppyrFRB,4),
     &   gjtable_suppyr   (totaxgj_suppyr,4),
! gjtable_suppyr for suppyrRS/suppyrFRB gj, with RS cell
! in col. 1 and FRB cell in col. 3
     &   gjtable_supbask  (totSDgj_supbask,4),
     &   gjtable_supng    (totSDgj_supng  ,4),
     &   gjtable_supaxax  (1              ,4),
     &   gjtable_supLTS   (totSDgj_supLTS,4),
     &   gjtable_spinstell(totaxgj_spinstell,4),
     &   gjtable_tuftIB   (totaxgj_tuftIB,4),
     &   gjtable_tuftRS   (totaxgj_tuftRS,4),
     &   gjtable_tuft     (totaxgj_tuft,4),
! gjtable_tuft for tuftIB/tuftRS gj, with IB cell
! in col. 1 and RS cell in col. 3.
     &   gjtable_nontuftRS(totaxgj_nontuftRS,4),
     &   gjtable_deepbask (totSDgj_deepbask,4),
     &   gjtable_deepng   (totSDgj_deepng  ,4),
     &   gjtable_deepaxax (1               ,4),
     &   gjtable_deepLTS  (totSDgj_deepLTS,4),
     &   gjtable_TCR      (totaxgj_TCR,4),
     &   gjtable_nRT      (totSDgj_nRT,4) 

c define compartments on which gj can form
       INTEGER
     &table_axgjcompallow_suppyrRS(num_axgjcompallow_suppyrRS)
     &          /74/,
c    &          /73/, ! 28 Nov. 2005, move proximally, to get more inhib. control.
     &table_axgjcompallow_suppyrFRB(num_axgjcompallow_suppyrFRB)
     &          /74/,
c    &          /73/, ! 28 Nov. 2005, move proximally
c Ectopics to superficial pyr. cells then go to #72, see
c   supergj.f
     &table_SDgjcompallow_supbask  (num_SDgjcompallow_supbask  )
     &          /3,4,16,17,29,30,42,43/,
     &table_SDgjcompallow_supng    (num_SDgjcompallow_supng    )
     &          /3,4,16,17,29,30,42,43/,
     &table_SDgjcompallow_supLTS   (num_SDgjcompallow_supLTS   )
     &          /3,4,16,17,29,30,42,43/,
     &table_axgjcompallow_spinstell(num_axgjcompallow_spinstell)
     &          /59/,
c Ectopics to spiny stellates then go to #57
     &table_axgjcompallow_tuftIB   (num_axgjcompallow_tuftIB   )
     &          /61/,
     &table_axgjcompallow_tuftRS   (num_axgjcompallow_tuftRS   )
     &          /61/,
c Ectopics to tufted pyr. cells then go to #60
     &table_axgjcompallow_nontuftRS(num_axgjcompallow_nontuftRS)
     &          /50/,
c Ectopics to nontufted deep pyr. cells then to #48
     &table_SDgjcompallow_deepbask (num_SDgjcompallow_deepbask )
     &          /3,4,16,17,29,30,42,43/,
     &table_SDgjcompallow_deepng   (num_SDgjcompallow_deepng   )
     &          /3,4,16,17,29,30,42,43/,
     &table_SDgjcompallow_deepLTS  (num_SDgjcompallow_deepLTS  )
     &          /3,4,16,17,29,30,42,43/,
     &table_axgjcompallow_TCR      (num_axgjcompallow_TCR      )
     &          /137/,
c Ectopics to TCR cells to #135
     &table_SDgjcompallow_nRT      (num_SDgjcompallow_nRT      )
     &          /3,4,16,17,29,30,42,43/

       integer ectr_suppyrRS, ectr_suppyrFRB, ectr_supbask,
     &  ectr_supng, ectr_deepng,
     &  ectr_supaxax, ectr_supLTS, ectr_spinstell,
     &  ectr_tuftIB, ectr_tuftRS, ectr_nontuftRS,
     &  ectr_deepbask, ectr_deepaxax, ectr_deepLTS,
     &  ectr_TCR, ectr_nRT

!      real*8 field_1mm_suppyrRS, field_1mm_suppyrFRB,
!    &  field_1mm_nontuftRS, field_1mm_tuftIB, field_1mm_tuftRS
!      real*8 field_2mm_suppyrRS, field_2mm_suppyrFRB,
!    &  field_2mm_nontuftRS, field_2mm_tuftIB, field_2mm_tuftRS
       real*8 field_1mm, field_2mm ! scalars to pass to subroutines
       real*8 field_1mm_local(1), field_2mm_local(1)  ! for mpi
       real*8 field_1mm_global(numnodes), field_2mm_global(numnodes) ! for mpi
       real*8 field_1mm_tot, field_2mm_tot  ! sums of global vectors

c Define tables used for computing dexp & GABA-B timecourse:
c dexptablesmall(i) = dexp(-z), i = int (z*1000.), 0<=z<=5.
c dexptablebig  (i) = dexp(-z), i = int (z*10.), 0<=z<=100.
        double precision:: dexptablesmall(0:5000)
        double precision::  dexptablebig  (0:1000)
        double precision:: otis_table (0:50000)
! if how_often = 50 and dt = .002, then otis_table structure
! corresponds to time steps of 0.1 ms, and it gives 5 s of data.

        real*8 noisepe_tuftIB  ! noisepe_tuftIB_save defined as parameter above
        real*8 noisepe_tuftRS  ! noisepe_tuftRS_save defined as parameter above

c Define arrays, constants, for voltages, applied currents,
c synaptic conductances, random numbers, etc.

       double precision::
     &  V_suppyrRS  (numcomp_suppyrRS, num_suppyrRS),
     &  V_suppyrFRB (numcomp_suppyrFRB,num_suppyrFRB), 
     &  V_supbask   (numcomp_supbask,  num_supbask),  
     &  V_supng     (numcomp_supng  ,  num_supng  ),  
     &  V_supaxax   (numcomp_supaxax,  num_supaxax), 
     &  V_supLTS    (numcomp_supLTS,   num_supLTS), 
     &  V_spinstell (numcomp_spinstell,num_spinstell),
     &  V_tuftIB    (numcomp_tuftIB,   num_tuftIB),  
     &  V_tuftRS    (numcomp_tuftRS,   num_tuftRS), 
     &  V_nontuftRS (numcomp_nontuftRS,num_nontuftRS),
     &  V_deepbask  (numcomp_deepbask, num_deepbask),
     &  V_deepng    (numcomp_deepng  , num_deepng  ),
     &  V_deepaxax  (numcomp_deepaxax, num_deepaxax),
     &  V_deepLTS   (numcomp_deepLTS,  num_deepLTS),
     &  V_TCR       (numcomp_TCR,      num_TCR),   
     &  V_nRT       (numcomp_nRT,      num_nRT) 

       double precision::
     &  curr_suppyrRS   (numcomp_suppyrRS, num_suppyrRS),
     &  curr_suppyrFRB  (numcomp_suppyrFRB,num_suppyrFRB), 
     &  curr_supbask    (numcomp_supbask,  num_supbask),  
     &  curr_supng      (numcomp_supng  ,  num_supng  ),  
     &  curr_supaxax    (numcomp_supaxax,  num_supaxax), 
     &  curr_supLTS     (numcomp_supLTS,   num_supLTS), 
     &  curr_spinstell  (numcomp_spinstell,num_spinstell),
     &  curr_tuftIB     (numcomp_tuftIB,   num_tuftIB),  
     &  curr_tuftRS     (numcomp_tuftRS,   num_tuftRS), 
     &  curr_nontuftRS  (numcomp_nontuftRS,num_nontuftRS),
     &  curr_deepbask   (numcomp_deepbask, num_deepbask),
     &  curr_deepng     (numcomp_deepng  , num_deepng  ),
     &  curr_deepaxax   (numcomp_deepaxax, num_deepaxax),
     &  curr_deepLTS    (numcomp_deepLTS,  num_deepLTS),
     &  curr_TCR        (numcomp_TCR,      num_TCR),   
     &  curr_nRT        (numcomp_nRT,      num_nRT) 

       double precision::
     & gAMPA_suppyrRS   (numcomp_suppyrRS, num_suppyrRS),
     & gAMPA_suppyrFRB  (numcomp_suppyrFRB,num_suppyrFRB), 
     & gAMPA_supbask    (numcomp_supbask,  num_supbask),  
     & gAMPA_supng      (numcomp_supng  ,  num_supng  ),  
     & gAMPA_supaxax    (numcomp_supaxax,  num_supaxax), 
     & gAMPA_supLTS     (numcomp_supLTS,   num_supLTS), 
     & gAMPA_spinstell  (numcomp_spinstell,num_spinstell),
     & gAMPA_tuftIB     (numcomp_tuftIB,   num_tuftIB),  
     & gAMPA_tuftRS     (numcomp_tuftRS,   num_tuftRS), 
     & gAMPA_nontuftRS  (numcomp_nontuftRS,num_nontuftRS),
     & gAMPA_deepbask   (numcomp_deepbask, num_deepbask),
     & gAMPA_deepng     (numcomp_deepng  , num_deepng  ),
     & gAMPA_deepaxax   (numcomp_deepaxax, num_deepaxax),
     & gAMPA_deepLTS    (numcomp_deepLTS,  num_deepLTS),
     & gAMPA_TCR        (numcomp_TCR,      num_TCR),   
     & gAMPA_nRT        (numcomp_nRT,      num_nRT) 

       double precision::
     & gNMDA_suppyrRS   (numcomp_suppyrRS, num_suppyrRS),
     & gNMDA_suppyrFRB  (numcomp_suppyrFRB,num_suppyrFRB), 
     & gNMDA_supbask    (numcomp_supbask,  num_supbask),  
     & gNMDA_supng      (numcomp_supng  ,  num_supng  ),  
     & gNMDA_supaxax    (numcomp_supaxax,  num_supaxax), 
     & gNMDA_supLTS     (numcomp_supLTS,   num_supLTS), 
     & gNMDA_spinstell  (numcomp_spinstell,num_spinstell),
     & gNMDA_tuftIB     (numcomp_tuftIB,   num_tuftIB),  
     & gNMDA_tuftRS     (numcomp_tuftRS,   num_tuftRS), 
     & gNMDA_nontuftRS  (numcomp_nontuftRS,num_nontuftRS),
     & gNMDA_deepbask   (numcomp_deepbask, num_deepbask),
     & gNMDA_deepng     (numcomp_deepng  , num_deepng  ),
     & gNMDA_deepaxax   (numcomp_deepaxax, num_deepaxax),
     & gNMDA_deepLTS    (numcomp_deepLTS,  num_deepLTS),
     & gNMDA_TCR        (numcomp_TCR,      num_TCR),   
     & gNMDA_nRT        (numcomp_nRT,      num_nRT) 

       double precision::
     & gGABA_A_suppyrRS (numcomp_suppyrRS, num_suppyrRS),
     & gGABA_A_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB), 
     & gGABA_A_supbask  (numcomp_supbask,  num_supbask),  
     & gGABA_A_supng    (numcomp_supng  ,  num_supng  ),  
     & gGABA_A_supaxax  (numcomp_supaxax,  num_supaxax), 
     & gGABA_A_supLTS   (numcomp_supLTS,   num_supLTS), 
     & gGABA_A_spinstell(numcomp_spinstell,num_spinstell),
     & gGABA_A_tuftIB   (numcomp_tuftIB,   num_tuftIB),  
     & gGABA_A_tuftRS   (numcomp_tuftRS,   num_tuftRS), 
     & gGABA_A_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     & gGABA_A_deepbask (numcomp_deepbask, num_deepbask),
     & gGABA_A_deepng   (numcomp_deepng  , num_deepng  ),
     & gGABA_A_deepaxax (numcomp_deepaxax, num_deepaxax),
     & gGABA_A_deepLTS  (numcomp_deepLTS,  num_deepLTS),
     & gGABA_A_TCR      (numcomp_TCR,      num_TCR),   
     & gGABA_A_nRT      (numcomp_nRT,      num_nRT) 

       double precision::
     & gGABA_B_suppyrRS (numcomp_suppyrRS, num_suppyrRS),
     & gGABA_B_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB), 
     & gGABA_B_spinstell(numcomp_spinstell,num_spinstell),
     & gGABA_B_tuftIB   (numcomp_tuftIB,   num_tuftIB),  
     & gGABA_B_tuftRS   (numcomp_tuftRS,   num_tuftRS), 
     & gGABA_B_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     & gGABA_B_TCR      (numcomp_TCR,      num_TCR),   
     & gGABA_B_nRT      (numcomp_nRT,      num_nRT) 

! define membrane and Ca state variables that must be passed
! to subroutines
       real*8  chi_suppyrRS(numcomp_suppyrRS,num_suppyrRS)
       real*8  mnaf_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     & mnap_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x hnaf_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mkdr_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mka_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x hka_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mk2_suppyrRS(numcomp_suppyrRS,num_suppyrRS), 
     x hk2_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mkm_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mkc_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mkahp_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mcat_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x hcat_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mcal_suppyrRS(numcomp_suppyrRS,num_suppyrRS),
     x mar_suppyrRS(numcomp_suppyrRS,num_suppyrRS)

       real*8  chi_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB)
       real*8  mnaf_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     & mnap_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x hnaf_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x mkdr_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x mka_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x hka_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x mk2_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB), 
     x hk2_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x mkm_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x mkc_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x mkahp_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x mcat_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x hcat_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x mcal_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB),
     x mar_suppyrFRB(numcomp_suppyrFRB,num_suppyrFRB)

       real*8  chi_supbask (numcomp_supbask ,num_supbask )
       real*8  mnaf_supbask (numcomp_supbask ,num_supbask ),
     & mnap_supbask (numcomp_supbask ,num_supbask ),
     x hnaf_supbask (numcomp_supbask ,num_supbask ),
     x mkdr_supbask (numcomp_supbask ,num_supbask ),
     x mka_supbask (numcomp_supbask ,num_supbask ),
     x hka_supbask (numcomp_supbask ,num_supbask ),
     x mk2_supbask (numcomp_supbask ,num_supbask ), 
     x hk2_supbask (numcomp_supbask ,num_supbask ),
     x mkm_supbask (numcomp_supbask ,num_supbask ),
     x mkc_supbask (numcomp_supbask ,num_supbask ),
     x mkahp_supbask (numcomp_supbask ,num_supbask ),
     x mcat_supbask (numcomp_supbask ,num_supbask ),
     x hcat_supbask (numcomp_supbask ,num_supbask ),
     x mcal_supbask (numcomp_supbask ,num_supbask ),
     x mar_supbask (numcomp_supbask ,num_supbask )

       real*8  chi_supng (numcomp_supng ,num_supng )
       real*8  mnaf_supng (numcomp_supng ,num_supng ),
     & mnap_supng (numcomp_supng ,num_supng ),
     x hnaf_supng (numcomp_supng ,num_supng ),
     x mkdr_supng (numcomp_supng ,num_supng ),
     x mka_supng (numcomp_supng ,num_supng ),
     x hka_supng (numcomp_supng ,num_supng ),
     x mk2_supng (numcomp_supng ,num_supng ), 
     x hk2_supng (numcomp_supng ,num_supng ),
     x mkm_supng (numcomp_supng ,num_supng ),
     x mkc_supng (numcomp_supng ,num_supng ),
     x mkahp_supng (numcomp_supng ,num_supng ),
     x mcat_supng (numcomp_supng ,num_supng ),
     x hcat_supng (numcomp_supng ,num_supng ),
     x mcal_supng (numcomp_supng ,num_supng ),
     x mar_supng (numcomp_supng ,num_supng )

       real*8  chi_supaxax (numcomp_supaxax ,num_supaxax )
       real*8  mnaf_supaxax (numcomp_supaxax ,num_supaxax ),
     & mnap_supaxax (numcomp_supaxax ,num_supaxax ),
     x hnaf_supaxax (numcomp_supaxax ,num_supaxax ),
     x mkdr_supaxax (numcomp_supaxax ,num_supaxax ),
     x mka_supaxax (numcomp_supaxax ,num_supaxax ),
     x hka_supaxax (numcomp_supaxax ,num_supaxax ),
     x mk2_supaxax (numcomp_supaxax ,num_supaxax ), 
     x hk2_supaxax (numcomp_supaxax ,num_supaxax ),
     x mkm_supaxax (numcomp_supaxax ,num_supaxax ),
     x mkc_supaxax (numcomp_supaxax ,num_supaxax ),
     x mkahp_supaxax (numcomp_supaxax ,num_supaxax ),
     x mcat_supaxax (numcomp_supaxax ,num_supaxax ),
     x hcat_supaxax (numcomp_supaxax ,num_supaxax ),
     x mcal_supaxax (numcomp_supaxax ,num_supaxax ),
     x mar_supaxax (numcomp_supaxax ,num_supaxax )

       real*8  chi_supLTS(numcomp_supLTS,num_supLTS)
       real*8  mnaf_supLTS(numcomp_supLTS,num_supLTS),
     & mnap_supLTS(numcomp_supLTS,num_supLTS),
     x hnaf_supLTS(numcomp_supLTS,num_supLTS),
     x mkdr_supLTS(numcomp_supLTS,num_supLTS),
     x mka_supLTS(numcomp_supLTS,num_supLTS),
     x hka_supLTS(numcomp_supLTS,num_supLTS),
     x mk2_supLTS(numcomp_supLTS,num_supLTS), 
     x hk2_supLTS(numcomp_supLTS,num_supLTS),
     x mkm_supLTS(numcomp_supLTS,num_supLTS),
     x mkc_supLTS(numcomp_supLTS,num_supLTS),
     x mkahp_supLTS(numcomp_supLTS,num_supLTS),
     x mcat_supLTS(numcomp_supLTS,num_supLTS),
     x hcat_supLTS(numcomp_supLTS,num_supLTS),
     x mcal_supLTS(numcomp_supLTS,num_supLTS),
     x mar_supLTS(numcomp_supLTS,num_supLTS)

      real*8  chi_spinstell(numcomp_spinstell,num_spinstell)
      real*8  mnaf_spinstell(numcomp_spinstell,num_spinstell),
     & mnap_spinstell(numcomp_spinstell,num_spinstell),
     x hnaf_spinstell(numcomp_spinstell,num_spinstell),
     x mkdr_spinstell(numcomp_spinstell,num_spinstell),
     x mka_spinstell(numcomp_spinstell,num_spinstell),
     x hka_spinstell(numcomp_spinstell,num_spinstell),
     x mk2_spinstell(numcomp_spinstell,num_spinstell), 
     x hk2_spinstell(numcomp_spinstell,num_spinstell),
     x mkm_spinstell(numcomp_spinstell,num_spinstell),
     x mkc_spinstell(numcomp_spinstell,num_spinstell),
     x mkahp_spinstell(numcomp_spinstell,num_spinstell),
     x mcat_spinstell(numcomp_spinstell,num_spinstell),
     x hcat_spinstell(numcomp_spinstell,num_spinstell),
     x mcal_spinstell(numcomp_spinstell,num_spinstell),
     x mar_spinstell(numcomp_spinstell,num_spinstell)


       real*8  chi_tuftIB(numcomp_tuftIB,num_tuftIB)
       real*8  mnaf_tuftIB(numcomp_tuftIB,num_tuftIB),
     & mnap_tuftIB(numcomp_tuftIB,num_tuftIB),
     x hnaf_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mkdr_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mka_tuftIB(numcomp_tuftIB,num_tuftIB),
     x hka_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mk2_tuftIB(numcomp_tuftIB,num_tuftIB), 
     x hk2_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mkm_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mkc_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mkahp_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mcat_tuftIB(numcomp_tuftIB,num_tuftIB),
     x hcat_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mcal_tuftIB(numcomp_tuftIB,num_tuftIB),
     x mar_tuftIB(numcomp_tuftIB,num_tuftIB)

       real*8  chi_tuftRS(numcomp_tuftRS,num_tuftRS)
       real*8  mnaf_tuftRS(numcomp_tuftRS,num_tuftRS),
     & mnap_tuftRS(numcomp_tuftRS,num_tuftRS),
     x hnaf_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mkdr_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mka_tuftRS(numcomp_tuftRS,num_tuftRS),
     x hka_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mk2_tuftRS(numcomp_tuftRS,num_tuftRS), 
     x hk2_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mkm_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mkc_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mkahp_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mcat_tuftRS(numcomp_tuftRS,num_tuftRS),
     x hcat_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mcal_tuftRS(numcomp_tuftRS,num_tuftRS),
     x mar_tuftRS(numcomp_tuftRS,num_tuftRS)

       real*8  chi_nontuftRS(numcomp_nontuftRS,num_nontuftRS)
       real*8  mnaf_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     & mnap_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x hnaf_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mkdr_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mka_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x hka_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mk2_nontuftRS(numcomp_nontuftRS,num_nontuftRS), 
     x hk2_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mkm_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mkc_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mkahp_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mcat_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x hcat_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mcal_nontuftRS(numcomp_nontuftRS,num_nontuftRS),
     x mar_nontuftRS(numcomp_nontuftRS,num_nontuftRS)

       real*8  chi_deepbask(numcomp_deepbask,num_deepbask)
       real*8  mnaf_deepbask(numcomp_deepbask,num_deepbask),
     & mnap_deepbask(numcomp_deepbask,num_deepbask),
     x hnaf_deepbask(numcomp_deepbask,num_deepbask),
     x mkdr_deepbask(numcomp_deepbask,num_deepbask),
     x mka_deepbask(numcomp_deepbask,num_deepbask),
     x hka_deepbask(numcomp_deepbask,num_deepbask),
     x mk2_deepbask(numcomp_deepbask,num_deepbask), 
     x hk2_deepbask(numcomp_deepbask,num_deepbask),
     x mkm_deepbask(numcomp_deepbask,num_deepbask),
     x mkc_deepbask(numcomp_deepbask,num_deepbask),
     x mkahp_deepbask(numcomp_deepbask,num_deepbask),
     x mcat_deepbask(numcomp_deepbask,num_deepbask),
     x hcat_deepbask(numcomp_deepbask,num_deepbask),
     x mcal_deepbask(numcomp_deepbask,num_deepbask),
     x mar_deepbask(numcomp_deepbask,num_deepbask)

       real*8  chi_deepng(numcomp_deepng,num_deepng)
       real*8  mnaf_deepng(numcomp_deepng,num_deepng),
     & mnap_deepng(numcomp_deepng,num_deepng),
     x hnaf_deepng(numcomp_deepng,num_deepng),
     x mkdr_deepng(numcomp_deepng,num_deepng),
     x mka_deepng(numcomp_deepng,num_deepng),
     x hka_deepng(numcomp_deepng,num_deepng),
     x mk2_deepng(numcomp_deepng,num_deepng), 
     x hk2_deepng(numcomp_deepng,num_deepng),
     x mkm_deepng(numcomp_deepng,num_deepng),
     x mkc_deepng(numcomp_deepng,num_deepng),
     x mkahp_deepng(numcomp_deepng,num_deepng),
     x mcat_deepng(numcomp_deepng,num_deepng),
     x hcat_deepng(numcomp_deepng,num_deepng),
     x mcal_deepng(numcomp_deepng,num_deepng),
     x mar_deepng(numcomp_deepng,num_deepng)

       real*8  chi_deepaxax(numcomp_deepaxax,num_deepaxax)
       real*8  mnaf_deepaxax(numcomp_deepaxax,num_deepaxax),
     & mnap_deepaxax(numcomp_deepaxax,num_deepaxax),
     x hnaf_deepaxax(numcomp_deepaxax,num_deepaxax),
     x mkdr_deepaxax(numcomp_deepaxax,num_deepaxax),
     x mka_deepaxax(numcomp_deepaxax,num_deepaxax),
     x hka_deepaxax(numcomp_deepaxax,num_deepaxax),
     x mk2_deepaxax(numcomp_deepaxax,num_deepaxax), 
     x hk2_deepaxax(numcomp_deepaxax,num_deepaxax),
     x mkm_deepaxax(numcomp_deepaxax,num_deepaxax),
     x mkc_deepaxax(numcomp_deepaxax,num_deepaxax),
     x mkahp_deepaxax(numcomp_deepaxax,num_deepaxax),
     x mcat_deepaxax(numcomp_deepaxax,num_deepaxax),
     x hcat_deepaxax(numcomp_deepaxax,num_deepaxax),
     x mcal_deepaxax(numcomp_deepaxax,num_deepaxax),
     x mar_deepaxax(numcomp_deepaxax,num_deepaxax)

       real*8  chi_deepLTS(numcomp_deepLTS,num_deepLTS)
       real*8  mnaf_deepLTS(numcomp_deepLTS,num_deepLTS),
     & mnap_deepLTS(numcomp_deepLTS,num_deepLTS),
     x hnaf_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkdr_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mka_deepLTS(numcomp_deepLTS,num_deepLTS),
     x hka_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mk2_deepLTS(numcomp_deepLTS,num_deepLTS), 
     x hk2_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkm_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkc_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkahp_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mcat_deepLTS(numcomp_deepLTS,num_deepLTS),
     x hcat_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mcal_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mar_deepLTS(numcomp_deepLTS,num_deepLTS)

       real*8  chi_tcr(numcomp_tcr,num_tcr)
       real*8  mnaf_tcr(numcomp_tcr,num_tcr),
     & mnap_tcr(numcomp_tcr,num_tcr),
     x hnaf_tcr(numcomp_tcr,num_tcr),
     x mkdr_tcr(numcomp_tcr,num_tcr),
     x mka_tcr(numcomp_tcr,num_tcr),
     x hka_tcr(numcomp_tcr,num_tcr),
     x mk2_tcr(numcomp_tcr,num_tcr), 
     x hk2_tcr(numcomp_tcr,num_tcr),
     x mkm_tcr(numcomp_tcr,num_tcr),
     x mkc_tcr(numcomp_tcr,num_tcr),
     x mkahp_tcr(numcomp_tcr,num_tcr),
     x mcat_tcr(numcomp_tcr,num_tcr),
     x hcat_tcr(numcomp_tcr,num_tcr),
     x mcal_tcr(numcomp_tcr,num_tcr),
     x mar_tcr(numcomp_tcr,num_tcr)

       real*8  chi_nRT(numcomp_nRT,num_nRT)
       real*8  mnaf_nRT(numcomp_nRT,num_nRT),
     & mnap_nRT(numcomp_nRT,num_nRT),
     x hnaf_nRT(numcomp_nRT,num_nRT),
     x mkdr_nRT(numcomp_nRT,num_nRT),
     x mka_nRT(numcomp_nRT,num_nRT),
     x hka_nRT(numcomp_nRT,num_nRT),
     x mk2_nRT(numcomp_nRT,num_nRT), 
     x hk2_nRT(numcomp_nRT,num_nRT),
     x mkm_nRT(numcomp_nRT,num_nRT),
     x mkc_nRT(numcomp_nRT,num_nRT),
     x mkahp_nRT(numcomp_nRT,num_nRT),
     x mcat_nRT(numcomp_nRT,num_nRT),
     x hcat_nRT(numcomp_nRT,num_nRT),
     x mcal_nRT(numcomp_nRT,num_nRT),
     x mar_nRT(numcomp_nRT,num_nRT)

       double precision
     &    ranvec_suppyrRS  (num_suppyrRS),
     &    ranvec_suppyrFRB (num_suppyrFRB), 
     &    ranvec_supbask   (num_supbask),  
     &    ranvec_supng     (num_supng  ),  
     &    ranvec_supaxax   (num_supaxax), 
     &    ranvec_supLTS    (num_supLTS), 
     &    ranvec_spinstell (num_spinstell),
     &    ranvec_tuftIB    (num_tuftIB),  
     &    ranvec_tuftRS    (num_tuftRS), 
     &    ranvec_nontuftRS (num_nontuftRS),
     &    ranvec_deepbask  (num_deepbask),
     &    ranvec_deepng    (num_deepng  ),
     &    ranvec_deepaxax  (num_deepaxax),
     &    ranvec_deepLTS   (num_deepLTS),
     &    ranvec_TCR       (num_TCR),   
     &    ranvec_nRT       (num_nRT),
     &    seed /137.d0/

c Define arrays for distal axon voltages which will be shared
c between nodes, and for axonal sites of possible gj
         double precision::
     &  distal_axon_suppyrRS  (maxcellspernode),
     &  ldistal_axon_suppyrRS (num_suppyrRS), ! use for outtime
     &      gj_axon_suppyrRS  (maxcellspernode),
     &  distal_axon_suppyrFRB (maxcellspernode),
     &  ldistal_axon_suppyrFRB(num_suppyrFRB),
     &      gj_axon_suppyrFRB (maxcellspernode),
     &  distal_axon_supbask   (maxcellspernode),
     &  ldistal_axon_supbask  (num_supbask),
     &      gj_axon_supbask   (maxcellspernode),
     &  distal_axon_supng     (maxcellspernode),
     &  ldistal_axon_supng    (num_supng  ),
     &      gj_axon_supng     (maxcellspernode),
     &  distal_axon_supaxax   (maxcellspernode),
     &  ldistal_axon_supaxax  (num_supaxax),
     &      gj_axon_supaxax   (maxcellspernode),
     &  distal_axon_supLTS    (maxcellspernode),
     &  ldistal_axon_supLTS   (num_supLTS),
     &      gj_axon_supLTS    (maxcellspernode) 
         double precision::
     &  distal_axon_spinstell (maxcellspernode),
     &  ldistal_axon_spinstell(num_spinstell),
     &      gj_axon_spinstell (maxcellspernode),
     &  distal_axon_tuftIB    (maxcellspernode),
     &  ldistal_axon_tuftIB   (num_tuftIB),
     &      gj_axon_tuftIB    (maxcellspernode),
     &  distal_axon_tuftRS    (maxcellspernode),
     &  ldistal_axon_tuftRS   (num_tuftRS),
     &      gj_axon_tuftRS    (maxcellspernode),
     &  distal_axon_nontuftRS (maxcellspernode),
     &  ldistal_axon_nontuftRS(num_nontuftRS),
     &      gj_axon_nontuftRS (maxcellspernode),
     &  distal_axon_deepbask  (maxcellspernode),
     &  ldistal_axon_deepbask (num_deepbask),
     &      gj_axon_deepbask  (maxcellspernode),
     &  distal_axon_deepng    (maxcellspernode),
     &  ldistal_axon_deepng   (num_deepng  ),
     &      gj_axon_deepng    (maxcellspernode),
     &  distal_axon_deepaxax  (maxcellspernode),
     &  ldistal_axon_deepaxax (num_deepaxax),
     &      gj_axon_deepaxax  (maxcellspernode),
     &  distal_axon_deepLTS   (maxcellspernode),
     &  ldistal_axon_deepLTS  (num_deepLTS),
     &      gj_axon_deepLTS   (maxcellspernode),
     &  distal_axon_TCR       (maxcellspernode),
     &  ldistal_axon_TCR      (num_TCR),
     &      gj_axon_TCR       (maxcellspernode),
     &  distal_axon_nRT       (maxcellspernode),
     &  ldistal_axon_nRT      (num_nRT),
     &      gj_axon_nRT       (maxcellspernode),
!    Communication will be complicated, however, because - say - a tuftIB
!   will have to communicate only the tuftIB axons it has integrated.
     &  distal_axon_global    (numnodes  * maxcellspernode),
     &      gj_axon_global    (numnodes  * maxcellspernode)
! distal_axon_global will be concatenation of individual
! distal_axon vectors       

! define arrays for axonal voltges, needed for mixed gj
! PROBABLY should organize this as per distal_axon...
! ALSO must consider how to handle interneuron gj, in case interneurons
! spread across multiple nodes...
c        double precision ::
c    &    vax_suppyrRS (num_suppyrRS), vax_suppyrFRB (num_suppyrFRB),
c    &    vax_tuftIB (num_tuftIB), vax_tuftRS (num_tuftRS)

         double precision::
     &  outtime_suppyrRS  (5000, num_suppyrRS),
     &  outtime_suppyrFRB (5000, num_suppyrFRB),
     &  outtime_supbask   (5000, num_supbask), 
     &  outtime_supng     (5000, num_supng  ), 
     &  outtime_supaxax   (5000, num_supaxax), 
     &  outtime_supLTS    (5000, num_supLTS),   
     &  outtime_spinstell (5000, num_spinstell), 
     &  outtime_tuftIB    (5000, num_tuftIB), 
     &  outtime_tuftRS    (5000, num_tuftRS),  
     &  outtime_nontuftRS (5000, num_nontuftRS),
     &  outtime_deepbask  (5000, num_deepbask),
     &  outtime_deepng    (5000, num_deepng  ),
     &  outtime_deepaxax  (5000, num_deepaxax),
     &  outtime_deepLTS   (5000, num_deepLTS), 
     &  outtime_TCR       (5000, num_TCR),      
     &  outtime_nRT       (5000, num_nRT)       

         INTEGER
     &  outctr_suppyrRS  (num_suppyrRS), 
     &  outctr_suppyrFRB (num_suppyrFRB),
     &  outctr_supbask   (num_supbask), 
     &  outctr_supng     (num_supng  ), 
     &  outctr_supaxax   (num_supaxax),
     &  outctr_supLTS    (num_supLTS),
     &  outctr_spinstell (num_spinstell),
     &  outctr_tuftIB    (num_tuftIB), 
     &  outctr_tuftRS    (num_tuftRS),
     &  outctr_nontuftRS (num_nontuftRS),
     &  outctr_deepbask  (num_deepbask),
     &  outctr_deepng    (num_deepng  ),
     &  outctr_deepaxax  (num_deepaxax),
     &  outctr_deepLTS   (num_deepLTS),
     &  outctr_TCR       (num_TCR), 
     &  outctr_nRT       (num_nRT)

        CHARACTER(LEN=9) nodecell(0:numnodes-1) ! will define which cell type is to be handled by each node

        INTEGER place(0:numnodes-1)  ! this will define whether a node is 1st, 2nd... in the set of nodes
! used by a given type of cell

        integer initialize, firstcell, lastcell ! used in integration calls 
        integer ictr, ioffset

       REAL*8 gettime, time1, time2, time, timtot
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20), z3, z4, z3a, z4a, z5, z6, z7
       REAL*8 z10, z11, z12, z13, z14, z10a, z10b
       INTEGER i, j, k, L, k0, m

       double precision scale_tuftIB_gNaP(61) 
       double precision scale_tuftIB_gKM(61), Mshift ! for shifting gKM rate functions.
       double precision scale_tuftIB_gCaL (num_tuftIB)
       double precision gCaL_tuftIB(numcomp_tuftIB, num_tuftIB)
! declare this gCaL here, because there are problems in integration routine in
! making sure this conductance gets saved
       double precision rel_axonshift_tuftIB, rel_axonshift_suppyrRS,
     &    rel_axonshift_suppyrFRB

        double precision scale_suppyrFRB_gNaP(74)

c START EXECUTION PHASE
          include 'mpif.h'
          call mpi_init (info)
          call mpi_comm_rank(mpi_comm_world, thisno, info)
          call mpi_comm_size(mpi_comm_world, nodes , info)
          time1 = gettime()

          do i = 1, 74 
           scale_suppyrFRB_gNaP(i) = 0.65d0
          end do

c Define gCaL scaling for tuftIB - depends on cell rather than compartment
       call durand(seed,num_tuftIB,ranvec_tuftIB)
          do L = 1, num_tuftIB
           scale_tuftIB_gCaL(L) =
c    &   0.20d0 + 0.10d0 * ranvec_tuftIB(L) 
     &   0.35d0 + 0.10d0 * ranvec_tuftIB(L) 
c    &   1.60d0 + 0.40d0 * ranvec_tuftIB(L) 
c    &   0.00d0 + 0.00d0 * ranvec_tuftIB(L) 
          end do

c Define variable for shifting gKM rate functions in tuftIB
          Mshift =  0.d0

c Define gKM scaling for tufted IB pyramids
         do i = 1, 55  ! soma  & dendrites
c          scale_tuftIB_gKM(i) = 1.00d0
!          scale_tuftIB_gKM(i) = 0.30d0
c          scale_tuftIB_gKM(i) = 0.40d0
c          scale_tuftIB_gKM(i) = 0.10d0
           scale_tuftIB_gKM(i) = 0.05d0
           scale_tuftIB_gNaP(i) = 0.0d0
         end do
         do i = 56, 61  ! axon
c          scale_tuftIB_gKM(i) = 1.00d0
!          scale_tuftIB_gKM(i) = 0.50d0
c          scale_tuftIB_gKM(i) = 0.40d0
           scale_tuftIB_gKM(i) = 0.10d0
c          scale_tuftIB_gKM(i) = 0.05d0
c          scale_tuftIB_gNaP(i) = 1.0d0
           scale_tuftIB_gNaP(i) = 0.0d0
         end do

c Define which cell type is handled by each processor
           nodecell(0) = 'suppyrRS '
           nodecell(1) = 'suppyrRS '
           nodecell(2) = 'suppyrRS '
           nodecell(3) = 'suppyrRS '
           nodecell(4) = 'suppyrFRB'
           nodecell(5) = 'supbask  '
           nodecell(6) = 'supng    '
           nodecell(7) = 'supaxax  '
           nodecell(8) = 'supLTS   '
           nodecell(9) = 'spinstell'
           nodecell(10) = 'tuftIB   '
           nodecell(11)= 'tuftIB   '
           nodecell(12)= 'tuftIB   '
           nodecell(13)= 'tuftIB   '
           nodecell(14)= 'tuftIB   '
           nodecell(15)= 'tuftIB   '
c          nodecell(16)= 'tuftIB   '
c          nodecell(17)= 'tuftIB   '
c          nodecell(18)= 'tuftRS   '
           nodecell(16)= 'nontuftRS'
           nodecell(17)= 'nontuftRS'
           nodecell(18)= 'deepbask '
           nodecell(19)= 'deepng   '
           nodecell(20)= 'deepaxax '
           nodecell(21)= 'deepLTS  '
           nodecell(22)= 'TCR      '
           nodecell(23)= 'nRT      '
          if (thisno.eq.0) then
            do i = 0, numnodes - 1
              write(6,786) i, nodecell(i)
786           format(i5,a10)
            end do
          end if

c Define "rank" of nodes assigned to each cell-type - will
c be used in figuring out how to partition the cells.
           place( 0) = 1  ! suppyrRS: 1
           place( 1) = 2  ! suppyrRS: 2
           place( 2) = 3  ! suppyrRS: 3
           place( 3) = 4  ! suppyrRS: 4
           place( 4) = 1  ! suppyrFRB 
           place( 5) = 1  ! supbask   
           place( 6) = 1  ! supng     
           place( 7) = 1  ! supaxax   
           place( 8) = 1  ! supLTS    
           place( 9) = 1  ! spinstell 
           place(10) = 1  ! tuftIB: 1 
           place(11) = 2  ! tuftIB: 2 
           place(12) = 3  ! tuftIB: 3 
           place(13) = 4  ! tuftIB: 4 
           place(14) = 5  ! tuftIB: 5 
c          place(15) = 6  ! tuftIB: 6 
c          place(16) = 7  ! tuftIB: 7   
c          place(17) = 8  ! tuftIB: 8   
           place(15) = 1  ! tuftRS      
           place(16) = 1  ! nontuftRS: 1
           place(17) = 2  ! nontuftRS: 2
           place(18) = 1  ! deepbask    
           place(19) = 1  ! deepng      
           place(20) = 1  ! deepaxax    
           place(21) = 1  ! deepLTS     
           place(22) = 1  ! TCR         
           place(23) = 1  ! nRT         

         do i = 1, 5000
           do j = 1, num_suppyrRS
        outtime_suppyrRS(i,j)             = -1.d5
           end do ! j
           do j = 1, num_suppyrFRB
        outtime_suppyrFRB(i,j)            = -1.d5 
           end do ! j
           do j = 1, num_supbask  
        outtime_supbask(i,j)              = -1.d5
           end do ! j
           do j = 1, num_supng    
        outtime_supng  (i,j)              = -1.d5
           end do ! j
           do j = 1, num_supaxax  
        outtime_supaxax(i,j)              = -1.d5
           end do ! j
           do j = 1, num_supLTS   
        outtime_supLTS(i,j)               = -1.d5
           end do ! j
           do j = 1, num_spinstell
        outtime_spinstell(i,j)            = -1.d5 
           end do ! j
           do j = 1, num_tuftIB   
        outtime_tuftIB(i,j)               = -1.d5
           end do ! j
           do j = 1, num_tuftRS   
        outtime_tuftRS(i,j)               = -1.d5
           end do ! j
           do j = 1, num_nontuftRS   
        outtime_nontuftRS(i,j)            = -1.d5
           end do ! j
           do j = 1, num_deepbask    
        outtime_deepbask(i,j)             = -1.d5
           end do ! j
           do j = 1, num_deepng      
        outtime_deepng  (i,j)             = -1.d5
           end do ! j
           do j = 1, num_deepaxax    
        outtime_deepaxax(i,j)             = -1.d5
           end do ! j
           do j = 1, num_deepLTS     
        outtime_deepLTS(i,j)              = -1.d5
           end do ! j
           do j = 1, num_TCR         
        outtime_TCR(i,j)                  = -1.d5
           end do ! j
           do j = 1, num_nRT         
        outtime_nRT(i,j)                  = -1.d5
           end do ! j
         end do ! do i

c         timtot = 1000.d0 
c         timtot = 1250.d0
          timtot = 5000.d0
c         timtot = 0.d0

c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
          call dexptablebig_setup   (dexptablebig)
          call otis_table_setup (otis_table,how_often,dt)

c Compartments contacted by axoaxonic interneurons are IS only
          do i = 1, num_suppyrRS 
          do j = 1, num_supaxax_to_suppyrRS 
             com_supaxax_to_suppyrRS (j,i) = 69
          end do
          end do
          do i = 1, num_suppyrFRB
          do j = 1, num_supaxax_to_suppyrFRB
             com_supaxax_to_suppyrFRB(j,i) = 69
          end do
          end do
          do i = 1, num_spinstell
          do j = 1, num_supaxax_to_spinstell
             com_supaxax_to_spinstell(j,i) = 54
          end do
          end do
          do i = 1, num_tuftIB   
          do j = 1, num_supaxax_to_tuftIB   
             com_supaxax_to_tuftIB   (j,i) = 56
          end do
          end do
          do i = 1, num_tuftRS   
          do j = 1, num_supaxax_to_tuftRS   
             com_supaxax_to_tuftRS   (j,i) = 56
          end do
          end do
          do i = 1, num_nontuftRS   
          do j = 1, num_supaxax_to_nontuftRS   
             com_supaxax_to_nontuftRS   (j,i) = 45
          end do
          end do
          do i = 1, num_suppyrRS    
          do j = 1, num_deepaxax_to_suppyrRS    
             com_deepaxax_to_suppyrRS    (j,i) = 69
          end do
          end do
          do i = 1, num_suppyrFRB   
          do j = 1, num_deepaxax_to_suppyrFRB   
             com_deepaxax_to_suppyrFRB   (j,i) = 69
          end do
          end do
          do i = 1, num_spinstell   
          do j = 1, num_deepaxax_to_spinstell   
             com_deepaxax_to_spinstell   (j,i) = 54
          end do
          end do
          do i = 1, num_tuftIB      
          do j = 1, num_deepaxax_to_tuftIB      
             com_deepaxax_to_tuftIB      (j,i) = 56 
          end do
          end do
          do i = 1, num_tuftRS      
          do j = 1, num_deepaxax_to_tuftRS      
             com_deepaxax_to_tuftRS      (j,i) = 56 
          end do
          end do
          do i = 1, num_nontuftRS      
          do j = 1, num_deepaxax_to_nontuftRS      
             com_deepaxax_to_nontuftRS      (j,i) = 45 
          end do
          end do
c End section on making axoaxonic cells connect to IS's

c Construct synaptic connectivity tables
                display = 0

          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_suppyrRS,           
     &     map_suppyrRS_to_suppyrRS,
     &     num_suppyrRS_to_suppyrRS,    display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_suppyrFRB,          
     &     map_suppyrRS_to_suppyrFRB,
     &     num_suppyrRS_to_suppyrFRB,   display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_supbask,            
     &     map_suppyrRS_to_supbask,  
     &     num_suppyrRS_to_supbask,     display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_supng  ,            
     &     map_suppyrRS_to_supng  ,  
     &     num_suppyrRS_to_supng  ,     display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_supaxax,            
     &     map_suppyrRS_to_supaxax,  
     &     num_suppyrRS_to_supaxax,     display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_supLTS,             
     &     map_suppyrRS_to_supLTS,   
     &     num_suppyrRS_to_supLTS,      display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_spinstell,          
     &     map_suppyrRS_to_spinstell,
     &     num_suppyrRS_to_spinstell,   display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_tuftIB,             
     &     map_suppyrRS_to_tuftIB,   
     &     num_suppyrRS_to_tuftIB,      display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_tuftRS,             
     &     map_suppyrRS_to_tuftRS,   
     &     num_suppyrRS_to_tuftRS,      display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_deepbask,           
     &     map_suppyrRS_to_deepbask, 
     &     num_suppyrRS_to_deepbask,    display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_deepaxax,           
     &     map_suppyrRS_to_deepaxax, 
     &     num_suppyrRS_to_deepaxax,    display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_deepLTS,            
     &     map_suppyrRS_to_deepLTS,  
     &     num_suppyrRS_to_deepLTS,     display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrRS, num_nontuftRS,          
     &     map_suppyrRS_to_nontuftRS,
     &     num_suppyrRS_to_nontuftRS,   display)

          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_suppyrRS,           
     &     map_suppyrFRB_to_suppyrRS, 
     &     num_suppyrFRB_to_suppyrRS,    display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_suppyrFRB,          
     &     map_suppyrFRB_to_suppyrFRB,
     &     num_suppyrFRB_to_suppyrFRB,   display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_supbask,            
     &     map_suppyrFRB_to_supbask,  
     &     num_suppyrFRB_to_supbask,     display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_supng  ,            
     &     map_suppyrFRB_to_supng  ,  
     &     num_suppyrFRB_to_supng  ,     display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_supaxax,            
     &     map_suppyrFRB_to_supaxax,  
     &     num_suppyrFRB_to_supaxax,     display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_supLTS,             
     &     map_suppyrFRB_to_supLTS,   
     &     num_suppyrFRB_to_supLTS,      display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_spinstell,          
     &     map_suppyrFRB_to_spinstell,
     &     num_suppyrFRB_to_spinstell,   display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_tuftIB,             
     &     map_suppyrFRB_to_tuftIB,   
     &     num_suppyrFRB_to_tuftIB,      display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_tuftRS,             
     &     map_suppyrFRB_to_tuftRS,   
     &     num_suppyrFRB_to_tuftRS,      display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_deepbask,           
     &     map_suppyrFRB_to_deepbask, 
     &     num_suppyrFRB_to_deepbask,    display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_deepaxax,           
     &     map_suppyrFRB_to_deepaxax, 
     &     num_suppyrFRB_to_deepaxax,    display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_deepLTS,            
     &     map_suppyrFRB_to_deepLTS,  
     &     num_suppyrFRB_to_deepLTS,     display)
          CALL synaptic_map_construct (thisno,
     &     num_suppyrFRB, num_nontuftRS,          
     &     map_suppyrFRB_to_nontuftRS,
     &     num_suppyrFRB_to_nontuftRS,   display)

          CALL synaptic_map_construct (thisno,
     &     num_supbask, num_suppyrRS,           
     &     map_supbask_to_suppyrRS, 
     &     num_supbask_to_suppyrRS,   display)
          CALL synaptic_map_construct (thisno,
     &     num_supbask, num_suppyrFRB,          
     &     map_supbask_to_suppyrFRB,
     &     num_supbask_to_suppyrFRB,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supbask, num_supbask,            
     &     map_supbask_to_supbask,  
     &     num_supbask_to_supbask,    display)
          CALL synaptic_map_construct (thisno,
     &     num_supbask, num_supng  ,            
     &     map_supbask_to_supng  ,  
     &     num_supbask_to_supng  ,    display)
          CALL synaptic_map_construct (thisno,
     &     num_supbask, num_supaxax,            
     &     map_supbask_to_supaxax,  
     &     num_supbask_to_supaxax,    display)
          CALL synaptic_map_construct (thisno,
     &     num_supbask, num_supLTS,             
     &     map_supbask_to_supLTS,   
     &     num_supbask_to_supLTS,     display)
          CALL synaptic_map_construct (thisno,
     &     num_supbask, num_spinstell,          
     &     map_supbask_to_spinstell,
     &     num_supbask_to_spinstell,  display)

          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_suppyrRS ,          
     &     map_supng_to_suppyrRS ,
     &     num_supng_to_suppyrRS ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_suppyrFRB,          
     &     map_supng_to_suppyrFRB,
     &     num_supng_to_suppyrFRB,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_nontuftRS,          
     &     map_supng_to_nontuftRS,
     &     num_supng_to_nontuftRS,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_tuftIB   ,          
     &     map_supng_to_tuftIB   ,
     &     num_supng_to_tuftIB   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_tuftRS   ,          
     &     map_supng_to_tuftRS   ,
     &     num_supng_to_tuftRS   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_supng    ,          
     &     map_supng_to_supng    ,
     &     num_supng_to_supng    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_supbask  ,          
     &     map_supng_to_supbask  ,
     &     num_supng_to_supbask  ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_supaxax, num_suppyrRS,           
     &     map_supaxax_to_suppyrRS, 
     &     num_supaxax_to_suppyrRS,   display)
          CALL synaptic_map_construct (thisno,
     &     num_supaxax, num_suppyrFRB,          
     &     map_supaxax_to_suppyrFRB,
     &     num_supaxax_to_suppyrFRB,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supaxax, num_spinstell,          
     &     map_supaxax_to_spinstell,
     &     num_supaxax_to_spinstell,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supaxax, num_tuftIB,             
     &     map_supaxax_to_tuftIB,   
     &     num_supaxax_to_tuftIB,     display)
          CALL synaptic_map_construct (thisno,
     &     num_supaxax, num_tuftRS,             
     &     map_supaxax_to_tuftRS,   
     &     num_supaxax_to_tuftRS,     display)
          CALL synaptic_map_construct (thisno,
     &     num_supaxax, num_nontuftRS,             
     &     map_supaxax_to_nontuftRS,   
     &     num_supaxax_to_nontuftRS,  display)

          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_suppyrRS,              
     &     map_supLTS_to_suppyrRS,    
     &     num_supLTS_to_suppyrRS ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_suppyrFRB,             
     &     map_supLTS_to_suppyrFRB,   
     &     num_supLTS_to_suppyrFRB,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_supbask,               
     &     map_supLTS_to_supbask,     
     &     num_supLTS_to_supbask,    display)
          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_supaxax,               
     &     map_supLTS_to_supaxax,     
     &     num_supLTS_to_supaxax,    display)
          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_supLTS,                
     &     map_supLTS_to_supLTS,      
     &     num_supLTS_to_supLTS,     display)
          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_spinstell,             
     &     map_supLTS_to_spinstell,   
     &     num_supLTS_to_spinstell,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_tuftIB,                
     &     map_supLTS_to_tuftIB,      
     &     num_supLTS_to_tuftIB   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_tuftRS,                
     &     map_supLTS_to_tuftRS,      
     &     num_supLTS_to_tuftRS   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_deepbask,              
     &     map_supLTS_to_deepbask,    
     &     num_supLTS_to_deepbask ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_deepaxax,              
     &     map_supLTS_to_deepaxax,    
     &     num_supLTS_to_deepaxax ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_deepLTS,               
     &     map_supLTS_to_deepLTS,     
     &     num_supLTS_to_deepLTS,    display)
          CALL synaptic_map_construct (thisno,
     &     num_supLTS,  num_nontuftRS,             
     &     map_supLTS_to_nontuftRS,   
     &     num_supLTS_to_nontuftRS,  display)

          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_suppyrRS,              
     &     map_spinstell_to_suppyrRS,    
     &     num_spinstell_to_suppyrRS,   display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_suppyrFRB,             
     &     map_spinstell_to_suppyrFRB,   
     &     num_spinstell_to_suppyrFRB,  display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_supbask,               
     &     map_spinstell_to_supbask,     
     &     num_spinstell_to_supbask,    display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_supaxax,               
     &     map_spinstell_to_supaxax,     
     &     num_spinstell_to_supaxax,    display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_supLTS,                
     &     map_spinstell_to_supLTS,      
     &     num_spinstell_to_supLTS,     display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_spinstell,             
     &     map_spinstell_to_spinstell,   
     &     num_spinstell_to_spinstell,  display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_tuftIB,                
     &     map_spinstell_to_tuftIB,      
     &     num_spinstell_to_tuftIB,     display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_tuftRS,                
     &     map_spinstell_to_tuftRS,      
     &     num_spinstell_to_tuftRS,     display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_deepbask,              
     &     map_spinstell_to_deepbask,    
     &     num_spinstell_to_deepbask,   display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_deepng  ,              
     &     map_spinstell_to_deepng  ,    
     &     num_spinstell_to_deepng  ,   display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_deepaxax,              
     &     map_spinstell_to_deepaxax,    
     &     num_spinstell_to_deepaxax,   display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_deepLTS,               
     &     map_spinstell_to_deepLTS,     
     &     num_spinstell_to_deepLTS,    display)
          CALL synaptic_map_construct (thisno,
     &     num_spinstell,  num_nontuftRS,             
     &     map_spinstell_to_nontuftRS,   
     &     num_spinstell_to_nontuftRS,  display)

          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_suppyrRS,              
     &     map_tuftIB_to_suppyrRS,    
     &     num_tuftIB_to_suppyrRS,   display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_suppyrFRB,             
     &     map_tuftIB_to_suppyrFRB,   
     &     num_tuftIB_to_suppyrFRB,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_supbask,               
     &     map_tuftIB_to_supbask,     
     &     num_tuftIB_to_supbask,    display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_supaxax,               
     &     map_tuftIB_to_supaxax,     
     &     num_tuftIB_to_supaxax,    display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_supLTS,                
     &     map_tuftIB_to_supLTS,      
     &     num_tuftIB_to_supLTS,     display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_spinstell,             
     &     map_tuftIB_to_spinstell,   
     &     num_tuftIB_to_spinstell,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_tuftIB   ,             
     &     map_tuftIB_to_tuftIB   ,   
     &     num_tuftIB_to_tuftIB   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_tuftRS   ,             
     &     map_tuftIB_to_tuftRS   ,   
     &     num_tuftIB_to_tuftRS   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_deepbask ,             
     &     map_tuftIB_to_deepbask ,   
     &     num_tuftIB_to_deepbask ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_deepng   ,             
     &     map_tuftIB_to_deepng   ,   
     &     num_tuftIB_to_deepng   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_deepaxax ,             
     &     map_tuftIB_to_deepaxax ,   
     &     num_tuftIB_to_deepaxax ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_deepLTS  ,             
     &     map_tuftIB_to_deepLTS  ,   
     &     num_tuftIB_to_deepLTS  ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftIB,  num_nontuftRS,             
     &     map_tuftIB_to_nontuftRS,   
     &     num_tuftIB_to_nontuftRS,  display)

          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_suppyrRS ,             
     &     map_tuftRS_to_suppyrRS ,   
     &     num_tuftRS_to_suppyrRS ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_suppyrFRB,             
     &     map_tuftRS_to_suppyrFRB,   
     &     num_tuftRS_to_suppyrFRB,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_supbask  ,             
     &     map_tuftRS_to_supbask  ,   
     &     num_tuftRS_to_supbask  ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_supaxax  ,             
     &     map_tuftRS_to_supaxax  ,   
     &     num_tuftRS_to_supaxax  ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_supLTS   ,             
     &     map_tuftRS_to_supLTS   ,   
     &     num_tuftRS_to_supLTS   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_spinstell,             
     &     map_tuftRS_to_spinstell,   
     &     num_tuftRS_to_spinstell,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_tuftIB   ,             
     &     map_tuftRS_to_tuftIB   ,   
     &     num_tuftRS_to_tuftIB   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_tuftRS   ,             
     &     map_tuftRS_to_tuftRS   ,   
     &     num_tuftRS_to_tuftRS   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_deepbask ,             
     &     map_tuftRS_to_deepbask ,   
     &     num_tuftRS_to_deepbask ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_deepng   ,             
     &     map_tuftRS_to_deepng   ,   
     &     num_tuftRS_to_deepng   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_deepaxax ,             
     &     map_tuftRS_to_deepaxax ,   
     &     num_tuftRS_to_deepaxax ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_deepLTS  ,             
     &     map_tuftRS_to_deepLTS  ,   
     &     num_tuftRS_to_deepLTS  ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_tuftRS,  num_nontuftRS,             
     &     map_tuftRS_to_nontuftRS,   
     &     num_tuftRS_to_nontuftRS,  display)

          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_spinstell,             
     &     map_deepbask_to_spinstell,   
     &     num_deepbask_to_spinstell,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_tuftIB   ,             
     &     map_deepbask_to_tuftIB   ,   
     &     num_deepbask_to_tuftIB   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_tuftRS   ,             
     &     map_deepbask_to_tuftRS   ,   
     &     num_deepbask_to_tuftRS   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_deepbask ,             
     &     map_deepbask_to_deepbask ,   
     &     num_deepbask_to_deepbask ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_deepng   ,             
     &     map_deepbask_to_deepng   ,   
     &     num_deepbask_to_deepng   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_deepaxax ,             
     &     map_deepbask_to_deepaxax ,   
     &     num_deepbask_to_deepaxax ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_deepLTS  ,             
     &     map_deepbask_to_deepLTS  ,   
     &     num_deepbask_to_deepLTS  ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_nontuftRS,             
     &     map_deepbask_to_nontuftRS,   
     &     num_deepbask_to_nontuftRS,  display)

          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_tuftIB   ,             
     &     map_deepng_to_tuftIB   ,   
     &     num_deepng_to_tuftIB   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_tuftRS   ,             
     &     map_deepng_to_tuftRS   ,   
     &     num_deepng_to_tuftRS   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_nontuftRS,             
     &     map_deepng_to_nontuftRS,   
     &     num_deepng_to_nontuftRS,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_spinstell,             
     &     map_deepng_to_spinstell,   
     &     num_deepng_to_spinstell,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_deepng   ,             
     &     map_deepng_to_deepng   ,   
     &     num_deepng_to_deepng   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_deepbask ,             
     &     map_deepng_to_deepbask ,   
     &     num_deepng_to_deepbask ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_deepaxax,  num_suppyrRS ,             
     &     map_deepaxax_to_suppyrRS ,   
     &     num_deepaxax_to_suppyrRS ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepaxax,  num_suppyrFRB,             
     &     map_deepaxax_to_suppyrFRB,   
     &     num_deepaxax_to_suppyrFRB,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepaxax,  num_spinstell,             
     &     map_deepaxax_to_spinstell,   
     &     num_deepaxax_to_spinstell,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepaxax,  num_tuftIB   ,             
     &     map_deepaxax_to_tuftIB   ,   
     &     num_deepaxax_to_tuftIB   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepaxax,  num_tuftRS   ,             
     &     map_deepaxax_to_tuftRS   ,   
     &     num_deepaxax_to_tuftRS   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepaxax,  num_nontuftRS   ,             
     &     map_deepaxax_to_nontuftRS   ,   
     &     num_deepaxax_to_nontuftRS   ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_suppyrRS    ,             
     &     map_deepLTS_to_suppyrRS    ,   
     &     num_deepLTS_to_suppyrRS    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_suppyrFRB   ,             
     &     map_deepLTS_to_suppyrFRB   ,   
     &     num_deepLTS_to_suppyrFRB   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_supbask     ,             
     &     map_deepLTS_to_supbask     ,   
     &     num_deepLTS_to_supbask     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_supaxax     ,             
     &     map_deepLTS_to_supaxax     ,   
     &     num_deepLTS_to_supaxax     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_supLTS      ,             
     &     map_deepLTS_to_supLTS      ,   
     &     num_deepLTS_to_supLTS      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_spinstell   ,             
     &     map_deepLTS_to_spinstell   ,   
     &     num_deepLTS_to_spinstell   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_tuftIB      ,             
     &     map_deepLTS_to_tuftIB      ,   
     &     num_deepLTS_to_tuftIB      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_tuftRS      ,             
     &     map_deepLTS_to_tuftRS      ,   
     &     num_deepLTS_to_tuftRS      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_deepbask    ,             
     &     map_deepLTS_to_deepbask    ,   
     &     num_deepLTS_to_deepbask    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_deepaxax    ,             
     &     map_deepLTS_to_deepaxax    ,   
     &     num_deepLTS_to_deepaxax    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_deepLTS     ,             
     &     map_deepLTS_to_deepLTS     ,   
     &     num_deepLTS_to_deepLTS     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS ,  num_nontuftRS   ,             
     &     map_deepLTS_to_nontuftRS   ,   
     &     num_deepLTS_to_nontuftRS   ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_suppyrRS    ,             
     &     map_TCR_to_suppyrRS    ,   
     &     num_TCR_to_suppyrRS    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_suppyrFRB   ,             
     &     map_TCR_to_suppyrFRB   ,   
     &     num_TCR_to_suppyrFRB   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_supbask     ,             
     &     map_TCR_to_supbask     ,   
     &     num_TCR_to_supbask     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_supng       ,             
     &     map_TCR_to_supng       ,   
     &     num_TCR_to_supng       ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_supaxax     ,             
     &     map_TCR_to_supaxax     ,   
     &     num_TCR_to_supaxax     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_spinstell   ,             
     &     map_TCR_to_spinstell   ,   
     &     num_TCR_to_spinstell   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_tuftIB      ,             
     &     map_TCR_to_tuftIB      ,   
     &     num_TCR_to_tuftIB      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_tuftRS      ,             
     &     map_TCR_to_tuftRS      ,   
     &     num_TCR_to_tuftRS      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_deepbask    ,             
     &     map_TCR_to_deepbask    ,   
     &     num_TCR_to_deepbask    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_deepng      ,             
     &     map_TCR_to_deepng      ,   
     &     num_TCR_to_deepng      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_deepaxax    ,             
     &     map_TCR_to_deepaxax    ,   
     &     num_TCR_to_deepaxax    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_nRT         ,             
     &     map_TCR_to_nRT         ,   
     &     num_TCR_to_nRT         ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_TCR ,  num_nontuftRS   ,             
     &     map_TCR_to_nontuftRS   ,   
     &     num_TCR_to_nontuftRS   ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_nRT ,  num_TCR         ,             
     &     map_nRT_to_TCR         ,   
     &     num_nRT_to_TCR         ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nRT ,  num_nRT         ,             
     &     map_nRT_to_nRT         ,   
     &     num_nRT_to_nRT         ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_suppyrRS    ,             
     &     map_nontuftRS_to_suppyrRS    ,   
     &     num_nontuftRS_to_suppyrRS    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_suppyrFRB   ,             
     &     map_nontuftRS_to_suppyrFRB   ,   
     &     num_nontuftRS_to_suppyrFRB   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_supbask     ,             
     &     map_nontuftRS_to_supbask     ,   
     &     num_nontuftRS_to_supbask     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_supaxax     ,             
     &     map_nontuftRS_to_supaxax     ,   
     &     num_nontuftRS_to_supaxax     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_supLTS      ,             
     &     map_nontuftRS_to_supLTS      ,   
     &     num_nontuftRS_to_supLTS      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_spinstell   ,             
     &     map_nontuftRS_to_spinstell   ,   
     &     num_nontuftRS_to_spinstell   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_tuftIB      ,             
     &     map_nontuftRS_to_tuftIB      ,   
     &     num_nontuftRS_to_tuftIB      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_tuftRS      ,             
     &     map_nontuftRS_to_tuftRS      ,   
     &     num_nontuftRS_to_tuftRS      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_deepbask    ,             
     &     map_nontuftRS_to_deepbask    ,   
     &     num_nontuftRS_to_deepbask    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_deepng      ,             
     &     map_nontuftRS_to_deepng      ,   
     &     num_nontuftRS_to_deepng      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_deepaxax    ,             
     &     map_nontuftRS_to_deepaxax    ,   
     &     num_nontuftRS_to_deepaxax    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_deepLTS     ,             
     &     map_nontuftRS_to_deepLTS     ,   
     &     num_nontuftRS_to_deepLTS     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_TCR         ,             
     &     map_nontuftRS_to_TCR         ,   
     &     num_nontuftRS_to_TCR         ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_nRT         ,             
     &     map_nontuftRS_to_nRT         ,   
     &     num_nontuftRS_to_nRT         ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_nontuftRS ,  num_nontuftRS   ,             
     &     map_nontuftRS_to_nontuftRS   ,   
     &     num_nontuftRS_to_nontuftRS   ,  display)
c Finish construction of synaptic connection tables.

c Construct synaptic compartment maps.  
                display = 0

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS, com_suppyrRS_to_suppyrRS,           
     &     num_suppyrRS_to_suppyrRS,
     &  ncompallow_suppyrRS_to_suppyrRS,
     &   compallow_suppyrRS_to_suppyrRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_suppyrRS_to_suppyrFRB,          
     &     num_suppyrRS_to_suppyrFRB,
     &    ncompallow_suppyrRS_to_suppyrFRB,
     &     compallow_suppyrRS_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supbask  , com_suppyrRS_to_supbask,            
     &     num_suppyrRS_to_supbask,
     &    ncompallow_suppyrRS_to_supbask,  
     &     compallow_suppyrRS_to_supbask,   display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supng    , com_suppyrRS_to_supng  ,            
     &     num_suppyrRS_to_supng  ,
     &    ncompallow_suppyrRS_to_supng  ,  
     &     compallow_suppyrRS_to_supng  ,   display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supaxax  , com_suppyrRS_to_supaxax,            
     &     num_suppyrRS_to_supaxax,  
     &    ncompallow_suppyrRS_to_supaxax,  
     &     compallow_suppyrRS_to_supaxax,   display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supLTS   , com_suppyrRS_to_supLTS,             
     &     num_suppyrRS_to_supLTS,   
     &    ncompallow_suppyrRS_to_supLTS,   
     &     compallow_suppyrRS_to_supLTS,    display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_suppyrRS_to_spinstell,          
     &     num_suppyrRS_to_spinstell,
     &    ncompallow_suppyrRS_to_spinstell,
     &     compallow_suppyrRS_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_suppyrRS_to_tuftIB   ,          
     &     num_suppyrRS_to_tuftIB   ,
     &    ncompallow_suppyrRS_to_tuftIB   ,
     &     compallow_suppyrRS_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_suppyrRS_to_tuftRS   ,          
     &     num_suppyrRS_to_tuftRS   ,
     &    ncompallow_suppyrRS_to_tuftRS   ,
     &     compallow_suppyrRS_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_suppyrRS_to_deepbask ,          
     &     num_suppyrRS_to_deepbask ,
     &    ncompallow_suppyrRS_to_deepbask ,
     &     compallow_suppyrRS_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepaxax , com_suppyrRS_to_deepaxax ,          
     &     num_suppyrRS_to_deepaxax ,
     &    ncompallow_suppyrRS_to_deepaxax ,
     &     compallow_suppyrRS_to_deepaxax , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS  , com_suppyrRS_to_deepLTS  ,          
     &     num_suppyrRS_to_deepLTS  ,
     &    ncompallow_suppyrRS_to_deepLTS  ,
     &     compallow_suppyrRS_to_deepLTS  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_suppyrRS_to_nontuftRS,          
     &     num_suppyrRS_to_nontuftRS,
     &    ncompallow_suppyrRS_to_nontuftRS,
     &     compallow_suppyrRS_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_suppyrFRB_to_suppyrRS ,          
     &     num_suppyrFRB_to_suppyrRS ,
     &    ncompallow_suppyrFRB_to_suppyrRS ,
     &     compallow_suppyrFRB_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_suppyrFRB_to_suppyrFRB,          
     &     num_suppyrFRB_to_suppyrFRB,
     &    ncompallow_suppyrFRB_to_suppyrFRB,
     &     compallow_suppyrFRB_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supbask  , com_suppyrFRB_to_supbask  ,          
     &     num_suppyrFRB_to_supbask  ,
     &    ncompallow_suppyrFRB_to_supbask  ,
     &     compallow_suppyrFRB_to_supbask  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supng    , com_suppyrFRB_to_supng    ,          
     &     num_suppyrFRB_to_supng    ,
     &    ncompallow_suppyrFRB_to_supng    ,
     &     compallow_suppyrFRB_to_supng    , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supaxax  , com_suppyrFRB_to_supaxax  ,          
     &     num_suppyrFRB_to_supaxax  ,
     &    ncompallow_suppyrFRB_to_supaxax  ,
     &     compallow_suppyrFRB_to_supaxax  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supLTS   , com_suppyrFRB_to_supLTS   ,          
     &     num_suppyrFRB_to_supLTS   ,
     &    ncompallow_suppyrFRB_to_supLTS   ,
     &     compallow_suppyrFRB_to_supLTS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_suppyrFRB_to_spinstell,          
     &     num_suppyrFRB_to_spinstell,
     &    ncompallow_suppyrFRB_to_spinstell,
     &     compallow_suppyrFRB_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_suppyrFRB_to_tuftIB   ,          
     &     num_suppyrFRB_to_tuftIB   ,
     &    ncompallow_suppyrFRB_to_tuftIB   ,
     &     compallow_suppyrFRB_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_suppyrFRB_to_tuftRS   ,          
     &     num_suppyrFRB_to_tuftRS   ,
     &    ncompallow_suppyrFRB_to_tuftRS   ,
     &     compallow_suppyrFRB_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_suppyrFRB_to_deepbask ,          
     &     num_suppyrFRB_to_deepbask ,
     &    ncompallow_suppyrFRB_to_deepbask ,
     &     compallow_suppyrFRB_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepaxax , com_suppyrFRB_to_deepaxax ,          
     &     num_suppyrFRB_to_deepaxax ,
     &    ncompallow_suppyrFRB_to_deepaxax ,
     &     compallow_suppyrFRB_to_deepaxax , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS  , com_suppyrFRB_to_deepLTS  ,          
     &     num_suppyrFRB_to_deepLTS  ,
     &    ncompallow_suppyrFRB_to_deepLTS  ,
     &     compallow_suppyrFRB_to_deepLTS  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_suppyrFRB_to_nontuftRS,          
     &     num_suppyrFRB_to_nontuftRS,
     &    ncompallow_suppyrFRB_to_nontuftRS,
     &     compallow_suppyrFRB_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_supbask_to_suppyrRS ,          
     &     num_supbask_to_suppyrRS ,
     &    ncompallow_supbask_to_suppyrRS ,
     &     compallow_supbask_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_supbask_to_suppyrFRB,          
     &     num_supbask_to_suppyrFRB,
     &    ncompallow_supbask_to_suppyrFRB,
     &     compallow_supbask_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supbask  , com_supbask_to_supbask  ,          
     &     num_supbask_to_supbask  ,
     &    ncompallow_supbask_to_supbask  ,
     &     compallow_supbask_to_supbask  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supng    , com_supbask_to_supng    ,          
     &     num_supbask_to_supng    ,
     &    ncompallow_supbask_to_supng    ,
     &     compallow_supbask_to_supng    , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supaxax  , com_supbask_to_supaxax  ,          
     &     num_supbask_to_supaxax  ,
     &    ncompallow_supbask_to_supaxax  ,
     &     compallow_supbask_to_supaxax  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supLTS   , com_supbask_to_supLTS   ,          
     &     num_supbask_to_supLTS   ,
     &    ncompallow_supbask_to_supLTS   ,
     &     compallow_supbask_to_supLTS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_supbask_to_spinstell,          
     &     num_supbask_to_spinstell,
     &    ncompallow_supbask_to_spinstell,
     &     compallow_supbask_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_supng_to_suppyrRS ,          
     &     num_supng_to_suppyrRS ,
     &    ncompallow_supng_to_suppyrRS ,
     &     compallow_supng_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_supng_to_suppyrFRB,          
     &     num_supng_to_suppyrFRB,
     &    ncompallow_supng_to_suppyrFRB,
     &     compallow_supng_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_supng_to_nontuftRS,          
     &     num_supng_to_nontuftRS,
     &    ncompallow_supng_to_nontuftRS,
     &     compallow_supng_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_supng_to_tuftIB   ,          
     &     num_supng_to_tuftIB   ,
     &    ncompallow_supng_to_tuftIB   ,
     &     compallow_supng_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_supng_to_tuftRS   ,          
     &     num_supng_to_tuftRS   ,
     &    ncompallow_supng_to_tuftRS   ,
     &     compallow_supng_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supng    , com_supng_to_supng    ,          
     &     num_supng_to_supng    ,
     &    ncompallow_supng_to_supng    ,
     &     compallow_supng_to_supng    , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supbask  , com_supng_to_supbask  ,          
     &     num_supng_to_supbask  ,
     &    ncompallow_supng_to_supbask  ,
     &     compallow_supng_to_supbask  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_supaxax_to_suppyrRS ,          
     &     num_supaxax_to_suppyrRS ,
     &    ncompallow_supaxax_to_suppyrRS ,
     &     compallow_supaxax_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_supaxax_to_suppyrFRB,          
     &     num_supaxax_to_suppyrFRB,
     &    ncompallow_supaxax_to_suppyrFRB,
     &     compallow_supaxax_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_supaxax_to_spinstell,          
     &     num_supaxax_to_spinstell,
     &    ncompallow_supaxax_to_spinstell,
     &     compallow_supaxax_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_supaxax_to_tuftIB   ,          
     &     num_supaxax_to_tuftIB   ,
     &    ncompallow_supaxax_to_tuftIB   ,
     &     compallow_supaxax_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_supaxax_to_tuftRS   ,          
     &     num_supaxax_to_tuftRS   ,
     &    ncompallow_supaxax_to_tuftRS   ,
     &     compallow_supaxax_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_supaxax_to_nontuftRS,          
     &     num_supaxax_to_nontuftRS,
     &    ncompallow_supaxax_to_nontuftRS,
     &     compallow_supaxax_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_supLTS_to_suppyrRS ,          
     &     num_supLTS_to_suppyrRS ,
     &    ncompallow_supLTS_to_suppyrRS ,
     &     compallow_supLTS_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_supLTS_to_suppyrFRB,          
     &     num_supLTS_to_suppyrFRB,
     &    ncompallow_supLTS_to_suppyrFRB,
     &     compallow_supLTS_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supbask  , com_supLTS_to_supbask  ,          
     &     num_supLTS_to_supbask  ,
     &    ncompallow_supLTS_to_supbask  ,
     &     compallow_supLTS_to_supbask  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supaxax  , com_supLTS_to_supaxax  ,          
     &     num_supLTS_to_supaxax  ,
     &    ncompallow_supLTS_to_supaxax  ,
     &     compallow_supLTS_to_supaxax  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supLTS   , com_supLTS_to_supLTS   ,          
     &     num_supLTS_to_supLTS   ,
     &    ncompallow_supLTS_to_supLTS   ,
     &     compallow_supLTS_to_supLTS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_supLTS_to_spinstell,          
     &     num_supLTS_to_spinstell,
     &    ncompallow_supLTS_to_spinstell,
     &     compallow_supLTS_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_supLTS_to_tuftIB   ,          
     &     num_supLTS_to_tuftIB   ,
     &    ncompallow_supLTS_to_tuftIB   ,
     &     compallow_supLTS_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_supLTS_to_tuftRS   ,          
     &     num_supLTS_to_tuftRS   ,
     &    ncompallow_supLTS_to_tuftRS   ,
     &     compallow_supLTS_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_supLTS_to_deepbask ,          
     &     num_supLTS_to_deepbask ,
     &    ncompallow_supLTS_to_deepbask ,
     &     compallow_supLTS_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepaxax , com_supLTS_to_deepaxax ,          
     &     num_supLTS_to_deepaxax ,
     &    ncompallow_supLTS_to_deepaxax ,
     &     compallow_supLTS_to_deepaxax , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS  , com_supLTS_to_deepLTS  ,          
     &     num_supLTS_to_deepLTS  ,
     &    ncompallow_supLTS_to_deepLTS  ,
     &     compallow_supLTS_to_deepLTS  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_supLTS_to_nontuftRS,          
     &     num_supLTS_to_nontuftRS,
     &    ncompallow_supLTS_to_nontuftRS,
     &     compallow_supLTS_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_spinstell_to_suppyrRS ,          
     &     num_spinstell_to_suppyrRS ,
     &    ncompallow_spinstell_to_suppyrRS ,
     &     compallow_spinstell_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_spinstell_to_suppyrFRB,          
     &     num_spinstell_to_suppyrFRB,
     &    ncompallow_spinstell_to_suppyrFRB,
     &     compallow_spinstell_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supbask  , com_spinstell_to_supbask  ,          
     &     num_spinstell_to_supbask  ,
     &    ncompallow_spinstell_to_supbask  ,
     &     compallow_spinstell_to_supbask  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supaxax  , com_spinstell_to_supaxax  ,          
     &     num_spinstell_to_supaxax  ,
     &    ncompallow_spinstell_to_supaxax  ,
     &     compallow_spinstell_to_supaxax  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supLTS   , com_spinstell_to_supLTS   ,          
     &     num_spinstell_to_supLTS   ,
     &    ncompallow_spinstell_to_supLTS   ,
     &     compallow_spinstell_to_supLTS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_spinstell_to_spinstell,          
     &     num_spinstell_to_spinstell,
     &    ncompallow_spinstell_to_spinstell,
     &     compallow_spinstell_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_spinstell_to_tuftIB   ,          
     &     num_spinstell_to_tuftIB   ,
     &    ncompallow_spinstell_to_tuftIB   ,
     &     compallow_spinstell_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_spinstell_to_tuftRS   ,          
     &     num_spinstell_to_tuftRS   ,
     &    ncompallow_spinstell_to_tuftRS   ,
     &     compallow_spinstell_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_spinstell_to_deepbask ,          
     &     num_spinstell_to_deepbask ,
     &    ncompallow_spinstell_to_deepbask ,
     &     compallow_spinstell_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_spinstell_to_deepng   ,          
     &     num_spinstell_to_deepng   ,
     &    ncompallow_spinstell_to_deepng   ,
     &     compallow_spinstell_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepaxax , com_spinstell_to_deepaxax ,          
     &     num_spinstell_to_deepaxax ,
     &    ncompallow_spinstell_to_deepaxax ,
     &     compallow_spinstell_to_deepaxax , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS  , com_spinstell_to_deepLTS  ,          
     &     num_spinstell_to_deepLTS  ,
     &    ncompallow_spinstell_to_deepLTS  ,
     &     compallow_spinstell_to_deepLTS  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_spinstell_to_nontuftRS,          
     &     num_spinstell_to_nontuftRS,
     &    ncompallow_spinstell_to_nontuftRS,
     &     compallow_spinstell_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_tuftIB_to_suppyrRS ,          
     &     num_tuftIB_to_suppyrRS ,
     &    ncompallow_tuftIB_to_suppyrRS ,
     &     compallow_tuftIB_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_tuftIB_to_suppyrFRB,          
     &     num_tuftIB_to_suppyrFRB,
     &    ncompallow_tuftIB_to_suppyrFRB,
     &     compallow_tuftIB_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supbask  , com_tuftIB_to_supbask  ,          
     &     num_tuftIB_to_supbask  ,
     &    ncompallow_tuftIB_to_supbask  ,
     &     compallow_tuftIB_to_supbask  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supaxax  , com_tuftIB_to_supaxax  ,          
     &     num_tuftIB_to_supaxax  ,
     &    ncompallow_tuftIB_to_supaxax  ,
     &     compallow_tuftIB_to_supaxax  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supLTS   , com_tuftIB_to_supLTS   ,          
     &     num_tuftIB_to_supLTS   ,
     &    ncompallow_tuftIB_to_supLTS   ,
     &     compallow_tuftIB_to_supLTS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_tuftIB_to_spinstell,          
     &     num_tuftIB_to_spinstell,
     &    ncompallow_tuftIB_to_spinstell,
     &     compallow_tuftIB_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_tuftIB_to_tuftIB   ,          
     &     num_tuftIB_to_tuftIB   ,
     &    ncompallow_tuftIB_to_tuftIB   ,
     &     compallow_tuftIB_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_tuftIB_to_tuftRS   ,          
     &     num_tuftIB_to_tuftRS   ,
     &    ncompallow_tuftIB_to_tuftRS   ,
     &     compallow_tuftIB_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_tuftIB_to_deepbask ,          
     &     num_tuftIB_to_deepbask ,
     &    ncompallow_tuftIB_to_deepbask ,
     &     compallow_tuftIB_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_tuftIB_to_deepng   ,          
     &     num_tuftIB_to_deepng   ,
     &    ncompallow_tuftIB_to_deepng   ,
     &     compallow_tuftIB_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepaxax , com_tuftIB_to_deepaxax ,          
     &     num_tuftIB_to_deepaxax ,
     &    ncompallow_tuftIB_to_deepaxax ,
     &     compallow_tuftIB_to_deepaxax , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS  , com_tuftIB_to_deepLTS  ,          
     &     num_tuftIB_to_deepLTS  ,
     &    ncompallow_tuftIB_to_deepLTS  ,
     &     compallow_tuftIB_to_deepLTS  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_tuftIB_to_nontuftRS,          
     &     num_tuftIB_to_nontuftRS,
     &    ncompallow_tuftIB_to_nontuftRS,
     &     compallow_tuftIB_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_tuftRS_to_suppyrRS ,          
     &     num_tuftRS_to_suppyrRS ,
     &    ncompallow_tuftRS_to_suppyrRS ,
     &     compallow_tuftRS_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_tuftRS_to_suppyrFRB,          
     &     num_tuftRS_to_suppyrFRB,
     &    ncompallow_tuftRS_to_suppyrFRB,
     &     compallow_tuftRS_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supbask  , com_tuftRS_to_supbask  ,          
     &     num_tuftRS_to_supbask  ,
     &    ncompallow_tuftRS_to_supbask  ,
     &     compallow_tuftRS_to_supbask  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supaxax  , com_tuftRS_to_supaxax  ,          
     &     num_tuftRS_to_supaxax  ,
     &    ncompallow_tuftRS_to_supaxax  ,
     &     compallow_tuftRS_to_supaxax  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supLTS   , com_tuftRS_to_supLTS   ,          
     &     num_tuftRS_to_supLTS   ,
     &    ncompallow_tuftRS_to_supLTS   ,
     &     compallow_tuftRS_to_supLTS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_tuftRS_to_spinstell,          
     &     num_tuftRS_to_spinstell,
     &    ncompallow_tuftRS_to_spinstell,
     &     compallow_tuftRS_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_tuftRS_to_tuftIB   ,          
     &     num_tuftRS_to_tuftIB   ,
     &    ncompallow_tuftRS_to_tuftIB   ,
     &     compallow_tuftRS_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_tuftRS_to_tuftRS   ,          
     &     num_tuftRS_to_tuftRS   ,
     &    ncompallow_tuftRS_to_tuftRS   ,
     &     compallow_tuftRS_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_tuftRS_to_deepbask ,          
     &     num_tuftRS_to_deepbask ,
     &    ncompallow_tuftRS_to_deepbask ,
     &     compallow_tuftRS_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_tuftRS_to_deepng   ,          
     &     num_tuftRS_to_deepng   ,
     &    ncompallow_tuftRS_to_deepng   ,
     &     compallow_tuftRS_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepaxax , com_tuftRS_to_deepaxax ,          
     &     num_tuftRS_to_deepaxax ,
     &    ncompallow_tuftRS_to_deepaxax ,
     &     compallow_tuftRS_to_deepaxax , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS  , com_tuftRS_to_deepLTS  ,          
     &     num_tuftRS_to_deepLTS  ,
     &    ncompallow_tuftRS_to_deepLTS  ,
     &     compallow_tuftRS_to_deepLTS  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_tuftRS_to_nontuftRS,          
     &     num_tuftRS_to_nontuftRS,
     &    ncompallow_tuftRS_to_nontuftRS,
     &     compallow_tuftRS_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_deepbask_to_spinstell,          
     &     num_deepbask_to_spinstell,
     &    ncompallow_deepbask_to_spinstell,
     &     compallow_deepbask_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_deepbask_to_tuftIB   ,          
     &     num_deepbask_to_tuftIB   ,
     &    ncompallow_deepbask_to_tuftIB   ,
     &     compallow_deepbask_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_deepbask_to_tuftRS   ,          
     &     num_deepbask_to_tuftRS   ,
     &    ncompallow_deepbask_to_tuftRS   ,
     &     compallow_deepbask_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_deepbask_to_deepbask ,          
     &     num_deepbask_to_deepbask ,
     &    ncompallow_deepbask_to_deepbask ,
     &     compallow_deepbask_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_deepbask_to_deepng   ,          
     &     num_deepbask_to_deepng   ,
     &    ncompallow_deepbask_to_deepng   ,
     &     compallow_deepbask_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepaxax , com_deepbask_to_deepaxax ,          
     &     num_deepbask_to_deepaxax ,
     &    ncompallow_deepbask_to_deepaxax ,
     &     compallow_deepbask_to_deepaxax , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS  , com_deepbask_to_deepLTS  ,          
     &     num_deepbask_to_deepLTS  ,
     &    ncompallow_deepbask_to_deepLTS  ,
     &     compallow_deepbask_to_deepLTS  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_deepbask_to_nontuftRS,          
     &     num_deepbask_to_nontuftRS,
     &    ncompallow_deepbask_to_nontuftRS,
     &     compallow_deepbask_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_deepng_to_tuftIB   ,          
     &     num_deepng_to_tuftIB   ,
     &    ncompallow_deepng_to_tuftIB   ,
     &     compallow_deepng_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_deepng_to_tuftRS   ,          
     &     num_deepng_to_tuftRS   ,
     &    ncompallow_deepng_to_tuftRS   ,
     &     compallow_deepng_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_deepng_to_nontuftRS,          
     &     num_deepng_to_nontuftRS,
     &    ncompallow_deepng_to_nontuftRS,
     &     compallow_deepng_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_deepng_to_spinstell,          
     &     num_deepng_to_spinstell,
     &    ncompallow_deepng_to_spinstell,
     &     compallow_deepng_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_deepng_to_deepng   ,          
     &     num_deepng_to_deepng   ,
     &    ncompallow_deepng_to_deepng   ,
     &     compallow_deepng_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_deepng_to_deepbask ,          
     &     num_deepng_to_deepbask ,
     &    ncompallow_deepng_to_deepbask ,
     &     compallow_deepng_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_deepaxax_to_suppyrRS ,          
     &     num_deepaxax_to_suppyrRS ,
     &    ncompallow_deepaxax_to_suppyrRS ,
     &     compallow_deepaxax_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_deepaxax_to_suppyrFRB,          
     &     num_deepaxax_to_suppyrFRB,
     &    ncompallow_deepaxax_to_suppyrFRB,
     &     compallow_deepaxax_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_deepaxax_to_spinstell,          
     &     num_deepaxax_to_spinstell,
     &    ncompallow_deepaxax_to_spinstell,
     &     compallow_deepaxax_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_deepaxax_to_tuftIB   ,          
     &     num_deepaxax_to_tuftIB   ,
     &    ncompallow_deepaxax_to_tuftIB   ,
     &     compallow_deepaxax_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_deepaxax_to_tuftRS   ,          
     &     num_deepaxax_to_tuftRS   ,
     &    ncompallow_deepaxax_to_tuftRS   ,
     &     compallow_deepaxax_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_deepaxax_to_nontuftRS,          
     &     num_deepaxax_to_nontuftRS,
     &    ncompallow_deepaxax_to_nontuftRS,
     &     compallow_deepaxax_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_deepLTS_to_suppyrRS ,          
     &     num_deepLTS_to_suppyrRS ,
     &    ncompallow_deepLTS_to_suppyrRS ,
     &     compallow_deepLTS_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_deepLTS_to_suppyrFRB,          
     &     num_deepLTS_to_suppyrFRB,
     &    ncompallow_deepLTS_to_suppyrFRB,
     &     compallow_deepLTS_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supbask  , com_deepLTS_to_supbask  ,          
     &     num_deepLTS_to_supbask  ,
     &    ncompallow_deepLTS_to_supbask  ,
     &     compallow_deepLTS_to_supbask  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supaxax  , com_deepLTS_to_supaxax  ,          
     &     num_deepLTS_to_supaxax  ,
     &    ncompallow_deepLTS_to_supaxax  ,
     &     compallow_deepLTS_to_supaxax  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supLTS   , com_deepLTS_to_supLTS   ,          
     &     num_deepLTS_to_supLTS   ,
     &    ncompallow_deepLTS_to_supLTS   ,
     &     compallow_deepLTS_to_supLTS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_deepLTS_to_spinstell,          
     &     num_deepLTS_to_spinstell,
     &    ncompallow_deepLTS_to_spinstell,
     &     compallow_deepLTS_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_deepLTS_to_tuftIB   ,          
     &     num_deepLTS_to_tuftIB   ,
     &    ncompallow_deepLTS_to_tuftIB   ,
     &     compallow_deepLTS_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_deepLTS_to_tuftRS   ,          
     &     num_deepLTS_to_tuftRS   ,
     &    ncompallow_deepLTS_to_tuftRS   ,
     &     compallow_deepLTS_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_deepLTS_to_deepbask ,          
     &     num_deepLTS_to_deepbask ,
     &    ncompallow_deepLTS_to_deepbask ,
     &     compallow_deepLTS_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepaxax , com_deepLTS_to_deepaxax ,          
     &     num_deepLTS_to_deepaxax ,
     &    ncompallow_deepLTS_to_deepaxax ,
     &     compallow_deepLTS_to_deepaxax , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS  , com_deepLTS_to_deepLTS  ,          
     &     num_deepLTS_to_deepLTS  ,
     &    ncompallow_deepLTS_to_deepLTS  ,
     &     compallow_deepLTS_to_deepLTS  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_deepLTS_to_nontuftRS,          
     &     num_deepLTS_to_nontuftRS,
     &    ncompallow_deepLTS_to_nontuftRS,
     &     compallow_deepLTS_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_TCR_to_suppyrRS ,          
     &     num_TCR_to_suppyrRS ,
     &    ncompallow_TCR_to_suppyrRS ,
     &     compallow_TCR_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_TCR_to_suppyrFRB,          
     &     num_TCR_to_suppyrFRB,
     &    ncompallow_TCR_to_suppyrFRB,
     &     compallow_TCR_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supbask  , com_TCR_to_supbask  ,          
     &     num_TCR_to_supbask  ,
     &    ncompallow_TCR_to_supbask  ,
     &     compallow_TCR_to_supbask  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supng    , com_TCR_to_supng    ,          
     &     num_TCR_to_supng    ,
     &    ncompallow_TCR_to_supng    ,
     &     compallow_TCR_to_supng    , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supaxax  , com_TCR_to_supaxax  ,          
     &     num_TCR_to_supaxax  ,
     &    ncompallow_TCR_to_supaxax  ,
     &     compallow_TCR_to_supaxax  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_TCR_to_spinstell,          
     &     num_TCR_to_spinstell,
     &    ncompallow_TCR_to_spinstell,
     &     compallow_TCR_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_TCR_to_tuftIB   ,          
     &     num_TCR_to_tuftIB   ,
     &    ncompallow_TCR_to_tuftIB   ,
     &     compallow_TCR_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_TCR_to_tuftRS   ,          
     &     num_TCR_to_tuftRS   ,
     &    ncompallow_TCR_to_tuftRS   ,
     &     compallow_TCR_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_TCR_to_deepbask ,          
     &     num_TCR_to_deepbask ,
     &    ncompallow_TCR_to_deepbask ,
     &     compallow_TCR_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_TCR_to_deepng   ,          
     &     num_TCR_to_deepng   ,
     &    ncompallow_TCR_to_deepng   ,
     &     compallow_TCR_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepaxax , com_TCR_to_deepaxax ,          
     &     num_TCR_to_deepaxax ,
     &    ncompallow_TCR_to_deepaxax ,
     &     compallow_TCR_to_deepaxax , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nRT      , com_TCR_to_nRT      ,          
     &     num_TCR_to_nRT      ,
     &    ncompallow_TCR_to_nRT      ,
     &     compallow_TCR_to_nRT      , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_TCR_to_nontuftRS,          
     &     num_TCR_to_nontuftRS,
     &    ncompallow_TCR_to_nontuftRS,
     &     compallow_TCR_to_nontuftRS, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_TCR      , com_nRT_to_TCR      ,          
     &     num_nRT_to_TCR      ,
     &    ncompallow_nRT_to_TCR      ,
     &     compallow_nRT_to_TCR      , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nRT      , com_nRT_to_nRT      ,          
     &     num_nRT_to_nRT      ,
     &    ncompallow_nRT_to_nRT      ,
     &     compallow_nRT_to_nRT      , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrRS , com_nontuftRS_to_suppyrRS ,          
     &     num_nontuftRS_to_suppyrRS ,
     &    ncompallow_nontuftRS_to_suppyrRS ,
     &     compallow_nontuftRS_to_suppyrRS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_suppyrFRB, com_nontuftRS_to_suppyrFRB,          
     &     num_nontuftRS_to_suppyrFRB,
     &    ncompallow_nontuftRS_to_suppyrFRB,
     &     compallow_nontuftRS_to_suppyrFRB, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supbask  , com_nontuftRS_to_supbask  ,          
     &     num_nontuftRS_to_supbask  ,
     &    ncompallow_nontuftRS_to_supbask  ,
     &     compallow_nontuftRS_to_supbask  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supaxax  , com_nontuftRS_to_supaxax  ,          
     &     num_nontuftRS_to_supaxax  ,
     &    ncompallow_nontuftRS_to_supaxax  ,
     &     compallow_nontuftRS_to_supaxax  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supLTS   , com_nontuftRS_to_supLTS   ,          
     &     num_nontuftRS_to_supLTS   ,
     &    ncompallow_nontuftRS_to_supLTS   ,
     &     compallow_nontuftRS_to_supLTS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_spinstell, com_nontuftRS_to_spinstell,          
     &     num_nontuftRS_to_spinstell,
     &    ncompallow_nontuftRS_to_spinstell,
     &     compallow_nontuftRS_to_spinstell, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftIB   , com_nontuftRS_to_tuftIB   ,          
     &     num_nontuftRS_to_tuftIB   ,
     &    ncompallow_nontuftRS_to_tuftIB   ,
     &     compallow_nontuftRS_to_tuftIB   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_tuftRS   , com_nontuftRS_to_tuftRS   ,          
     &     num_nontuftRS_to_tuftRS   ,
     &    ncompallow_nontuftRS_to_tuftRS   ,
     &     compallow_nontuftRS_to_tuftRS   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_nontuftRS_to_deepbask ,          
     &     num_nontuftRS_to_deepbask ,
     &    ncompallow_nontuftRS_to_deepbask ,
     &     compallow_nontuftRS_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_nontuftRS_to_deepng   ,          
     &     num_nontuftRS_to_deepng   ,
     &    ncompallow_nontuftRS_to_deepng   ,
     &     compallow_nontuftRS_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepaxax , com_nontuftRS_to_deepaxax ,          
     &     num_nontuftRS_to_deepaxax ,
     &    ncompallow_nontuftRS_to_deepaxax ,
     &     compallow_nontuftRS_to_deepaxax , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS  , com_nontuftRS_to_deepLTS  ,          
     &     num_nontuftRS_to_deepLTS  ,
     &    ncompallow_nontuftRS_to_deepLTS  ,
     &     compallow_nontuftRS_to_deepLTS  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_TCR      , com_nontuftRS_to_TCR      ,          
     &     num_nontuftRS_to_TCR      ,
     &    ncompallow_nontuftRS_to_TCR      ,
     &     compallow_nontuftRS_to_TCR      , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nRT      , com_nontuftRS_to_nRT      ,          
     &     num_nontuftRS_to_nRT      ,
     &    ncompallow_nontuftRS_to_nRT      ,
     &     compallow_nontuftRS_to_nRT      , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_nontuftRS, com_nontuftRS_to_nontuftRS,          
     &     num_nontuftRS_to_nontuftRS,
     &    ncompallow_nontuftRS_to_nontuftRS,
     &     compallow_nontuftRS_to_nontuftRS, display)

c Finish construction of synaptic compartment maps. 

c Display TCR-to-spinstell connectivity: input cell, compartment
c        write(6,988)
988      format('  TCR to spinstell connectivity')
         do i = 1, num_spinstell
c          write(6,989) i
989        format(' spinstell number ',i4)
           do j = 1, num_TCR_to_spinstell
c           write(6,990) j, map_TCR_to_spinstell(j,i),
c    &             com_TCR_to_spinstell(j,i)
990         format(i5,i5,i5)
           end do
          end do


c Construct gap-junction tables
! axax interneurons a special case
           gjtable_supaxax(1,1) = 1
           gjtable_supaxax(1,2) = 12
           gjtable_supaxax(1,3) = 2
           gjtable_supaxax(1,4) = 12
           gjtable_deepaxax(1,1) = 1
           gjtable_deepaxax(1,2) = 12
           gjtable_deepaxax(1,3) = 2
           gjtable_deepaxax(1,4) = 12

      CALL groucho_gapbld (thisno, num_suppyrRS,
     & totaxgj_suppyrRS  , gjtable_suppyrRS,
     & table_axgjcompallow_suppyrRS, 
     & num_axgjcompallow_suppyrRS, 0) 

      CALL groucho_gapbld (thisno, num_suppyrFRB,
     & totaxgj_suppyrFRB , gjtable_suppyrFRB,
     & table_axgjcompallow_suppyrFRB,
     & num_axgjcompallow_suppyrFRB,  0) 

      call GROUCHO_gapbld_mix (thisno, num_suppyrRS,        
     & num_suppyrFRB, totaxgj_suppyr, gjtable_suppyr,
     & table_axgjcompallow_suppyrRS,
     & num_axgjcompallow_suppyrRS, 0)
      
      CALL groucho_gapbld (thisno, num_spinstell,
     & totaxgj_spinstell , gjtable_spinstell,
     & table_axgjcompallow_spinstell,
     & num_axgjcompallow_spinstell,  0) 

      CALL groucho_gapbld (thisno, num_tuftIB,   
     & totaxgj_tuftIB    , gjtable_tuftIB   ,
     & table_axgjcompallow_tuftIB   ,
     & num_axgjcompallow_tuftIB   ,  0) 

      CALL groucho_gapbld (thisno, num_tuftRS,   
     & totaxgj_tuftRS    , gjtable_tuftRS   ,
     & table_axgjcompallow_tuftRS   ,
     & num_axgjcompallow_tuftRS   ,  0) 

      call GROUCHO_gapbld_mix (thisno, num_tuftIB,        
     & num_tuftRS, totaxgj_tuft, gjtable_tuft,
     & table_axgjcompallow_tuftIB,
     & num_axgjcompallow_tuftIB, 0)
      
      CALL groucho_gapbld (thisno, num_nontuftRS,   
     & totaxgj_nontuftRS    , gjtable_nontuftRS   ,
     & table_axgjcompallow_nontuftRS   ,
     & num_axgjcompallow_nontuftRS   ,  0) 

      CALL groucho_gapbld (thisno, num_supbask  ,   
     & totSDgj_supbask      , gjtable_supbask     ,
     & table_SDgjcompallow_supbask     ,
     & num_SDgjcompallow_supbask     ,  0) 

      CALL groucho_gapbld (thisno, num_supng    ,   
     & totSDgj_supng        , gjtable_supng       ,
     & table_SDgjcompallow_supng       ,
     & num_SDgjcompallow_supng       ,  0) 

      CALL groucho_gapbld (thisno, num_supLTS   ,   
     & totSDgj_supLTS       , gjtable_supLTS      ,
     & table_SDgjcompallow_supLTS      ,
     & num_SDgjcompallow_supLTS      ,  0) 

      CALL groucho_gapbld (thisno, num_deepbask ,   
     & totSDgj_deepbask     , gjtable_deepbask    ,
     & table_SDgjcompallow_deepbask    ,
     & num_SDgjcompallow_deepbask    ,  0) 

      CALL groucho_gapbld (thisno, num_deepng   ,   
     & totSDgj_deepng       , gjtable_deepng      ,
     & table_SDgjcompallow_deepng      ,
     & num_SDgjcompallow_deepng      ,  0) 

      CALL groucho_gapbld (thisno, num_deepLTS  ,   
     & totSDgj_deepLTS      , gjtable_deepLTS     ,
     & table_SDgjcompallow_deepLTS     ,
     & num_SDgjcompallow_deepLTS     ,  0) 

      CALL groucho_gapbld (thisno, num_TCR      ,   
     & totaxgj_TCR          , gjtable_TCR         ,
     & table_axgjcompallow_TCR         ,
     & num_axgjcompallow_TCR         ,  0) 

      CALL groucho_gapbld (thisno, num_nRT      ,   
     & totSDgj_nRT          , gjtable_nRT         ,
     & table_SDgjcompallow_nRT         ,
     & num_SDgjcompallow_nRT         ,  0) 

! Define spread of values for gGABA_nRT_to_TCR
       call durand(seed,num_nRT,ranvec_nRT)
       do L = 1, num_nRT
        gGABA_nRT_to_TCR(L) = 0.7d-3 + 1.4d-3 * ranvec_nRT(L)
       end do

! Define tonic currents to different cell types
       call durand(seed,num_suppyrRS ,ranvec_suppyrRS )
       do L = 1, num_suppyrRS 
       do i = 69, 74  ! axonal compartments
c       curr_suppyrRS  (i,L) =  0.030d0
        curr_suppyrRS  (i,L) =  0.045d0
       end do
c       curr_suppyrRS  (1,L) =  0.00d0 + 0.01d0 *
        curr_suppyrRS  (1,L) =  0.35d0 + 0.01d0 *
     &    ranvec_suppyrRS (L)
       do i = 41, 42  ! apical dendrite
         curr_suppyrRS  (i,L) = -0.25d0
       end do
       end do
c Possible bias to apical dendrites of one cell
       curr_suppyrRS (41, 2) = -0.15d0
       curr_suppyrRS (42, 2) = -0.15d0
c Possible bias to apical dendrites of all cells
       do L = 1, num_suppyrRS
         curr_suppyrRS (41, L) = -0.15d0
         curr_suppyrRS (42, L) = -0.15d0
       end do

       call durand(seed,num_suppyrFRB,ranvec_suppyrFRB)
       do L = 1, num_suppyrFRB
       do i = 69, 74  ! axonal compartments
        curr_suppyrFRB (i,L) = 0.250d0 + 0.005d0 *
c       curr_suppyrFRB (1,L) = -0.500d0 + 0.05d0 *
     &    ranvec_suppyrFRB(L)
       end do
       end do

       call durand(seed,num_supbask  ,ranvec_supbask  )
       do L = 1, num_supbask    
c       curr_supbask   (1,L) = -0.10d0 + 0.02d0 *
        curr_supbask   (1,L) = -0.06d0 + 0.01d0 *
     &    ranvec_supbask  (L)
       end do

       call durand(seed,num_supaxax  ,ranvec_supaxax  )
       do L = 1, num_supaxax    
c       curr_supaxax   (1,L) = -0.10d0 + 0.02d0 *
        curr_supaxax   (1,L) = -0.06d0 + 0.01d0 *
     &    ranvec_supaxax  (L)
       end do

       call durand(seed,num_deepbask  ,ranvec_deepbask  )
       do L = 1, num_deepbask    
        curr_deepbask   (1,L) = -0.10d0 + 0.02d0 *
     &    ranvec_deepbask  (L)
       end do


       call durand(seed,num_spinstell,ranvec_spinstell)
       do L = 1, num_spinstell  
c       curr_spinstell (1,L) = -0.10d0 + 0.05d0 *
c       curr_spinstell (1,L) = -0.25d0 + 0.05d0 *
c       curr_spinstell (1,L) =  0.05d0 + 0.03d0 *
        curr_spinstell (1,L) =  0.07d0 + 0.03d0 *
     &    ranvec_spinstell(L)
       end do

       call durand(seed,num_deepLTS   ,ranvec_deepLTS   )
       do L = 1, num_deepLTS
         curr_deepLTS (1,L) = 0.28d0 + 0.01d0 *
     &      ranvec_deepLTS (L)
       end do

       call durand(seed,num_tuftIB   ,ranvec_tuftIB   )
       do L = 1, num_tuftIB    
!       curr_tuftIB    (1,L) = 0.10d0 + 0.1d0 *
!       curr_tuftIB    (1,L) = 0.00d0 + 0.1d0 *
!       curr_tuftIB    (1,L) = -1.00d0 + 0.1d0 * ! ? suppress intrinsic bursting
!       curr_tuftIB    (1,L) = 1.00d0 + 0.1d0 *  !  ? wake up cortex ?
!       curr_tuftIB    (1,L) = 1.50d0 + 0.1d0 *  !  ? REALLY wake up cortex ?
          do i = 2, 34
c       curr_tuftIB    (i,L) = 0.00d0 + 0.040d0 *  ! current to basal/oblique dendrites 
c       curr_tuftIB    (i,L) =-0.02d0 + 0.040d0 *  ! current to basal/oblique dendrites 
! -0.02 + 0.04...  used in interactB5
c       curr_tuftIB    (i,L) =-0.01d0 + 0.01d0 *  ! current to basal/oblique dendrites 
        curr_tuftIB    (i,L) = 0.020d0 + 0.01d0 *  ! current to basal/oblique dendrites 
     &    ranvec_tuftIB   (L)
          end do
        curr_tuftIB (57,L) = -0.015d0 ! axon
        curr_tuftIB (58,L) = -0.015d0 ! axon
        curr_tuftIB (59,L) = -0.015d0 ! axon
        curr_tuftIB (60,L) = -0.015d0 ! axon
        curr_tuftIB (61,L) = -0.015d0 ! axon
       end do

       call durand(seed,num_tuftRS   ,ranvec_tuftRS   )
       do L = 1, num_tuftRS    
!       curr_tuftRS    (1,L) = 0.10d0 + 0.1d0 *
!       curr_tuftRS    (1,L) = 0.00d0 + 0.1d0 *
!       curr_tuftRS    (1,L) = 1.00d0 + 0.1d0 *  !  ? wake up cortex?
          do i = 2, 34
        curr_tuftRS    (i,L) = -0.05d0 + 0.01d0 *  ! current to basal/oblique dendrites 
     &    ranvec_tuftRS   (L)
          end do
        curr_tuftRS (57,L) = -0.04d0 ! axon
        curr_tuftRS (58,L) = -0.04d0 ! axon
        curr_tuftRS (59,L) = -0.04d0 ! axon
        curr_tuftRS (60,L) = -0.04d0 ! axon
        curr_tuftRS (61,L) = -0.04d0 ! axon
       end do

       do L = 1, num_supng
          curr_supng (1,L) = -0.03d0 ! to suppress spontaneous firing
       end do

       do L = 1, num_deepng
c         curr_deepng (1,L) = -0.03d0 ! to suppress spontaneous firing
c         curr_deepng (1,L) = -0.025d0 ! to suppress spontaneous firing
          curr_deepng (1,L) = -0.045d0 ! to suppress spontaneous firing
c         curr_deepng (1,L) = -0.06d0 ! increase the hyperpol. curr for delta90
       end do

       call durand(seed,num_nontuftRS ,ranvec_nontuftRS )
       do L = 1, num_nontuftRS    
          if (L.le.(num_nontuftRS - nontuftRS_nFRB)) then
             z = 0.60d0
!            z = 0.50d0
          else
c            z = 0.0d0 ! 27 Jan 2005
!            z = 0.4d0
             z = 1.0d0
          endif
c       curr_nontuftRS  (1,L) = z + 0.1d0 *
c       curr_nontuftRS  (1,L) = z - 0.4d0 * ! decrease spont. firing
        curr_nontuftRS  (1,L) = z - 0.2d0 * ! decrease spont. firing
!       curr_nontuftRS  (1,L) = 0.75d0 + 0.1d0 *  ! ? wake up cortex
     &    ranvec_nontuftRS (L)
       end do

       call durand(seed,num_nRT       ,ranvec_nRT       )
       do L = 1, num_nRT          
        curr_nRT        (1,L) = -0.05d0 + 0.05d0 *
     &    ranvec_nRT       (L)
c        if (thisno.eq.0) then
c        write(6,8782) L, curr_nRT(1,L)
c        endif
8782     format(i4,3x,f8.3)
       end do

! During sz, curr to TCR can be zero
       call durand(seed,num_TCR       ,ranvec_TCR       )
       do L = 1, num_TCR          
!       curr_TCR        (1,L) = 1.40d0 + 0.01d0 *
        curr_TCR        (1,L) = 0.00d0 + 0.01d0 *
     &    ranvec_TCR       (L)
       end do

! ? remove from the picture layers 2/3 and 4, by hyperpolarizing the respective axons
           go to 9901
             do L = 1, num_suppyrRS
              curr_suppyrRS(numcomp_suppyrRS,L) = -0.25d0
             end do

             do L = 1, num_suppyrFRB
              curr_suppyrFRB(numcomp_suppyrFRB,L) = -0.25d0
             end do

             do L = 1, num_supbask  
              curr_supbask  (numcomp_supbask  ,L) = -0.25d0
             end do

             do L = 1, num_supng    
              curr_supng    (numcomp_supng    ,L) = -0.25d0
             end do

             do L = 1, num_supaxax  
              curr_supaxax  (numcomp_supaxax  ,L) = -0.25d0
             end do

             do L = 1, num_supLTS   
              curr_supLTS   (numcomp_supLTS   ,L) = -0.25d0
             end do

             do L = 1, num_spinstell
              curr_spinstell(numcomp_spinstell,L) = -0.25d0
             end do
9901       continue

! ? remove from the picture layers 5 and 6, by hyperpolarizing the respective axons
           go to 9902
             do L = 1, num_tuftIB   
              curr_tuftIB   (numcomp_tuftIB   ,L) = -0.25d0
             end do

             do L = 1, num_tuftRS   
              curr_tuftRS   (numcomp_tuftRS   ,L) = -0.25d0
             end do

             do L = 1, num_nontuftRS
              curr_nontuftRS(numcomp_nontuftRS,L) = -0.25d0
             end do

             do L = 1, num_deepbask 
              curr_deepbask (numcomp_deepbask ,L) = -0.25d0
             end do

             do L = 1, num_deepng   
              curr_deepng   (numcomp_deepng   ,L) = -0.25d0
             end do

             do L = 1, num_deepaxax 
              curr_deepaxax (numcomp_deepaxax ,L) = -0.25d0
             end do

             do L = 1, num_deepLTS  
              curr_deepLTS  (numcomp_deepLTS  ,L) = -0.25d0
             end do
9902       continue

       seed = 137.d0

       O = 0
       time = 0.d0

c CODE BELOW FOR "PICROTOXIN": scale all GABA-A
!        GOTO 30
c        z1 = 1.25d0  ! for intracortical IPSCs ! groucho84 uses 1.25
!        z1 = 1.00d0  ! for intracortical IPSCs
c        z1ap = 1.50d0  ! for intracortical IPSCs, basket, ng, and axoaxonal->princ.
         z1ap = 0.50d0  ! for intracortical IPSCs, basket, ng, and axoaxonal->princ.
c        z1ap = 3.00d0  ! for intracortical IPSCs, basket, ng, and axoaxonal->princ.
c        z1ai = 1.00d0  ! for intracortical IPSCs, basket and ng->inh.
         z1ai = 0.10d0  ! for intracortical IPSCs, basket and ng->inh.
!        z1b = 1.00d0  ! for intracortical IPSCs, LTS 
         z1bp =  0.05d0  ! for intracortical IPSCs, LTS -> princ. ! groucho84 in effect uses 1.25 
         z1bi = 1.00d0  ! for intracortical IPSCs, LTS -> inh.
         z2 = 1.00d0  ! for intrathalamic IPSCs, usual 1.00
      gGABA_supbask_to_suppyrRS   =  z1ap * gGABA_supbask_to_suppyrRS
      gGABA_supbask_to_suppyrFRB  =  z1ap * gGABA_supbask_to_suppyrFRB
      gGABA_supbask_to_supbask    =  z1ai * gGABA_supbask_to_supbask
      gGABA_supbask_to_supng      =  z1ai * gGABA_supbask_to_supng  
      gGABA_supbask_to_supaxax    =  z1ai * gGABA_supbask_to_supaxax
      gGABA_supbask_to_supLTS     =  z1ai * gGABA_supbask_to_supLTS
      gGABA_supbask_to_spinstell  =  z1ap * gGABA_supbask_to_spinstell

      gGABA_supng_to_suppyrRS   =  z1ap * gGABA_supng_to_suppyrRS
      gGABA_supng_to_suppyrFRB  =  z1ap * gGABA_supng_to_suppyrFRB
      gGABA_supng_to_nontuftRS  =  z1ap * gGABA_supng_to_nontuftRS
      gGABA_supng_to_tuftIB     =  z1ap * gGABA_supng_to_tuftIB   
      gGABA_supng_to_tuftRS     =  z1ap * gGABA_supng_to_tuftRS   
      gGABA_supng_to_supng      =  z1ai * gGABA_supng_to_supng    
      gGABA_supng_to_supbask    =  z1ai * gGABA_supng_to_supbask  

      gGABA_supaxax_to_suppyrRS   =  z1ap * gGABA_supaxax_to_suppyrRS
      gGABA_supaxax_to_suppyrFRB  =  z1ap * gGABA_supaxax_to_suppyrFRB
      gGABA_supaxax_to_spinstell  =  z1ap * gGABA_supaxax_to_spinstell
      gGABA_supaxax_to_tuftIB     =  z1ap * gGABA_supaxax_to_tuftIB
      gGABA_supaxax_to_tuftRS     =  z1ap * gGABA_supaxax_to_tuftRS
      gGABA_supaxax_to_nontuftRS  =  z1ap * gGABA_supaxax_to_nontuftRS

      gGABA_supLTS_to_suppyrRS    =  z1bp * gGABA_supLTS_to_suppyrRS
      gGABA_supLTS_to_suppyrFRB   =  z1bp * gGABA_supLTS_to_suppyrFRB
      gGABA_supLTS_to_supbask     =  z1bi * gGABA_supLTS_to_supbask
      gGABA_supLTS_to_supaxax     =  z1bi * gGABA_supLTS_to_supaxax
      gGABA_supLTS_to_supLTS      =  z1bi * gGABA_supLTS_to_supLTS
      gGABA_supLTS_to_spinstell   =  z1bp * gGABA_supLTS_to_spinstell
      gGABA_supLTS_to_tuftIB      =  z1bp * gGABA_supLTS_to_tuftIB
      gGABA_supLTS_to_tuftRS      =  z1bp * gGABA_supLTS_to_tuftRS
      gGABA_supLTS_to_deepbask    =  z1bi * gGABA_supLTS_to_deepbask
      gGABA_supLTS_to_deepaxax    =  z1bi * gGABA_supLTS_to_deepaxax
      gGABA_supLTS_to_deepLTS     =  z1bi * gGABA_supLTS_to_deepLTS
      gGABA_supLTS_to_nontuftRS   =  z1bp * gGABA_supLTS_to_nontuftRS

      gGABA_deepbask_to_spinstell =  z1ap * gGABA_deepbask_to_spinstell
      gGABA_deepbask_to_tuftIB    =  z1ap * gGABA_deepbask_to_tuftIB
      gGABA_deepbask_to_tuftRS    =  z1ap * gGABA_deepbask_to_tuftRS
      gGABA_deepbask_to_deepbask  =  z1ai * gGABA_deepbask_to_deepbask
      gGABA_deepbask_to_deepng    =  z1ai * gGABA_deepbask_to_deepng  
      gGABA_deepbask_to_deepaxax  =  z1ai * gGABA_deepbask_to_deepaxax
      gGABA_deepbask_to_deepLTS   =  z1ai * gGABA_deepbask_to_deepLTS
      gGABA_deepbask_to_nontuftRS =  z1ap * gGABA_deepbask_to_nontuftRS

      gGABA_deepng_to_tuftIB    =  z1ap * gGABA_deepng_to_tuftIB
      gGABA_deepng_to_tuftRS    =  z1ap * gGABA_deepng_to_tuftRS
      gGABA_deepng_to_nontuftRS =  z1ap * gGABA_deepng_to_nontuftRS
      gGABA_deepng_to_spinstell =  z1ap * gGABA_deepng_to_spinstell
      gGABA_deepng_to_deepng    =  z1ai * gGABA_deepng_to_deepng   
      gGABA_deepng_to_deepbask  =  z1ai * gGABA_deepng_to_deepbask 

      gGABA_deepaxax_to_suppyrRS   = z1ap * gGABA_deepaxax_to_suppyrRS
      gGABA_deepaxax_to_suppyrFRB  = z1ap * gGABA_deepaxax_to_suppyrFRB
      gGABA_deepaxax_to_spinstell  = z1ap * gGABA_deepaxax_to_spinstell
      gGABA_deepaxax_to_tuftIB     = z1ap * gGABA_deepaxax_to_tuftIB
      gGABA_deepaxax_to_tuftRS     = z1ap * gGABA_deepaxax_to_tuftRS
      gGABA_deepaxax_to_nontuftRS  = z1ap * gGABA_deepaxax_to_nontuftRS

      gGABA_deepLTS_to_suppyrRS    = z1bp * gGABA_deepLTS_to_suppyrRS
      gGABA_deepLTS_to_suppyrFRB   = z1bp * gGABA_deepLTS_to_suppyrFRB
      gGABA_deepLTS_to_supbask     = z1bi * gGABA_deepLTS_to_supbask
      gGABA_deepLTS_to_supaxax     = z1bi * gGABA_deepLTS_to_supaxax
      gGABA_deepLTS_to_supLTS      = z1bi * gGABA_deepLTS_to_supLTS
      gGABA_deepLTS_to_spinstell   = z1bp * gGABA_deepLTS_to_spinstell
      gGABA_deepLTS_to_tuftIB      = z1bp * gGABA_deepLTS_to_tuftIB
      gGABA_deepLTS_to_tuftRS      = z1bp * gGABA_deepLTS_to_tuftRS
      gGABA_deepLTS_to_deepbask    = z1bi * gGABA_deepLTS_to_deepbask
      gGABA_deepLTS_to_deepaxax    = z1bi * gGABA_deepLTS_to_deepaxax
      gGABA_deepLTS_to_deepLTS     = z1bi * gGABA_deepLTS_to_deepLTS
      gGABA_deepLTS_to_nontuftRS   = z1bp * gGABA_deepLTS_to_nontuftRS

        do L = 1, num_nRT
      gGABA_nRT_to_TCR(L)          = z2 * gGABA_nRT_to_TCR(L) 
        end do
      gGABA_nRT_to_nRT             = z2 * gGABA_nRT_to_nRT
30          CONTINUE
c End "PICROTOXIN" code

c "CGP" section: scale GABA-B conductances.
c            z7 = 1.5d0
             z7 = 3.0d0
!            z7 = 0.d0
       gGABAB_supng_to_suppyrRS  = z7*gGABAB_supng_to_suppyrRS
       gGABAB_supng_to_suppyrFRB = z7*gGABAB_supng_to_suppyrFRB
       gGABAB_supng_to_nontuftRS = z7*gGABAB_supng_to_nontuftRS
       gGABAB_supng_to_tuftIB    = z7*gGABAB_supng_to_tuftIB   
       gGABAB_supng_to_tuftRS    = z7*gGABAB_supng_to_tuftRS   

       gGABAB_deepng_to_tuftIB   = z7*gGABAB_deepng_to_tuftIB
       gGABAB_deepng_to_tuftRS   = z7*gGABAB_deepng_to_tuftRS
       gGABAB_deepng_to_nontuftRS= z7*gGABAB_deepng_to_nontuftRS
       gGABAB_deepng_to_spinstell= z7*gGABAB_deepng_to_spinstell
31           CONTINUE
c End CGP section

! Code below is "NBQX": scale all AMPA; see also below for possibility
! of further additional scaling of connections between layers
!        GOTO 35
!          z1 = 1.00d0  ! intracortical e/i
c          z1 = 2.00d0  ! intracortical e/i ! usual 1.00; use 2.0 for delta78
           z1 = 0.10d0  ! intracortical e/i ! usual 1.00; use 2.0 for delta78
           z3 = 0.00d0  ! TCR -> cortical i ! usual 1.0
           z4 = 1.00d0  ! TCR -> nRT & nontuftRS ->nRT ! usual 1.00
c          z5 = 0.25d0  ! spinstell -> spinstell; may reduce as in spindle series, 8 May'04
           z5 = 0.10d0  ! spinstell -> spinstell; may reduce as in spindle series, 8 May'04
!          z5 = 0.00d0  ! spinstell -> spinstell; may reduce as in spindle series, 8 May'04
c          z6 = 0.10d0    ! layer 5 tuftIB or RS -> layer 5 tuftIB or RS
           z6 = 0.05d0    ! layer 5 tuftIB or RS -> layer 5 tuftIB or RS; 16 May '07, further reduction
c          z6 = 1.00d0    ! layer 5 tuftIB or RS -> layer 5 tuftIB or RS
!          z2 = 2.00d0  ! everything else; note that this may be INCREASED, usual 1.0
c In groucho84, z2 = 2.00d0
           z2 = 0.20d0  ! everything else; note that this may be INCREASED, usual 1.0

      gAMPA_suppyrRS_to_suppyrRS= z2 * gAMPA_suppyrRS_to_suppyrRS
      gAMPA_suppyrRS_to_suppyrFRB= z2 * gAMPA_suppyrRS_to_suppyrFRB
      gAMPA_suppyrRS_to_supbask  = z1 * gAMPA_suppyrRS_to_supbask
      gAMPA_suppyrRS_to_supng    = z1 * gAMPA_suppyrRS_to_supng  
      gAMPA_suppyrRS_to_supaxax  = z1 * gAMPA_suppyrRS_to_supaxax
      gAMPA_suppyrRS_to_supLTS   = z1 * gAMPA_suppyrRS_to_supLTS
      gAMPA_suppyrRS_to_spinstell= z2 * gAMPA_suppyrRS_to_spinstell
      gAMPA_suppyrRS_to_tuftIB   = z2 * gAMPA_suppyrRS_to_tuftIB
      gAMPA_suppyrRS_to_tuftRS   = z2 * gAMPA_suppyrRS_to_tuftRS
      gAMPA_suppyrRS_to_deepbask = z1 * gAMPA_suppyrRS_to_deepbask
      gAMPA_suppyrRS_to_deepaxax = z1 * gAMPA_suppyrRS_to_deepaxax
      gAMPA_suppyrRS_to_deepLTS  = z1 * gAMPA_suppyrRS_to_deepLTS
      gAMPA_suppyrRS_to_nontuftRS= z2 * gAMPA_suppyrRS_to_nontuftRS

      gAMPA_suppyrFRB_to_suppyrRS= z2 * gAMPA_suppyrFRB_to_suppyrRS
      gAMPA_suppyrFRB_to_suppyrFRB=z2 * gAMPA_suppyrFRB_to_suppyrFRB
      gAMPA_suppyrFRB_to_supbask  =z1 * gAMPA_suppyrFRB_to_supbask
      gAMPA_suppyrFRB_to_supng    =z1 * gAMPA_suppyrFRB_to_supng  
      gAMPA_suppyrFRB_to_supaxax  =z1 * gAMPA_suppyrFRB_to_supaxax
      gAMPA_suppyrFRB_to_supLTS   =z1 * gAMPA_suppyrFRB_to_supLTS
      gAMPA_suppyrFRB_to_spinstell=z2 * gAMPA_suppyrFRB_to_spinstell
      gAMPA_suppyrFRB_to_tuftIB   =z2 * gAMPA_suppyrFRB_to_tuftIB
      gAMPA_suppyrFRB_to_tuftRS   =z2 * gAMPA_suppyrFRB_to_tuftRS
      gAMPA_suppyrFRB_to_deepbask =z1 * gAMPA_suppyrFRB_to_deepbask
      gAMPA_suppyrFRB_to_deepaxax =z1 * gAMPA_suppyrFRB_to_deepaxax
      gAMPA_suppyrFRB_to_deepLTS  =z1 * gAMPA_suppyrFRB_to_deepLTS
      gAMPA_suppyrFRB_to_nontuftRS=z2 * gAMPA_suppyrFRB_to_nontuftRS

      gAMPA_spinstell_to_suppyrRS = z2 * gAMPA_spinstell_to_suppyrRS
      gAMPA_spinstell_to_suppyrFRB= z2 * gAMPA_spinstell_to_suppyrFRB
      gAMPA_spinstell_to_supbask  = z1 * gAMPA_spinstell_to_supbask
      gAMPA_spinstell_to_supaxax  = z1 * gAMPA_spinstell_to_supaxax
      gAMPA_spinstell_to_supLTS   = z1 * gAMPA_spinstell_to_supLTS
      gAMPA_spinstell_to_spinstell= z5 * gAMPA_spinstell_to_spinstell
      gAMPA_spinstell_to_tuftIB   = z2 * gAMPA_spinstell_to_tuftIB
      gAMPA_spinstell_to_tuftRS   = z2 * gAMPA_spinstell_to_tuftRS
      gAMPA_spinstell_to_deepbask = z1 * gAMPA_spinstell_to_deepbask
      gAMPA_spinstell_to_deepng   = z1 * gAMPA_spinstell_to_deepng  
      gAMPA_spinstell_to_deepaxax = z1 * gAMPA_spinstell_to_deepaxax
      gAMPA_spinstell_to_deepLTS  = z1 * gAMPA_spinstell_to_deepLTS
      gAMPA_spinstell_to_nontuftRS= z2 * gAMPA_spinstell_to_nontuftRS

      gAMPA_tuftIB_to_suppyrRS    = z2 * gAMPA_tuftIB_to_suppyrRS
      gAMPA_tuftIB_to_suppyrFRB   = z2 * gAMPA_tuftIB_to_suppyrFRB
      gAMPA_tuftIB_to_supbask     = z1 * gAMPA_tuftIB_to_supbask
      gAMPA_tuftIB_to_supaxax     = z1 * gAMPA_tuftIB_to_supaxax
      gAMPA_tuftIB_to_supLTS      = z1 * gAMPA_tuftIB_to_supLTS
      gAMPA_tuftIB_to_spinstell   = z2 * gAMPA_tuftIB_to_spinstell
      gAMPA_tuftIB_to_tuftIB      = z6 * gAMPA_tuftIB_to_tuftIB
      gAMPA_tuftIB_to_tuftRS      = z6 * gAMPA_tuftIB_to_tuftRS
      gAMPA_tuftIB_to_deepbask    = z1 * gAMPA_tuftIB_to_deepbask
      gAMPA_tuftIB_to_deepng      = z1 * gAMPA_tuftIB_to_deepng  
      gAMPA_tuftIB_to_deepaxax    = z1 * gAMPA_tuftIB_to_deepaxax
      gAMPA_tuftIB_to_deepLTS     = z1 * gAMPA_tuftIB_to_deepLTS
      gAMPA_tuftIB_to_nontuftRS   = z2 * gAMPA_tuftIB_to_nontuftRS

      gAMPA_tuftRS_to_suppyrRS    = z2 * gAMPA_tuftRS_to_suppyrRS
      gAMPA_tuftRS_to_suppyrFRB   = z2 * gAMPA_tuftRS_to_suppyrFRB
      gAMPA_tuftRS_to_supbask     = z1 * gAMPA_tuftRS_to_supbask 
      gAMPA_tuftRS_to_supaxax     = z1 * gAMPA_tuftRS_to_supaxax
      gAMPA_tuftRS_to_supLTS      = z1 * gAMPA_tuftRS_to_supLTS
      gAMPA_tuftRS_to_spinstell   = z2 * gAMPA_tuftRS_to_spinstell
      gAMPA_tuftRS_to_tuftIB      = z6 * gAMPA_tuftRS_to_tuftIB
      gAMPA_tuftRS_to_tuftRS      = z6 * gAMPA_tuftRS_to_tuftRS
      gAMPA_tuftRS_to_deepbask    = z1 * gAMPA_tuftRS_to_deepbask
      gAMPA_tuftRS_to_deepng      = z1 * gAMPA_tuftRS_to_deepng  
      gAMPA_tuftRS_to_deepaxax    = z1 * gAMPA_tuftRS_to_deepaxax
      gAMPA_tuftRS_to_deepLTS     = z1 * gAMPA_tuftRS_to_deepLTS
      gAMPA_tuftRS_to_nontuftRS   = z2 * gAMPA_tuftRS_to_nontuftRS

      gAMPA_TCR_to_suppyrRS        = z2 * gAMPA_TCR_to_suppyrRS
      gAMPA_TCR_to_suppyrFRB       = z2 * gAMPA_TCR_to_suppyrFRB
      gAMPA_TCR_to_supbask         = z3 * gAMPA_TCR_to_supbask
      gAMPA_TCR_to_supng           = z3 * gAMPA_TCR_to_supng   
      gAMPA_TCR_to_supaxax         = z3 * gAMPA_TCR_to_supaxax
      gAMPA_TCR_to_spinstell       = z2 * gAMPA_TCR_to_spinstell
      gAMPA_TCR_to_tuftIB          = z2 * gAMPA_TCR_to_tuftIB
      gAMPA_TCR_to_tuftRS          = z2 * gAMPA_TCR_to_tuftRS
      gAMPA_TCR_to_deepbask        = z3 * gAMPA_TCR_to_deepbask
      gAMPA_TCR_to_deepng          = z3 * gAMPA_TCR_to_deepng  
      gAMPA_TCR_to_deepaxax        = z3 * gAMPA_TCR_to_deepaxax
      gAMPA_TCR_to_nRT             = z4 * gAMPA_TCR_to_nRT
      gAMPA_TCR_to_nontuftRS       = z2 * gAMPA_TCR_to_nontuftRS

      gAMPA_nontuftRS_to_suppyrRS  = z2 * gAMPA_nontuftRS_to_suppyrRS
      gAMPA_nontuftRS_to_suppyrFRB = z2 * gAMPA_nontuftRS_to_suppyrFRB
      gAMPA_nontuftRS_to_supbask   = z1 * gAMPA_nontuftRS_to_supbask
      gAMPA_nontuftRS_to_supaxax   = z1 * gAMPA_nontuftRS_to_supaxax
      gAMPA_nontuftRS_to_supLTS    = z1 * gAMPA_nontuftRS_to_supLTS
      gAMPA_nontuftRS_to_spinstell = z2 * gAMPA_nontuftRS_to_spinstell
      gAMPA_nontuftRS_to_tuftIB    = z2 * gAMPA_nontuftRS_to_tuftIB
      gAMPA_nontuftRS_to_tuftRS    = z2 * gAMPA_nontuftRS_to_tuftRS
      gAMPA_nontuftRS_to_deepbask  = z1 * gAMPA_nontuftRS_to_deepbask
      gAMPA_nontuftRS_to_deepng    = z1 * gAMPA_nontuftRS_to_deepng  
      gAMPA_nontuftRS_to_deepaxax  = z1 * gAMPA_nontuftRS_to_deepaxax
      gAMPA_nontuftRS_to_deepLTS   = z1 * gAMPA_nontuftRS_to_deepLTS
      gAMPA_nontuftRS_to_TCR       = z2 * gAMPA_nontuftRS_to_TCR
      gAMPA_nontuftRS_to_nRT       = z4 * gAMPA_nontuftRS_to_nRT
      gAMPA_nontuftRS_to_nontuftRS = z2 * gAMPA_nontuftRS_to_nontuftRS

! Code below: allows for
! further additional scaling of connections between layers
      z10a = 1.d0 ! scales deep pyramids to superficial excitatory cells
      z10b = 0.5d0 ! scales deep pyramids to superficial inhibitory cells
      z12 = 1.d0 ! scales superficial pyramids to deep cells & spinstell

      gAMPA_suppyrRS_to_spinstell= z12* gAMPA_suppyrRS_to_spinstell
      gAMPA_suppyrRS_to_tuftIB   = z12* gAMPA_suppyrRS_to_tuftIB
      gAMPA_suppyrRS_to_tuftRS   = z12* gAMPA_suppyrRS_to_tuftRS
      gAMPA_suppyrRS_to_deepbask = z12* gAMPA_suppyrRS_to_deepbask
      gAMPA_suppyrRS_to_deepaxax = z12* gAMPA_suppyrRS_to_deepaxax
      gAMPA_suppyrRS_to_deepLTS  = z12* gAMPA_suppyrRS_to_deepLTS
      gAMPA_suppyrRS_to_nontuftRS= z12* gAMPA_suppyrRS_to_nontuftRS

      gAMPA_suppyrFRB_to_spinstell=z12* gAMPA_suppyrFRB_to_spinstell
      gAMPA_suppyrFRB_to_tuftIB   =z12* gAMPA_suppyrFRB_to_tuftIB
      gAMPA_suppyrFRB_to_tuftRS   =z12* gAMPA_suppyrFRB_to_tuftRS
      gAMPA_suppyrFRB_to_deepbask =z12* gAMPA_suppyrFRB_to_deepbask
      gAMPA_suppyrFRB_to_deepaxax =z12* gAMPA_suppyrFRB_to_deepaxax
      gAMPA_suppyrFRB_to_deepLTS  =z12* gAMPA_suppyrFRB_to_deepLTS
      gAMPA_suppyrFRB_to_nontuftRS=z12* gAMPA_suppyrFRB_to_nontuftRS


      gAMPA_tuftIB_to_suppyrRS    = z10a* gAMPA_tuftIB_to_suppyrRS
      gAMPA_tuftIB_to_suppyrFRB   = z10a* gAMPA_tuftIB_to_suppyrFRB
      gAMPA_tuftIB_to_supbask     = z10b* gAMPA_tuftIB_to_supbask
      gAMPA_tuftIB_to_supaxax     = z10b* gAMPA_tuftIB_to_supaxax
      gAMPA_tuftIB_to_supLTS      = z10b* gAMPA_tuftIB_to_supLTS

      gAMPA_tuftRS_to_suppyrRS    = z10a* gAMPA_tuftRS_to_suppyrRS
      gAMPA_tuftRS_to_suppyrFRB   = z10a* gAMPA_tuftRS_to_suppyrFRB
      gAMPA_tuftRS_to_supbask     = z10b* gAMPA_tuftRS_to_supbask 
      gAMPA_tuftRS_to_supaxax     = z10b* gAMPA_tuftRS_to_supaxax
      gAMPA_tuftRS_to_supLTS      = z10b* gAMPA_tuftRS_to_supLTS

      gAMPA_nontuftRS_to_suppyrRS  = z10a* gAMPA_nontuftRS_to_suppyrRS
      gAMPA_nontuftRS_to_suppyrFRB = z10a* gAMPA_nontuftRS_to_suppyrFRB
      gAMPA_nontuftRS_to_supbask   = z10b* gAMPA_nontuftRS_to_supbask
      gAMPA_nontuftRS_to_supaxax   = z10b* gAMPA_nontuftRS_to_supaxax
      gAMPA_nontuftRS_to_supLTS    = z10b* gAMPA_nontuftRS_to_supLTS
35         CONTINUE
c End "NBQX" section.

c Code below scales TCR output to cortex (not to nRT), AMPA & NMDA
!     goto 60
c      z = 1.d0
       z = 0.d0

      gAMPA_TCR_to_suppyrRS = z * gAMPA_TCR_to_suppyrRS
      gNMDA_TCR_to_suppyrRS = z * gNMDA_TCR_to_suppyrRS
      gAMPA_TCR_to_suppyrFRB = z * gAMPA_TCR_to_suppyrFRB
      gNMDA_TCR_to_suppyrFRB = z * gNMDA_TCR_to_suppyrFRB
      gAMPA_TCR_to_supbask = z * gAMPA_TCR_to_supbask
      gNMDA_TCR_to_supbask = z * gNMDA_TCR_to_supbask
      gAMPA_TCR_to_supaxax = z * gAMPA_TCR_to_supaxax
      gNMDA_TCR_to_supaxax = z * gNMDA_TCR_to_supaxax
      gAMPA_TCR_to_spinstell = z * gAMPA_TCR_to_spinstell
      gNMDA_TCR_to_spinstell = z * gNMDA_TCR_to_spinstell
      gAMPA_TCR_to_tuftIB = z * gAMPA_TCR_to_tuftIB
      gNMDA_TCR_to_tuftIB = z * gNMDA_TCR_to_tuftIB
      gAMPA_TCR_to_tuftRS = z * gAMPA_TCR_to_tuftRS
      gNMDA_TCR_to_tuftRS = z * gNMDA_TCR_to_tuftRS
      gAMPA_TCR_to_deepbask = z * gAMPA_TCR_to_deepbask
      gNMDA_TCR_to_deepbask = z * gNMDA_TCR_to_deepbask
      gAMPA_TCR_to_deepaxax = z * gAMPA_TCR_to_deepaxax
      gNMDA_TCR_to_deepaxax = z * gNMDA_TCR_to_deepaxax
      gAMPA_TCR_to_nontuftRS = z * gAMPA_TCR_to_nontuftRS
      gNMDA_TCR_to_nontuftRS = z * gNMDA_TCR_to_nontuftRS

60          CONTINUE

c Code below scales some/all NMDA conductances. APV.
!        GOTO 40
c        z1 = 0.0d0 ! to interneurons
         z1 = 0.5d0 ! to interneurons
c        z1 = 0.0d0 ! to interneurons; = 1 in delta85
! Usual scaling of NMDA to princ. cells, including FRB, is 0.5
         z2 = 0.25d0 ! to  cort. principal cells, except FRB
c        z2 = 2.0d0 ! to  cort. principal cells, except FRB; = 1 for delta79
c        z2 = 0.0d0 ! to  cort. principal cells, except FRB; = 1 for delta79,= 2 for delta85
         z3 = 0.0d0 ! to suppyrFRB
         z4 = 0.0d0  ! to TCR and nRT and from TCR to cort. princ.
         z5 = 1.05d0  ! spinstell to spinstell
c        z5 = 0.00d0  ! spinstell to spinstell
      gNMDA_suppyrRS_to_suppyrRS= z2 *
     &  gNMDA_suppyrRS_to_suppyrRS
      gNMDA_suppyrRS_to_suppyrFRB= z3 *
     &  gNMDA_suppyrRS_to_suppyrFRB
      gNMDA_suppyrRS_to_supbask  = z1 *
     &  gNMDA_suppyrRS_to_supbask
      gNMDA_suppyrRS_to_supng    = z1 *
     &  gNMDA_suppyrRS_to_supng  
      gNMDA_suppyrRS_to_supaxax  = z1 *
     &  gNMDA_suppyrRS_to_supaxax
      gNMDA_suppyrRS_to_supLTS   = z1 *
     &  gNMDA_suppyrRS_to_supLTS   
      gNMDA_suppyrRS_to_spinstell= z2 *
     &  gNMDA_suppyrRS_to_spinstell
      gNMDA_suppyrRS_to_tuftIB   = z2 *
     &  gNMDA_suppyrRS_to_tuftIB
      gNMDA_suppyrRS_to_tuftRS   = z2 *
     &  gNMDA_suppyrRS_to_tuftRS  
      gNMDA_suppyrRS_to_deepbask = z1 *
     &  gNMDA_suppyrRS_to_deepbask 
      gNMDA_suppyrRS_to_deepaxax = z1 *
     &  gNMDA_suppyrRS_to_deepaxax
      gNMDA_suppyrRS_to_deepLTS  = z1 *
     &  gNMDA_suppyrRS_to_deepLTS 
      gNMDA_suppyrRS_to_nontuftRS= z2 *
     &  gNMDA_suppyrRS_to_nontuftRS
 
      gNMDA_suppyrFRB_to_suppyrRS= z2 *
     &  gNMDA_suppyrFRB_to_suppyrRS
      gNMDA_suppyrFRB_to_suppyrFRB= z3 *
     &  gNMDA_suppyrFRB_to_suppyrFRB
      gNMDA_suppyrFRB_to_supbask  = z1 *
     &  gNMDA_suppyrFRB_to_suppyrFRB
      gNMDA_suppyrFRB_to_supng    = z1 *
     &  gNMDA_suppyrFRB_to_supng   
      gNMDA_suppyrFRB_to_supaxax  = z1 *
     &  gNMDA_suppyrFRB_to_supaxax 
      gNMDA_suppyrFRB_to_supLTS   = z1 *
     &  gNMDA_suppyrFRB_to_supLTS
      gNMDA_suppyrFRB_to_spinstell= z2 *
     &  gNMDA_suppyrFRB_to_spinstell
      gNMDA_suppyrFRB_to_tuftIB   = z2 *
     &  gNMDA_suppyrFRB_to_tuftIB
      gNMDA_suppyrFRB_to_tuftRS   = z2 *
     &  gNMDA_suppyrFRB_to_tuftRS   
      gNMDA_suppyrFRB_to_deepbask = z1 *
     &  gNMDA_suppyrFRB_to_deepbask
      gNMDA_suppyrFRB_to_deepaxax = z1 *
     &  gNMDA_suppyrFRB_to_deepaxax
      gNMDA_suppyrFRB_to_deepLTS  = z1 *
     &  gNMDA_suppyrFRB_to_deepLTS
      gNMDA_suppyrFRB_to_nontuftRS= z2 *
     &  gNMDA_suppyrFRB_to_nontuftRS

      gNMDA_spinstell_to_suppyrRS = z2 *
     &  gNMDA_spinstell_to_suppyrRS
      gNMDA_spinstell_to_suppyrFRB= z3 *
     &  gNMDA_spinstell_to_suppyrFRB
      gNMDA_spinstell_to_supbask  = z1 *
     &  gNMDA_spinstell_to_supbask
      gNMDA_spinstell_to_supaxax  = z1 *
     &  gNMDA_spinstell_to_supaxax 
      gNMDA_spinstell_to_supLTS   = z1 *
     &  gNMDA_spinstell_to_supLTS
      gNMDA_spinstell_to_spinstell= z5 *
     &  gNMDA_spinstell_to_spinstell
      gNMDA_spinstell_to_tuftIB   = z2 *
     &  gNMDA_spinstell_to_tuftIB 
      gNMDA_spinstell_to_tuftRS   = z2 *
     &  gNMDA_spinstell_to_tuftRS 
      gNMDA_spinstell_to_deepbask = z1 *
     &  gNMDA_spinstell_to_deepbask 
      gNMDA_spinstell_to_deepng   = z1 *
     &  gNMDA_spinstell_to_deepng   
      gNMDA_spinstell_to_deepaxax = z1 *
     &  gNMDA_spinstell_to_deepaxax
      gNMDA_spinstell_to_deepLTS  = z1 *
     &  gNMDA_spinstell_to_deepLTS 
      gNMDA_spinstell_to_nontuftRS= z2 *
     &  gNMDA_spinstell_to_nontuftRS

      gNMDA_tuftIB_to_suppyrRS    = z2 *
     &  gNMDA_tuftIB_to_suppyrRS 
      gNMDA_tuftIB_to_suppyrFRB   = z3 *
     &  gNMDA_tuftIB_to_suppyrFRB
      gNMDA_tuftIB_to_supbask     = z1 *
     &  gNMDA_tuftIB_to_supbask 
      gNMDA_tuftIB_to_supaxax     = z1 *
     &  gNMDA_tuftIB_to_supaxax 
      gNMDA_tuftIB_to_supLTS      = z1 *
     &  gNMDA_tuftIB_to_supLTS 
      gNMDA_tuftIB_to_spinstell   = z2 *
     &  gNMDA_tuftIB_to_spinstell 
      gNMDA_tuftIB_to_tuftIB      = z2 *
     &  gNMDA_tuftIB_to_tuftIB 
      gNMDA_tuftIB_to_tuftRS      = z2 *
     &  gNMDA_tuftIB_to_tuftRS  
      gNMDA_tuftIB_to_deepbask    = z1 *
     &  gNMDA_tuftIB_to_deepbask 
      gNMDA_tuftIB_to_deepng      = z1 *
     &  gNMDA_tuftIB_to_deepng   
      gNMDA_tuftIB_to_deepaxax    = z1 *
     &  gNMDA_tuftIB_to_deepaxax
      gNMDA_tuftIB_to_deepLTS     = z1 *
     &  gNMDA_tuftIB_to_deepLTS   
      gNMDA_tuftIB_to_nontuftRS   = z2 *
     &  gNMDA_tuftIB_to_nontuftRS  

      gNMDA_tuftRS_to_suppyrRS    = z2 *
     &  gNMDA_tuftRS_to_suppyrRS
      gNMDA_tuftRS_to_suppyrFRB   = z3 *
     &  gNMDA_tuftRS_to_suppyrFRB
      gNMDA_tuftRS_to_supbask     = z1 *
     &  gNMDA_tuftRS_to_supbask 
      gNMDA_tuftRS_to_supaxax     = z1 *
     &  gNMDA_tuftRS_to_supaxax 
      gNMDA_tuftRS_to_supLTS      = z1 *
     &  gNMDA_tuftRS_to_supLTS   
      gNMDA_tuftRS_to_spinstell   = z2 *
     &  gNMDA_tuftRS_to_spinstell  
      gNMDA_tuftRS_to_tuftIB      = z2 *
     &  gNMDA_tuftRS_to_tuftIB
      gNMDA_tuftRS_to_tuftRS      = z2 *
     &  gNMDA_tuftRS_to_tuftRS 
      gNMDA_tuftRS_to_deepbask    = z1 *
     &  gNMDA_tuftRS_to_deepbask
      gNMDA_tuftRS_to_deepng      = z1 *
     &  gNMDA_tuftRS_to_deepng  
      gNMDA_tuftRS_to_deepaxax    = z1 *
     &  gNMDA_tuftRS_to_deepaxax
      gNMDA_tuftRS_to_deepLTS     = z1 *
     &  gNMDA_tuftRS_to_deepLTS   
      gNMDA_tuftRS_to_nontuftRS   = z2 *
     &  gNMDA_tuftRS_to_nontuftRS 

      gNMDA_TCR_to_suppyrRS        = z4 *
     &  gNMDA_TCR_to_suppyrRS 
      gNMDA_TCR_to_suppyrFRB       = z4 *
     &  gNMDA_TCR_to_suppyrFRB 
      gNMDA_TCR_to_supbask         = z1 *
     &  gNMDA_TCR_to_supbask
      gNMDA_TCR_to_supng           = z1 *
     &  gNMDA_TCR_to_supng  
      gNMDA_TCR_to_supaxax         = z1 *
     &  gNMDA_TCR_to_supaxax 
      gNMDA_TCR_to_spinstell       = z4 *
     &  gNMDA_TCR_to_spinstell 
      gNMDA_TCR_to_tuftIB          = z4 *
     &  gNMDA_TCR_to_tuftIB   
      gNMDA_TCR_to_tuftRS          = z4 *
     &  gNMDA_TCR_to_tuftRS 
      gNMDA_TCR_to_deepbask        = z1 *
     &  gNMDA_TCR_to_deepbask 
      gNMDA_TCR_to_deepng          = z1 *
     &  gNMDA_TCR_to_deepng   
      gNMDA_TCR_to_deepaxax        = z1 *
     &  gNMDA_TCR_to_deepaxax 
      gNMDA_TCR_to_nRT             = z1 *
     &  gNMDA_TCR_to_nRT  
      gNMDA_TCR_to_nontuftRS       = z4 *
     &  gNMDA_TCR_to_nontuftRS  

      gNMDA_nontuftRS_to_suppyrRS  = z2 *
     &  gNMDA_nontuftRS_to_suppyrRS
      gNMDA_nontuftRS_to_suppyrFRB = z3 *
     & gNMDA_nontuftRS_to_suppyrFRB 
      gNMDA_nontuftRS_to_supbask   = z1 *
     &  gNMDA_nontuftRS_to_supbask 
      gNMDA_nontuftRS_to_supaxax   = z1 *
     &  gNMDA_nontuftRS_to_supaxax  
      gNMDA_nontuftRS_to_supLTS    = z1 *
     &  gNMDA_nontuftRS_to_supLTS 
      gNMDA_nontuftRS_to_spinstell = z2 *
     &  gNMDA_nontuftRS_to_spinstell 
      gNMDA_nontuftRS_to_tuftIB    = z2 *
     &  gNMDA_nontuftRS_to_tuftIB 
      gNMDA_nontuftRS_to_tuftRS    = z2 *
     &  gNMDA_nontuftRS_to_tuftRS 
      gNMDA_nontuftRS_to_deepbask  = z1 *
     & gNMDA_nontuftRS_to_deepbask
      gNMDA_nontuftRS_to_deepaxax  = z1 *
     &  gNMDA_nontuftRS_to_deepaxax 
      gNMDA_nontuftRS_to_deepLTS   = z1 *
     &  gNMDA_nontuftRS_to_deepLTS
      gNMDA_nontuftRS_to_TCR       = z4 *
     &  gNMDA_nontuftRS_to_TCR 
      gNMDA_nontuftRS_to_nRT       = z4 *
     &  gNMDA_nontuftRS_to_nRT    
      gNMDA_nontuftRS_to_nontuftRS = z2 *
     &  gNMDA_nontuftRS_to_nontuftRS 
40    CONTINUE
c End section scaling all NMDA conductances.       

c INITIALIZE ALL THE INTEGRATION SUBROUTINES
        initialize = 0
        firstcell = 1
        lastcell =  1
      IF (nodecell(thisno).eq.'suppyrRS ') then
       CALL INTEGRATE_suppyrRSXPB (O, time, num_suppyrRS,
     &    V_suppyrRS, curr_suppyrRS,
     &    initialize, firstcell, lastcell,
     & gAMPA_suppyrRS, gNMDA_suppyrRS, gGABA_A_suppyrRS,
     & gGABA_B_suppyrRS, Mg, 
     & gapcon_suppyrRS  ,totaxgj_suppyrRS   ,gjtable_suppyrRS, dt,
     & totaxgj_suppyr, gjtable_suppyr, num_suppyrFRB,
!    & vax_suppyrFRB,
     & v_suppyrFRB,
     &  chi_suppyrRS,mnaf_suppyrRS,mnap_suppyrRS,
     &  hnaf_suppyrRS,mkdr_suppyrRS,mka_suppyrRS,
     &  hka_suppyrRS,mk2_suppyrRS,hk2_suppyrRS,
     &  mkm_suppyrRS,mkc_suppyrRS,mkahp_suppyrRS,
     &  mcat_suppyrRS,hcat_suppyrRS,mcal_suppyrRS,
     &  mar_suppyrRS,field_1mm   ,field_2mm,rel_axonshift_suppyrRS)

      ELSE if (nodecell(thisno).eq.'suppyrFRB') then
       CALL INTEGRATE_suppyrFRBxPB (O, time, num_suppyrFRB,
     &    V_suppyrFRB, curr_suppyrFRB,
     &    initialize, firstcell, lastcell, 
     & gAMPA_suppyrFRB, gNMDA_suppyrFRB, gGABA_A_suppyrFRB,
     & gGABA_B_suppyrFRB, Mg, 
     & gapcon_suppyrFRB ,totaxgj_suppyrFRB  ,gjtable_suppyrFRB, dt,
     & totaxgj_suppyr, gjtable_suppyr, num_suppyrRS,
!    & vax_suppyrRS,
     & v_suppyrRS,
     &  chi_suppyrFRB,mnaf_suppyrFRB,mnap_suppyrFRB,
     &  hnaf_suppyrFRB,mkdr_suppyrFRB,mka_suppyrFRB,
     &  hka_suppyrFRB,mk2_suppyrFRB,hk2_suppyrFRB,
     &  mkm_suppyrFRB,mkc_suppyrFRB,mkahp_suppyrFRB,
     &  mcat_suppyrFRB,hcat_suppyrFRB,mcal_suppyrFRB,
     &  mar_suppyrFRB,field_1mm  ,field_2mm,rel_axonshift_suppyrFRB,
     & scale_suppyrFRB_gNaP)

      ELSE if (nodecell(thisno).eq.'supbask  ') then
       CALL INTEGRATE_supbaskx (O, time, num_supbask ,
     &    V_supbask , curr_supbask ,
     $    initialize, firstcell, lastcell,
     & gAMPA_supbask , gNMDA_supbask , gGABA_A_supbask ,
     & Mg, 
     & gapcon_supbask   ,totSDgj_supbask    ,gjtable_supbask , dt,
     &  chi_supbask,mnaf_supbask,mnap_supbask,
     &  hnaf_supbask,mkdr_supbask,mka_supbask,
     &  hka_supbask,mk2_supbask,hk2_supbask,
     &  mkm_supbask,mkc_supbask,mkahp_supbask,
     &  mcat_supbask,hcat_supbask,mcal_supbask,
     &  mar_supbask)

      ELSE if (nodecell(thisno).eq.'supng    ') then
       CALL INTEGRATE_supng    (O, time, num_supng   ,
     &    V_supng   , curr_supng   ,
     $    initialize, firstcell, lastcell,
     & gAMPA_supng   , gNMDA_supng   , gGABA_A_supng   ,
     & Mg, 
     & gapcon_supng     ,totSDgj_supng      ,gjtable_supng   , dt,
     &  chi_supng  ,mnaf_supng  ,mnap_supng  ,
     &  hnaf_supng  ,mkdr_supng  ,mka_supng  ,
     &  hka_supng  ,mk2_supng  ,hk2_supng  ,
     &  mkm_supng  ,mkc_supng  ,mkahp_supng  ,
     &  mcat_supng  ,hcat_supng  ,mcal_supng  ,
     &  mar_supng  )

      ELSE if (nodecell(thisno).eq.'supaxax  ') then
       CALL INTEGRATE_supaxaxx (O, time, num_supaxax ,
     &    V_supaxax , curr_supaxax ,
     &    initialize, firstcell, lastcell,
     & gAMPA_supaxax , gNMDA_supaxax , gGABA_A_supaxax ,
     & Mg, 
     & gapcon_supaxax   ,totSDgj_supaxax    ,gjtable_supaxax , dt,
     &  chi_supaxax,mnaf_supaxax,mnap_supaxax,
     &  hnaf_supaxax,mkdr_supaxax,mka_supaxax,
     &  hka_supaxax,mk2_supaxax,hk2_supaxax,
     &  mkm_supaxax,mkc_supaxax,mkahp_supaxax,
     &  mcat_supaxax,hcat_supaxax,mcal_supaxax,
     &  mar_supaxax)

      ELSE if (nodecell(thisno).eq.'supLTS   ') then
       CALL INTEGRATE_supLTSx  (O, time, num_supLTS  ,
     &    V_supLTS  , curr_supLTS  ,
     &    initialize, firstcell, lastcell,
     & gAMPA_supLTS  , gNMDA_supLTS  , gGABA_A_supLTS  ,
     & Mg, 
     & gapcon_supLTS    ,totSDgj_supLTS     ,gjtable_supLTS  , dt,
     &  chi_supLTS,mnaf_supLTS,mnap_supLTS,
     &  hnaf_supLTS,mkdr_supLTS,mka_supLTS,
     &  hka_supLTS,mk2_supLTS,hk2_supLTS,
     &  mkm_supLTS,mkc_supLTS,mkahp_supLTS,
     &  mcat_supLTS,hcat_supLTS,mcal_supLTS,
     &  mar_supLTS)

      ELSE if (nodecell(thisno).eq.'spinstell') then
       CALL INTEGRATE_spinstelldiegoxB (O, time, num_spinstell,
     &    V_spinstell, curr_spinstell,
     &    initialize, firstcell, lastcell,
     & gAMPA_spinstell, gNMDA_spinstell, gGABA_A_spinstell,
     & gGABA_B_spinstell, Mg, 
     & gapcon_spinstell,totaxgj_spinstell,gjtable_spinstell, dt,
     &  chi_spinstell,mnaf_spinstell,mnap_spinstell,
     &  hnaf_spinstell,mkdr_spinstell,mka_spinstell,
     &  hka_spinstell,mk2_spinstell,hk2_spinstell,
     &  mkm_spinstell,mkc_spinstell,mkahp_spinstell,
     &  mcat_spinstell,hcat_spinstell,mcal_spinstell,
     &  mar_spinstell)

      ELSE if (nodecell(thisno).eq.'tuftIB   ') then
       CALL INTEGRATE_tuftIBVx3B (O, time, num_tuftIB,
     &    V_tuftIB, curr_tuftIB,
     &  initialize, firstcell, lastcell,
     & gAMPA_tuftIB, gNMDA_tuftIB, gGABA_A_tuftIB,
     & gGABA_B_tuftIB, Mg, 
     & gapcon_tuftIB,totaxgj_tuftIB,gjtable_tuftIB, dt,
     & totaxgj_tuft  , gjtable_tuft  , num_tuftRS   ,
c    & vax_tuftRS   ,
     & v_tuftRS   ,
     &  chi_tuftIB,mnaf_tuftIB,mnap_tuftIB,
     &  hnaf_tuftIB,mkdr_tuftIB,mka_tuftIB,
     &  hka_tuftIB,mk2_tuftIB,hk2_tuftIB,
     &  mkm_tuftIB,mkc_tuftIB,mkahp_tuftIB,
     &  mcat_tuftIB,hcat_tuftIB,mcal_tuftIB,
     &  mar_tuftIB,field_1mm       ,field_2mm       ,
     &  scale_tuftIB_gKAHP, scale_tuftIB_gNaP,
     &  scale_tuftIB_gKM  , scale_tuftIB_gKA, 
     &  scale_tuftIB_gCaL, scale_tuftIB_gKC,
     & rel_axonshift_tuftIB,gCaL_tuftIB,Mshift,
     & scale_tuftIB_gAR)

      ELSE if (nodecell(thisno).eq.'tuftRS   ') then
       CALL INTEGRATE_tuftRSXXB (O, time, num_tuftRS,
     &    V_tuftRS, curr_tuftRS,
     & initialize, firstcell, lastcell,
     & gAMPA_tuftRS, gNMDA_tuftRS, gGABA_A_tuftRS,
     & gGABA_B_tuftRS, Mg, 
     & gapcon_tuftRS,totaxgj_tuftRS,gjtable_tuftRS, dt,
     & totaxgj_tuft  , gjtable_tuft  , num_tuftIB   ,
c    & vax_tuftIB   ,
     & v_tuftIB   ,
     &  chi_tuftRS,mnaf_tuftRS,mnap_tuftRS,
     &  hnaf_tuftRS,mkdr_tuftRS,mka_tuftRS,
     &  hka_tuftRS,mk2_tuftRS,hk2_tuftRS,
     &  mkm_tuftRS,mkc_tuftRS,mkahp_tuftRS,
     &  mcat_tuftRS,hcat_tuftRS,mcal_tuftRS,
     &  mar_tuftRS,field_1mm       ,field_2mm       )

      ELSE if (nodecell(thisno).eq.'nontuftRS') then
       CALL INTEGRATE_nontuftRSXXB (O, time, num_nontuftRS,
     &    V_nontuftRS, curr_nontuftRS,
     &  initialize, firstcell, lastcell,
     & gAMPA_nontuftRS, gNMDA_nontuftRS, gGABA_A_nontuftRS,
     & gGABA_B_nontuftRS, Mg, 
     & gapcon_nontuftRS,totaxgj_nontuftRS,gjtable_nontuftRS, dt,
     &  chi_nontuftRS,mnaf_nontuftRS,mnap_nontuftRS,
     &  hnaf_nontuftRS,mkdr_nontuftRS,mka_nontuftRS,
     &  hka_nontuftRS,mk2_nontuftRS,hk2_nontuftRS,
     &  mkm_nontuftRS,mkc_nontuftRS,mkahp_nontuftRS,
     &  mcat_nontuftRS,hcat_nontuftRS,mcal_nontuftRS,
     &  mar_nontuftRS,field_1mm          ,field_2mm          ,
     &  nontuftRS_nFRB)

      ELSE if (nodecell(thisno).eq.'deepbask ') then
       CALL INTEGRATE_deepbaskx  (O, time, num_deepbask ,
     &    V_deepbask , curr_deepbask ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepbask, gNMDA_deepbask, gGABA_A_deepbask,
     & Mg, 
     & gapcon_deepbask  ,totSDgj_deepbask   ,gjtable_deepbask, dt,
     &  chi_deepbask,mnaf_deepbask,mnap_deepbask,
     &  hnaf_deepbask,mkdr_deepbask,mka_deepbask,
     &  hka_deepbask,mk2_deepbask,hk2_deepbask,
     &  mkm_deepbask,mkc_deepbask,mkahp_deepbask,
     &  mcat_deepbask,hcat_deepbask,mcal_deepbask,
     &  mar_deepbask)

      ELSE if (nodecell(thisno).eq.'deepng   ') then
       CALL INTEGRATE_deepng     (O, time, num_deepng   ,
     &    V_deepng   , curr_deepng   ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepng  , gNMDA_deepng  , gGABA_A_deepng  ,
     & Mg, 
     & gapcon_deepng    ,totSDgj_deepng     ,gjtable_deepng  , dt,
     &  chi_deepng  ,mnaf_deepng  ,mnap_deepng  ,
     &  hnaf_deepng  ,mkdr_deepng  ,mka_deepng  ,
     &  hka_deepng  ,mk2_deepng  ,hk2_deepng  ,
     &  mkm_deepng  ,mkc_deepng  ,mkahp_deepng  ,
     &  mcat_deepng  ,hcat_deepng  ,mcal_deepng  ,
     &  mar_deepng  )

      ELSE if (nodecell(thisno).eq.'deepaxax ') then
       CALL INTEGRATE_deepaxaxx (O, time, num_deepaxax ,
     &    V_deepaxax , curr_deepaxax ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepaxax, gNMDA_deepaxax, gGABA_A_deepaxax,
     & Mg, 
     & gapcon_deepaxax  ,totSDgj_deepaxax   ,gjtable_deepaxax, dt,
     &  chi_deepaxax,mnaf_deepaxax,mnap_deepaxax,
     &  hnaf_deepaxax,mkdr_deepaxax,mka_deepaxax,
     &  hka_deepaxax,mk2_deepaxax,hk2_deepaxax,
     &  mkm_deepaxax,mkc_deepaxax,mkahp_deepaxax,
     &  mcat_deepaxax,hcat_deepaxax,mcal_deepaxax,
     &  mar_deepaxax)

      ELSE if (nodecell(thisno).eq.'deepLTS  ') then
       CALL INTEGRATE_deepLTSx  (O, time, num_deepLTS  ,
     &    V_deepLTS  , curr_deepLTS  ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepLTS  , gNMDA_deepLTS  , gGABA_A_deepLTS  ,
     & Mg, 
     & gapcon_deepLTS  ,totSDgj_deepLTS  ,gjtable_deepLTS  , dt,
     &  chi_deepLTS,mnaf_deepLTS,mnap_deepLTS,
     &  hnaf_deepLTS,mkdr_deepLTS,mka_deepLTS,
     &  hka_deepLTS,mk2_deepLTS,hk2_deepLTS,
     &  mkm_deepLTS,mkc_deepLTS,mkahp_deepLTS,
     &  mcat_deepLTS,hcat_deepLTS,mcal_deepLTS,
     &  mar_deepLTS)

      ELSE if (nodecell(thisno).eq.'TCR      ') then
       CALL INTEGRATE_tcrxB      (O, time, num_tcr      ,
     &    V_tcr      , curr_tcr      ,
     & initialize, firstcell, lastcell,
     & gAMPA_tcr      , gNMDA_tcr      , gGABA_A_tcr      ,
     & gGABA_B_tcr, Mg, 
     & gapcon_tcr      ,totaxgj_tcr      ,gjtable_tcr      , dt,
     &  chi_tcr,mnaf_tcr,mnap_tcr,
     &  hnaf_tcr,mkdr_tcr,mka_tcr,
     &  hka_tcr,mk2_tcr,hk2_tcr,
     &  mkm_tcr,mkc_tcr,mkahp_tcr,
     &  mcat_tcr,hcat_tcr,mcal_tcr,
     &  mar_tcr)

      ELSE if (nodecell(thisno).eq.'nRT      ') then
       CALL INTEGRATE_nRTxB      (O, time, num_nRT      ,
     &    V_nRT      , curr_nRT      ,
     & initialize, firstcell, lastcell,
     & gAMPA_nRT      , gNMDA_nRT      , gGABA_A_nRT      ,
     & gGABA_B_nRT, Mg, 
     & gapcon_nRT      ,totSDgj_nRT      ,gjtable_nRT      , dt,
     &  chi_nRT,mnaf_nRT,mnap_nRT,
     &  hnaf_nRT,mkdr_nRT,mka_nRT,
     &  hka_nRT,mk2_nRT,hk2_nRT,
     &  mkm_nRT,mkc_nRT,mkahp_nRT,
     &  mcat_nRT,hcat_nRT,mcal_nRT,
     &  mar_nRT)

      ENDIF
c END INITIALIZATION OF INTEGRATION SUBROUTINES


c BEGIN guts of main program.
c Each node takes care of all the cells of a particular type.
c On a node: enumerate the cells of its type; calculate their
c  synaptic inputs; set applied currents, including those
c  required by ectopic generation; call the numerical integration
c  subroutine; set up the distal_axon vector.  Each node 
c  broadcasts its own distal_axon vector to all the others, and also
c  receives distal_axon vectors from all the others.
c Then, update outtime array and outctr vector.  Repeat.

1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000

c Define shift of tuftIB axonal gNa rate functions, & other axon shifts
       rel_axonshift_tuftIB = 5.0d0 + 0.0d0 * time/timtot
c      rel_axonshift_suppyrRS = 5.d0                      
       rel_axonshift_suppyrRS = 2.d0                      
       rel_axonshift_suppyrFRB= 5.d0                      

! Calculate node-average of total tuftIB GABA-B per cell, & use it to
! control noisepe_tuftIB
        IF (nodecell(thisno) .eq. 'tuftIB   ') THEN
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_tuftIB  
          lastcell = firstcell - 1 + ncellspernode_tuftIB  
           z = 0.d0
           do L = firstcell, lastcell
           do i = 1, numcomp_tuftIB
            z = z + gGABA_B_tuftIB (i,L)
           end do
           end do
             z = z / dble(ncellspernode_tuftIB)
          if (z.le.10.d-3) then
            noisepe_tuftIB = 0.d0 * noisepe_tuftIB_save
          else if (z.le.20.d-3) then
      noisepe_tuftIB=0.d0 *noisepe_tuftIB_save*(1.d0-(z-10.d-3)/10.d-3)
          else
            noisepe_tuftIB = 0.d0
          endif
c          if (mod(O,500).eq.0) then
c            write(6,4908) thisno, noisepe_tuftIB
4908         format('thisno = ',i3,2x,F6.5)
c          endif
        ENDIF   

! Calculate node-average of total tuftRS GABA-B per cell, & use it to
! control noisepe_tuftRS
        IF (nodecell(thisno) .eq. 'tuftRS   ') THEN
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_tuftRS  
          lastcell = firstcell - 1 + ncellspernode_tuftRS  
           z = 0.d0
           do L = firstcell, lastcell
           do i = 1, numcomp_tuftRS
            z = z + gGABA_B_tuftRS (i,L)
           end do
           end do
             z = z / dble(ncellspernode_tuftRS)
          if (z.le.10.d-3) then
            noisepe_tuftRS = 0.d0 * noisepe_tuftRS_save
          else if (z.le.20.d-3) then
      noisepe_tuftRS=0.d0 *noisepe_tuftRS_save*(1.d0-(z-10.d-3)/10.d-3)
          else
            noisepe_tuftRS = 0.d0
          endif
        ENDIF   

       initialize = 1  ! so integration subroutines actually integrate

c Current pulses to selected TCR cells
!      if ((time.gt.300.d0).and.(time.le.310.d0)) then
!       curr_TCR(1,33) = 1.d0
!       curr_TCR(1,42) = 1.d0
!       curr_TCR(1,94) = 1.d0
!       curr_TCR(1,40) = 1.d0
!       curr_TCR(1,29) = 1.d0
!      else
!       curr_TCR(1,33) = 0.d0
!       curr_TCR(1,42) = 0.d0
!       curr_TCR(1,94) = 0.d0
!       curr_TCR(1,40) = 0.d0
!       curr_TCR(1,29) = 0.d0
!      endif

c      IF (THISNO.EQ.0) THEN
       IF (nodecell(thisno) .eq. 'suppyrRS ') THEN
c suppyrRS

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_suppyrRS
          lastcell = firstcell - 1 + ncellspernode_suppyrRS

          IF (MOD(O,how_often).eq.0) then
c 1st set suppyrRS synaptic conductances to 0:

          do i = 1, numcomp_suppyrRS
!         do j = 1, num_suppyrRS
          do j = firstcell, lastcell ! Note
         gAMPA_suppyrRS(i,j)   = 0.d0
         gNMDA_suppyrRS(i,j)   = 0.d0
         gGABA_A_suppyrRS(i,j) = 0.d0
         gGABA_B_suppyrRS(i,j) = 0.d0
          end do
          end do

!        do L = 1, num_suppyrRS
         do L = firstcell, lastcell  ! Note

c Handle suppyrRS   -> suppyrRS
      do i = 1, num_suppyrRS_to_suppyrRS
       j = map_suppyrRS_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrRS(k,L)  = gAMPA_suppyrRS(k,L) +
     &  gAMPA_suppyrRS_to_suppyrRS * delta * z      
! end AMPA part

! NMDA part: 25 Oct 2015, alter code below
! NMDA part: 28 Oct 2015, alter code below
c       if (delta.le.5.d0) then
c       if (delta.le.2.d0) then
        if (delta.le.1.d0) then
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
c    &  gNMDA_suppyrRS_to_suppyrRS * delta * 0.2d0
c    &  gNMDA_suppyrRS_to_suppyrRS * delta * 0.5d0
     &  gNMDA_suppyrRS_to_suppyrRS * delta * 1.0d0
        else
c      dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_suppyrRS
c      dexparg = (delta - 2.d0)/tauNMDA_suppyrRS_to_suppyrRS
       dexparg = (delta - 1.d0)/tauNMDA_suppyrRS_to_suppyrRS
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_suppyrRS_to_suppyrRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_suppyrRS
       if (gNMDA_suppyrRS(k,L).gt.z)
     &  gNMDA_suppyrRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i

c Handle suppyrFRB  -> suppyrRS
      do i = 1, num_suppyrFRB_to_suppyrRS
       j = map_suppyrFRB_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrRS(k,L)  = gAMPA_suppyrRS(k,L) +
     &  gAMPA_suppyrFRB_to_suppyrRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_suppyrFRB_to_suppyrRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_suppyrRS
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_suppyrFRB_to_suppyrRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_suppyrRS
       if (gNMDA_suppyrRS(k,L).gt.z)
     &  gNMDA_suppyrRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supbask    -> suppyrRS
      do i = 1, num_supbask_to_suppyrRS
       j = map_supbask_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_supbask_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_supbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supbask_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrRS(k,L)  = gGABA_A_suppyrRS(k,L) +
     &  gGABA_supbask_to_suppyrRS * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle supng      -> suppyrRS
      do i = 1, num_supng_to_suppyrRS
       j = map_supng_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_supng_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrRS(k,L)  = gGABA_A_suppyrRS(k,L) +
     &  gGABA_supng_to_suppyrRS * z      
! end GABA-A part

      gGABA_B_suppyrRS(k,L) = gGABA_B_suppyrRS(k,L) +
     &   gGABAB_supng_to_suppyrRS * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle supaxax    -> suppyrRS
      do i = 1, num_supaxax_to_suppyrRS
       j = map_supaxax_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_supaxax_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_supaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supaxax_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrRS(k,L)  = gGABA_A_suppyrRS(k,L) +
     &  gGABA_supaxax_to_suppyrRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supLTS     -> suppyrRS
      do i = 1, num_supLTS_to_suppyrRS
       j = map_supLTS_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_supLTS_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrRS(k,L)  = gGABA_A_suppyrRS(k,L) +
     &  gGABA_supLTS_to_suppyrRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle spinstell  -> suppyrRS
      do i = 1, num_spinstell_to_suppyrRS
       j = map_spinstell_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_spinstell_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrRS(k,L)  = gAMPA_suppyrRS(k,L) +
     &  gAMPA_spinstell_to_suppyrRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_spinstell_to_suppyrRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_suppyrRS
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_spinstell_to_suppyrRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_suppyrRS
       if (gNMDA_suppyrRS(k,L).gt.z)
     &  gNMDA_suppyrRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB     -> suppyrRS
      do i = 1, num_tuftIB_to_suppyrRS
       j = map_tuftIB_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrRS(k,L)  = gAMPA_suppyrRS(k,L) +
     &  gAMPA_tuftIB_to_suppyrRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_tuftIB_to_suppyrRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_suppyrRS
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_tuftIB_to_suppyrRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_suppyrRS
       if (gNMDA_suppyrRS(k,L).gt.z)
     &  gNMDA_suppyrRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS     -> suppyrRS
      do i = 1, num_tuftRS_to_suppyrRS
       j = map_tuftRS_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrRS(k,L)  = gAMPA_suppyrRS(k,L) +
     &  gAMPA_tuftRS_to_suppyrRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_tuftRS_to_suppyrRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_suppyrRS
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_tuftRS_to_suppyrRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_suppyrRS
       if (gNMDA_suppyrRS(k,L).gt.z)
     &  gNMDA_suppyrRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepaxax   -> suppyrRS
      do i = 1, num_deepaxax_to_suppyrRS
       j = map_deepaxax_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_deepaxax_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepaxax_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrRS(k,L)  = gGABA_A_suppyrRS(k,L) +
     &  gGABA_deepaxax_to_suppyrRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepLTS    -> suppyrRS
      do i = 1, num_deepLTS_to_suppyrRS
       j = map_deepLTS_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrRS(k,L)  = gGABA_A_suppyrRS(k,L) +
     &  gGABA_deepLTS_to_suppyrRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle TCR        -> suppyrRS
      do i = 1, num_TCR_to_suppyrRS
       j = map_TCR_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_TCR_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrRS(k,L)  = gAMPA_suppyrRS(k,L) +
     &  gAMPA_TCR_to_suppyrRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_TCR_to_suppyrRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_suppyrRS
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_TCR_to_suppyrRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_suppyrRS
       if (gNMDA_suppyrRS(k,L).gt.z)
     &  gNMDA_suppyrRS(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle nontuftRS  -> suppyrRS
      do i = 1, num_nontuftRS_to_suppyrRS
       j = map_nontuftRS_to_suppyrRS(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_suppyrRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_suppyrRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrRS(k,L)  = gAMPA_suppyrRS(k,L) +
     &  gAMPA_nontuftRS_to_suppyrRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_nontuftRS_to_suppyrRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_suppyrRS
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrRS(k,L) = gNMDA_suppyrRS(k,L) +
     &  gNMDA_nontuftRS_to_suppyrRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_suppyrRS
       if (gNMDA_suppyrRS(k,L).gt.z)
     &  gNMDA_suppyrRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i

         end do
c End enumeration of suppyrRS
       ENDIF ! if (mod(O,how_often).eq.0)...

! Define phasic currents to suppyrRS cells, ectopic spikes,
! tonic synaptic conductances

      if (mod(O,200).eq.0) then
       call durand(seed,num_suppyrRS,ranvec_suppyrRS) 
!       do L = 1, num_suppyrRS
        do L = firstcell, lastcell  ! Note
         if ((ranvec_suppyrRS(L).gt.0.d0).and.
     &     (ranvec_suppyrRS(L).le.noisepe_suppyrRS)) then
          curr_suppyrRS(72,L) = 0.4d0
          ectr_suppyrRS = ectr_suppyrRS + 1
         else
          curr_suppyrRS(72,L) = 0.d0
         endif 
        end do
      endif


! Call integration routine for suppyrRS cells
       CALL INTEGRATE_suppyrRSXPB (O, time, num_suppyrRS,
     &    V_suppyrRS, curr_suppyrRS,
     &    initialize, firstcell, lastcell,
     & gAMPA_suppyrRS, gNMDA_suppyrRS, gGABA_A_suppyrRS,
     & gGABA_B_suppyrRS, Mg, 
     & gapcon_suppyrRS  ,totaxgj_suppyrRS   ,gjtable_suppyrRS, dt,
     & totaxgj_suppyr, gjtable_suppyr, num_suppyrFRB,
c    & vax_suppyrFRB,
     & v_suppyrFRB,
     &  chi_suppyrRS,mnaf_suppyrRS,mnap_suppyrRS,
     &  hnaf_suppyrRS,mkdr_suppyrRS,mka_suppyrRS,
     &  hka_suppyrRS,mk2_suppyrRS,hk2_suppyrRS,
     &  mkm_suppyrRS,mkc_suppyrRS,mkahp_suppyrRS,
     &  mcat_suppyrRS,hcat_suppyrRS,mcal_suppyrRS,
     &  mar_suppyrRS,field_1mm ,field_2mm, rel_axonshift_suppyrRS)

       IF (mod(O,5).eq.0) then
! Set up axonal gj voltage array and broadcast it to node 1
! (FRB cells) and receive FRB array - for mixed gj.
!      do L = 1, num_suppyrRS
!       vax_suppyrRS (L) = V_suppyrRS (74,L)
!      end do
!c RE-DO ALL THIS!!
!     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
!    &    mpi_double_precision, 0, mpi_comm_world, info)
!     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
!    &    mpi_double_precision, 1, mpi_comm_world, info)
!     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
!    &    mpi_double_precision, 6, mpi_comm_world, info)
!     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
!    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
            k = table_axgjcompallow_suppyrRS(1)
            gj_axon_suppyrRS (L-firstcell+1) = V_suppyrRS (k,L)
       end do
  
           call mpi_allgather (    gj_axon_suppyrRS,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      ENDIF ! vax set-up, broadcasting and receiving
  

       IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
! also field data                                     
c      do L = 1, num_suppyrRS
       do L = firstcell, lastcell
        distal_axon_suppyrRS (L-firstcell+1) = V_suppyrRS (72,L)
       end do
  
           call mpi_allgather (distal_axon_suppyrRS,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = field_1mm
        field_2mm_local(1) = field_2mm
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        ENDIF ! if (mod(O,how_often).eq.0) ....

          if (mod(O,250).eq.0) then ! write voltages to sysprint
c     write(6,213) time,
c    & ldistal_axon_suppyrRS  (1), 
c    & ldistal_axon_suppyrRS (num_suppyrRS),
c    & ldistal_axon_suppyrFRB (1),
c    & ldistal_axon_suppyrFRB(num_suppyrFRB),
c    & ldistal_axon_supbask   (1),
c    & ldistal_axon_supbask  (num_supbask ),
c    & ldistal_axon_supaxax   (1),
c    & ldistal_axon_supaxax  (num_supaxax  ) 
c     write(6,213) time,
c    & ldistal_axon_supLTS    (1),
c    & ldistal_axon_supLTS   (num_supLTS   ),
c    & ldistal_axon_spinstell (1),
c    & ldistal_axon_spinstell(num_spinstell),
c    & ldistal_axon_tuftIB    (1),
c    & ldistal_axon_tuftIB   (num_tuftIB   ) 
c     write(6,213) time,
c    & ldistal_axon_tuftRS    (1),
c    & ldistal_axon_tuftRS   (num_tuftRS   ),
c    & ldistal_axon_nontuftRS (1),
c    & ldistal_axon_nontuftRS(num_nontuftRS),
c    & ldistal_axon_deepbask  (1),
c    & ldistal_axon_deepbask (num_deepbask ) 
c     write(6,213) time,
c    & ldistal_axon_deepaxax  (1),
c    & ldistal_axon_deepaxax (num_deepaxax ),
c    & ldistal_axon_deepLTS   (1),
c    & ldistal_axon_deepLTS  (num_deepLTS  ),
c    & ldistal_axon_TCR       (1),
c    & ldistal_axon_TCR      (num_TCR      ),
c    & ldistal_axon_nRT       (1),
c    & ldistal_axon_nRT      (num_nRT      ) 
213    FORMAT(f7.2,8f5.0)
          endif ! end writing to sysprint



! END thisno for suppyrRS

c      ELSE IF (THISNO.EQ.1) THEN
       ELSE IF (nodecell(thisno) .eq. 'suppyrFRB') THEN
c suppyrFRB

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_suppyrFRB
          lastcell = firstcell - 1 + ncellspernode_suppyrFRB

         IF (MOD(O,how_often).eq.0) then
c 1st set suppyrFRB synaptic conductances to 0:

          do i = 1, numcomp_suppyrFRB
c         do j = 1, num_suppyrFRB
          do j = firstcell, lastcell
         gAMPA_suppyrFRB(i,j)   = 0.d0 
         gNMDA_suppyrFRB(i,j)   = 0.d0
         gGABA_A_suppyrFRB(i,j) = 0.d0
         gGABA_B_suppyrFRB(i,j) = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS   -> suppyrFRB
      do i = 1, num_suppyrRS_to_suppyrFRB
       j = map_suppyrRS_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrFRB(k,L)  = gAMPA_suppyrFRB(k,L) +
     &  gAMPA_suppyrRS_to_suppyrFRB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_suppyrRS_to_suppyrFRB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_suppyrFRB
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_suppyrRS_to_suppyrFRB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_suppyrFRB
       if (gNMDA_suppyrFRB(k,L).gt.z)
     &  gNMDA_suppyrFRB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB  -> suppyrFRB
      do i = 1, num_suppyrFRB_to_suppyrFRB
       j = map_suppyrFRB_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrFRB(k,L)  = gAMPA_suppyrFRB(k,L) +
     &  gAMPA_suppyrFRB_to_suppyrFRB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_suppyrFRB_to_suppyrFRB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_suppyrFRB
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_suppyrFRB_to_suppyrFRB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_suppyrFRB
       if (gNMDA_suppyrFRB(k,L).gt.z)
     &  gNMDA_suppyrFRB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supbask    -> suppyrFRB
      do i = 1, num_supbask_to_suppyrFRB
       j = map_supbask_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_supbask_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_supbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supbask_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrFRB(k,L)  = gGABA_A_suppyrFRB(k,L) +
     &  gGABA_supbask_to_suppyrFRB * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle supng      -> suppyrFRB
      do i = 1, num_supng_to_suppyrFRB
       j = map_supng_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_supng_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrFRB(k,L)  = gGABA_A_suppyrFRB(k,L) +
     &  gGABA_supng_to_suppyrFRB * z      
! end GABA-A part

      gGABA_B_suppyrFRB(k,L) = gGABA_B_suppyrFRB(k,L) +
     &   gGABAB_supng_to_suppyrFRB * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i



c Handle supaxax    -> suppyrFRB
      do i = 1, num_supaxax_to_suppyrFRB
       j = map_supaxax_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_supaxax_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_supaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supaxax_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrFRB(k,L)  = gGABA_A_suppyrFRB(k,L) +
     &  gGABA_supaxax_to_suppyrFRB * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supLTS     -> suppyrFRB
      do i = 1, num_supLTS_to_suppyrFRB
       j = map_supLTS_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_supLTS_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrFRB(k,L)  = gGABA_A_suppyrFRB(k,L) +
     &  gGABA_supLTS_to_suppyrFRB * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle spinstell  -> suppyrFRB
      do i = 1, num_spinstell_to_suppyrFRB
       j = map_spinstell_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_spinstell_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrFRB(k,L)  = gAMPA_suppyrFRB(k,L) +
     &  gAMPA_spinstell_to_suppyrFRB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_spinstell_to_suppyrFRB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_suppyrFRB
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_spinstell_to_suppyrFRB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_suppyrFRB
       if (gNMDA_suppyrFRB(k,L).gt.z)
     &  gNMDA_suppyrFRB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB     -> suppyrFRB
      do i = 1, num_tuftIB_to_suppyrFRB
       j = map_tuftIB_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrFRB(k,L)  = gAMPA_suppyrFRB(k,L) +
     &  gAMPA_tuftIB_to_suppyrFRB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_tuftIB_to_suppyrFRB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_suppyrFRB
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_tuftIB_to_suppyrFRB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_suppyrFRB
       if (gNMDA_suppyrFRB(k,L).gt.z)
     &  gNMDA_suppyrFRB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS     -> suppyrFRB
      do i = 1, num_tuftRS_to_suppyrFRB
       j = map_tuftRS_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrFRB(k,L)  = gAMPA_suppyrFRB(k,L) +
     &  gAMPA_tuftRS_to_suppyrFRB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_tuftRS_to_suppyrFRB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_suppyrFRB
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_tuftRS_to_suppyrFRB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_suppyrFRB
       if (gNMDA_suppyrFRB(k,L).gt.z)
     &  gNMDA_suppyrFRB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepaxax   -> suppyrFRB
      do i = 1, num_deepaxax_to_suppyrFRB
       j = map_deepaxax_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_deepaxax_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepaxax_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrFRB(k,L)  = gGABA_A_suppyrFRB(k,L) +
     &  gGABA_deepaxax_to_suppyrFRB * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepLTS    -> suppyrFRB
      do i = 1, num_deepLTS_to_suppyrFRB
       j = map_deepLTS_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_suppyrFRB(k,L)  = gGABA_A_suppyrFRB(k,L) +
     &  gGABA_deepLTS_to_suppyrFRB * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle TCR        -> suppyrFRB
      do i = 1, num_TCR_to_suppyrFRB
       j = map_TCR_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_TCR_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrFRB(k,L)  = gAMPA_suppyrFRB(k,L) +
     &  gAMPA_TCR_to_suppyrFRB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_TCR_to_suppyrFRB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_suppyrFRB
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_TCR_to_suppyrFRB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_suppyrFRB
       if (gNMDA_suppyrFRB(k,L).gt.z)
     &  gNMDA_suppyrFRB(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle nontuftRS  -> suppyrFRB
      do i = 1, num_nontuftRS_to_suppyrFRB
       j = map_nontuftRS_to_suppyrFRB(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_suppyrFRB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_suppyrFRB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_suppyrFRB(k,L)  = gAMPA_suppyrFRB(k,L) +
     &  gAMPA_nontuftRS_to_suppyrFRB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_nontuftRS_to_suppyrFRB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_suppyrFRB
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_suppyrFRB(k,L) = gNMDA_suppyrFRB(k,L) +
     &  gNMDA_nontuftRS_to_suppyrFRB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_suppyrFRB
       if (gNMDA_suppyrFRB(k,L).gt.z)
     &  gNMDA_suppyrFRB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of suppyrFRB
         ENDIF ! if (mod(O,how_often).eq.0)...

! Define currents to suppyrFRB cells, ectopic spikes,
! tonic synaptic conductances

      if (mod(O,200).eq.0) then
       call durand(seed,num_suppyrFRB,ranvec_suppyrFRB) 
        do L = firstcell, lastcell
         if ((ranvec_suppyrFRB(L).gt.0.d0).and.
     &     (ranvec_suppyrFRB(L).le.noisepe_suppyrFRB)) then
          curr_suppyrFRB(72,L) = 0.4d0
          ectr_suppyrFRB = ectr_suppyrFRB + 1
         else
          curr_suppyrFRB(72,L) = 0.d0
         endif 
        end do
      endif

! Call integration routine for suppyrFRB cells
       CALL INTEGRATE_suppyrFRBxPB (O, time, num_suppyrFRB,
     &    V_suppyrFRB, curr_suppyrFRB,
     &    initialize, firstcell, lastcell, 
     & gAMPA_suppyrFRB, gNMDA_suppyrFRB, gGABA_A_suppyrFRB,
     & gGABA_B_suppyrFRB, Mg, 
     & gapcon_suppyrFRB ,totaxgj_suppyrFRB  ,gjtable_suppyrFRB, dt,
     & totaxgj_suppyr, gjtable_suppyr, num_suppyrRS,
c    & vax_suppyrRS,
     & v_suppyrRS,
     &  chi_suppyrFRB,mnaf_suppyrFRB,mnap_suppyrFRB,
     &  hnaf_suppyrFRB,mkdr_suppyrFRB,mka_suppyrFRB,
     &  hka_suppyrFRB,mk2_suppyrFRB,hk2_suppyrFRB,
     &  mkm_suppyrFRB,mkc_suppyrFRB,mkahp_suppyrFRB,
     &  mcat_suppyrFRB,hcat_suppyrFRB,mcal_suppyrFRB,
     &  mar_suppyrFRB,field_1mm ,field_2mm, rel_axonshift_suppyrFRB,
     &  scale_suppyrFRB_gNaP)

       IF (mod(O,5).eq.0) then
! Set up axonal gj voltage array and broadcast it to node 0
! (RS cells) and receive RS array - for mixed gj.
!      do L = 1, num_suppyrFRB
!       vax_suppyrFRB (L) = V_suppyrFRB (74,L)
!      end do

!     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
!    &    mpi_double_precision, 0, mpi_comm_world, info)
!     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
!    &    mpi_double_precision, 1, mpi_comm_world, info)
!     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
!    &    mpi_double_precision, 6, mpi_comm_world, info)
!     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
!    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
            k = table_axgjcompallow_suppyrFRB(1)
            gj_axon_suppyrFRB(L-firstcell+1) = V_suppyrFRB(k,L)
       end do
  
           call mpi_allgather (    gj_axon_suppyrFRB,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      ENDIF ! vax set-up, broadcasting and receiving
  

       IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_suppyrFRB
       do L = firstcell, lastcell
        distal_axon_suppyrFRB (L-firstcell+1) = V_suppyrFRB (72,L)
       end do
  
           call mpi_allgather (distal_axon_suppyrFRB,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = field_1mm
        field_2mm_local(1) = field_2mm
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
                        
        ENDIF ! if (mod(O,how_often).eq.0) ....

! END thisno for suppyrFRB

c      ELSE IF (THISNO.EQ.2) THEN
       ELSE IF (nodecell(thisno) .eq. 'supbask  ') THEN
c supbask

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_supbask 
          lastcell = firstcell - 1 + ncellspernode_supbask 

        IF (mod(O,how_often).eq.0) then
c 1st set supbask synaptic conductances to 0:

          do i = 1, numcomp_supbask
          do j = firstcell, lastcell
         gAMPA_supbask(i,j)     = 0.d0
         gNMDA_supbask(i,j)     = 0.d0
         gGABA_A_supbask(i,j)   = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS   -> supbask
      do i = 1, num_suppyrRS_to_supbask  
       j = map_suppyrRS_to_supbask(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_supbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_supbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supbask(k,L)  = gAMPA_supbask(k,L) +
     &  gAMPA_suppyrRS_to_supbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_suppyrRS_to_supbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_supbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_suppyrRS_to_supbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_supbask  
       if (gNMDA_supbask(k,L).gt.z)
     &  gNMDA_supbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB  -> supbask
      do i = 1, num_suppyrFRB_to_supbask  
       j = map_suppyrFRB_to_supbask(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_supbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_supbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supbask(k,L)  = gAMPA_supbask(k,L) +
     &  gAMPA_suppyrFRB_to_supbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_suppyrFRB_to_supbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_supbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_suppyrFRB_to_supbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_supbask  
       if (gNMDA_supbask(k,L).gt.z)
     &  gNMDA_supbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supbask    -> supbask
      do i = 1, num_supbask_to_supbask  
       j = map_supbask_to_supbask(i,L) ! j = presynaptic cell
       k = com_supbask_to_supbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_supbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supbask_to_supbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supbask(k,L)  = gGABA_A_supbask(k,L) +
     &  gGABA_supbask_to_supbask * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle supng      -> supbask 
      do i = 1, num_supng_to_supbask 
       j = map_supng_to_supbask (i,L) ! j = presynaptic cell
       k = com_supng_to_supbask (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_supbask 
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supbask (k,L)  = gGABA_A_supbask (k,L) +
     &  gGABA_supng_to_supbask  * z      
! end GABA-A part

       end do ! m
      end do ! i



c Handle supLTS     -> supbask
      do i = 1, num_supLTS_to_supbask  
       j = map_supLTS_to_supbask(i,L) ! j = presynaptic cell
       k = com_supLTS_to_supbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_supbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supbask(k,L)  = gGABA_A_supbask(k,L) +
     &  gGABA_supLTS_to_supbask * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle spinstell  -> supbask
      do i = 1, num_spinstell_to_supbask  
       j = map_spinstell_to_supbask(i,L) ! j = presynaptic cell
       k = com_spinstell_to_supbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_supbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supbask(k,L)  = gAMPA_supbask(k,L) +
     &  gAMPA_spinstell_to_supbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_spinstell_to_supbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_supbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_spinstell_to_supbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_supbask  
       if (gNMDA_supbask(k,L).gt.z)
     &  gNMDA_supbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB     -> supbask
      do i = 1, num_tuftIB_to_supbask  
       j = map_tuftIB_to_supbask(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_supbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_supbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supbask(k,L)  = gAMPA_supbask(k,L) +
     &  gAMPA_tuftIB_to_supbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_tuftIB_to_supbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_supbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_tuftIB_to_supbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_supbask  
       if (gNMDA_supbask(k,L).gt.z)
     &  gNMDA_supbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS     -> supbask
      do i = 1, num_tuftRS_to_supbask  
       j = map_tuftRS_to_supbask(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_supbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_supbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supbask(k,L)  = gAMPA_supbask(k,L) +
     &  gAMPA_tuftRS_to_supbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_tuftRS_to_supbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_supbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_tuftRS_to_supbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_supbask  
       if (gNMDA_supbask(k,L).gt.z)
     &  gNMDA_supbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepLTS    -> supbask
      do i = 1, num_deepLTS_to_supbask  
       j = map_deepLTS_to_supbask(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_supbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_supbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supbask(k,L)  = gGABA_A_supbask(k,L) +
     &  gGABA_deepLTS_to_supbask * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepTCR    -> supbask
      do i = 1, num_TCR_to_supbask  
       j = map_TCR_to_supbask(i,L) ! j = presynaptic cell
       k = com_TCR_to_supbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_supbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supbask(k,L)  = gAMPA_supbask(k,L) +
     &  gAMPA_TCR_to_supbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_TCR_to_supbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_supbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_TCR_to_supbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_supbask  
       if (gNMDA_supbask(k,L).gt.z)
     &  gNMDA_supbask(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle nontuftRS  -> supbask
      do i = 1, num_nontuftRS_to_supbask  
       j = map_nontuftRS_to_supbask(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_supbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_supbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supbask(k,L)  = gAMPA_supbask(k,L) +
     &  gAMPA_nontuftRS_to_supbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_nontuftRS_to_supbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_supbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supbask(k,L) = gNMDA_supbask(k,L) +
     &  gNMDA_nontuftRS_to_supbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_supbask  
       if (gNMDA_supbask(k,L).gt.z)
     &  gNMDA_supbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of supbask  
        ENDIF  ! if (mod(O,how_often).eq.0) ....

! Define currents to supbask   cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for supbask   cells
       CALL INTEGRATE_supbaskx (O, time, num_supbask ,
     &    V_supbask , curr_supbask ,
     $    initialize, firstcell, lastcell,
     & gAMPA_supbask , gNMDA_supbask , gGABA_A_supbask ,
     & Mg, 
     & gapcon_supbask   ,totSDgj_supbask    ,gjtable_supbask , dt,
     &  chi_supbask,mnaf_supbask,mnap_supbask,
     &  hnaf_supbask,mkdr_supbask,mka_supbask,
     &  hka_supbask,mk2_supbask,hk2_supbask,
     &  mkm_supbask,mkc_supbask,mkahp_supbask,
     &  mcat_supbask,hcat_supbask,mcal_supbask,
     &  mar_supbask)

      if (mod(O,5).eq.0) then

!     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
!    &    mpi_double_precision, 0, mpi_comm_world, info)
!     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
!    &    mpi_double_precision, 1, mpi_comm_world, info)
!     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
!    &    mpi_double_precision, 6, mpi_comm_world, info)
!     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
!    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
c           k = table_axgjcompallow_suppyrRS(1)
            k = numcomp_supbask                   
            gj_axon_supbask  (L-firstcell+1) = V_supbask  (k,L)
       end do
  
           call mpi_allgather (    gj_axon_supbask ,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

      IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_supbask  
       do L = firstcell, lastcell
        distal_axon_supbask   (L-firstcell+1) = V_supbask   (59,L)
       end do
  
           call mpi_allgather (distal_axon_supbask,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = 0.d0     
        field_2mm_local(1) = 0.d0     
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
           ENDIF  ! if (mod(O,how_often).eq.0) ....

! END thisno for supbask

       ELSE IF (nodecell(thisno) .eq. 'supng    ') THEN
c supng  

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_supng   
          lastcell = firstcell - 1 + ncellspernode_supng   

        IF (mod(O,how_often).eq.0) then
c 1st set supng   synaptic conductances to 0:

          do i = 1, numcomp_supbask
          do j = firstcell, lastcell
         gAMPA_supng  (i,j)     = 0.d0
         gNMDA_supng  (i,j)     = 0.d0
         gGABA_A_supng  (i,j)   = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS   -> supng  
      do i = 1, num_suppyrRS_to_supng    
       j = map_suppyrRS_to_supng  (i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_supng  (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_supng    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supng  (k,L)  = gAMPA_supng  (k,L) +
     &  gAMPA_suppyrRS_to_supng   * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supng  (k,L) = gNMDA_supng  (k,L) +
     &  gNMDA_suppyrRS_to_supng   * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_supng    
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supng  (k,L) = gNMDA_supng  (k,L) +
     &  gNMDA_suppyrRS_to_supng   * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_supng    
       if (gNMDA_supng  (k,L).gt.z)
     &  gNMDA_supng  (k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB  -> supng  
      do i = 1, num_suppyrFRB_to_supng    
       j = map_suppyrFRB_to_supng  (i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_supng  (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_supng    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supng  (k,L)  = gAMPA_supng  (k,L) +
     &  gAMPA_suppyrFRB_to_supng   * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supng  (k,L) = gNMDA_supng  (k,L) +
     &  gNMDA_suppyrFRB_to_supng   * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_supng    
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supng  (k,L) = gNMDA_supng  (k,L) +
     &  gNMDA_suppyrFRB_to_supng   * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_supng    
       if (gNMDA_supng  (k,L).gt.z)
     &  gNMDA_supng  (k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supbask    -> supng  
      do i = 1, num_supbask_to_supng    
       j = map_supbask_to_supng  (i,L) ! j = presynaptic cell
       k = com_supbask_to_supng  (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_supbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supbask_to_supng    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supng  (k,L)  = gGABA_A_supng  (k,L) +
     &  gGABA_supbask_to_supng   * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle supng      -> supng   
      do i = 1, num_supng_to_supng   
       j = map_supng_to_supng   (i,L) ! j = presynaptic cell
       k = com_supng_to_supng   (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_supng   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supng   (k,L)  = gGABA_A_supng   (k,L) +
     &  gGABA_supng_to_supng    * z      
! end GABA-A part

       end do ! m
      end do ! i



c Handle deepTCR    -> supng  
      do i = 1, num_TCR_to_supng    
       j = map_TCR_to_supng  (i,L) ! j = presynaptic cell
       k = com_TCR_to_supng  (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_supng    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supng  (k,L)  = gAMPA_supng  (k,L) +
     &  gAMPA_TCR_to_supng   * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supng  (k,L) = gNMDA_supng  (k,L) +
     &  gNMDA_TCR_to_supng   * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_supng    
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supng  (k,L) = gNMDA_supng  (k,L) +
     &  gNMDA_TCR_to_supng   * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_supng    
       if (gNMDA_supng  (k,L).gt.z)
     &  gNMDA_supng  (k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


         end do
c End enumeration of supng    
        ENDIF  ! if (mod(O,how_often).eq.0) ....

! Define currents to supng     cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for supng     cells
       CALL INTEGRATE_supng    (O, time, num_supng   ,
     &    V_supng   , curr_supng   ,
     $    initialize, firstcell, lastcell,
     & gAMPA_supng   , gNMDA_supng   , gGABA_A_supng   ,
     & Mg, 
     & gapcon_supng     ,totSDgj_supng      ,gjtable_supng   , dt,
     &  chi_supng  ,mnaf_supng  ,mnap_supng  ,
     &  hnaf_supng  ,mkdr_supng  ,mka_supng  ,
     &  hka_supng  ,mk2_supng  ,hk2_supng  ,
     &  mkm_supng  ,mkc_supng  ,mkahp_supng  ,
     &  mcat_supng  ,hcat_supng  ,mcal_supng  ,
     &  mar_supng  )

      if (mod(O,5).eq.0) then

       do L = firstcell, lastcell
c           k = table_axgjcompallow_suppyrRS(1)
            k = numcomp_supng                     
            gj_axon_supng    (L-firstcell+1) = V_supng    (k,L)
       end do
  
           call mpi_allgather (    gj_axon_supng   ,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

      IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_supbask  
       do L = firstcell, lastcell
        distal_axon_supng     (L-firstcell+1) = V_supng     (59,L)
       end do
  
           call mpi_allgather (distal_axon_supng  ,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = 0.d0     
        field_2mm_local(1) = 0.d0     
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
           ENDIF  ! if (mod(O,how_often).eq.0) ....

! END thisno for supng  

c      ELSE IF (THISNO.EQ.3) THEN
       ELSE IF (nodecell(thisno) .eq. 'supaxax  ') THEN
c supaxax

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_supaxax 
          lastcell = firstcell - 1 + ncellspernode_supaxax 

         IF (mod(O,how_often).eq.0) then
c 1st set supaxax synaptic conductances to 0:

          do i = 1, numcomp_supaxax
          do j = firstcell, lastcell
         gAMPA_supaxax(i,j)     = 0.d0
         gNMDA_supaxax(i,j)     = 0.d0
         gGABA_A_supaxax(i,j)   = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS   -> supaxax
      do i = 1, num_suppyrRS_to_supaxax  
       j = map_suppyrRS_to_supaxax(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_supaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_supaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supaxax(k,L)  = gAMPA_supaxax(k,L) +
     &  gAMPA_suppyrRS_to_supaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_suppyrRS_to_supaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_supaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_suppyrRS_to_supaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_supaxax  
       if (gNMDA_supaxax(k,L).gt.z)
     &  gNMDA_supaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB  -> supaxax
      do i = 1, num_suppyrFRB_to_supaxax  
       j = map_suppyrFRB_to_supaxax(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_supaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_supaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supaxax(k,L)  = gAMPA_supaxax(k,L) +
     &  gAMPA_suppyrFRB_to_supaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_suppyrFRB_to_supaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_supaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_suppyrFRB_to_supaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_supaxax  
       if (gNMDA_supaxax(k,L).gt.z)
     &  gNMDA_supaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supbask    -> supaxax
      do i = 1, num_supbask_to_supaxax  
       j = map_supbask_to_supaxax(i,L) ! j = presynaptic cell
       k = com_supbask_to_supaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_supbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supbask_to_supaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supaxax(k,L)  = gGABA_A_supaxax(k,L) +
     &  gGABA_supbask_to_supaxax * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supLTS     -> supaxax
      do i = 1, num_supLTS_to_supaxax  
       j = map_supLTS_to_supaxax(i,L) ! j = presynaptic cell
       k = com_supLTS_to_supaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_supaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supaxax(k,L)  = gGABA_A_supaxax(k,L) +
     &  gGABA_supLTS_to_supaxax * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle spinstell  -> supaxax
      do i = 1, num_spinstell_to_supaxax  
       j = map_spinstell_to_supaxax(i,L) ! j = presynaptic cell
       k = com_spinstell_to_supaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_supaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supaxax(k,L)  = gAMPA_supaxax(k,L) +
     &  gAMPA_spinstell_to_supaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_spinstell_to_supaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_supaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_spinstell_to_supaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_supaxax  
       if (gNMDA_supaxax(k,L).gt.z)
     &  gNMDA_supaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB     -> supaxax
      do i = 1, num_tuftIB_to_supaxax  
       j = map_tuftIB_to_supaxax(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_supaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_supaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supaxax(k,L)  = gAMPA_supaxax(k,L) +
     &  gAMPA_tuftIB_to_supaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_tuftIB_to_supaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_supaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_tuftIB_to_supaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_supaxax  
       if (gNMDA_supaxax(k,L).gt.z)
     &  gNMDA_supaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS     -> supaxax
      do i = 1, num_tuftRS_to_supaxax  
       j = map_tuftRS_to_supaxax(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_supaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_supaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supaxax(k,L)  = gAMPA_supaxax(k,L) +
     &  gAMPA_tuftRS_to_supaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_tuftRS_to_supaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_supaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_tuftRS_to_supaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_supaxax  
       if (gNMDA_supaxax(k,L).gt.z)
     &  gNMDA_supaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepLTS    -> supaxax
      do i = 1, num_deepLTS_to_supaxax  
       j = map_deepLTS_to_supaxax(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_supaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_supaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supaxax(k,L)  = gGABA_A_supaxax(k,L) +
     &  gGABA_deepLTS_to_supaxax * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle TCR        -> supaxax
      do i = 1, num_TCR_to_supaxax  
       j = map_TCR_to_supaxax(i,L) ! j = presynaptic cell
       k = com_TCR_to_supaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_supaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supaxax(k,L)  = gAMPA_supaxax(k,L) +
     &  gAMPA_TCR_to_supaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_TCR_to_supaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_supaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_TCR_to_supaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_supaxax  
       if (gNMDA_supaxax(k,L).gt.z)
     &  gNMDA_supaxax(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle nontuftRS  -> supaxax
      do i = 1, num_nontuftRS_to_supaxax  
       j = map_nontuftRS_to_supaxax(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_supaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_supaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supaxax(k,L)  = gAMPA_supaxax(k,L) +
     &  gAMPA_nontuftRS_to_supaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_nontuftRS_to_supaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_supaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supaxax(k,L) = gNMDA_supaxax(k,L) +
     &  gNMDA_nontuftRS_to_supaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_supaxax  
       if (gNMDA_supaxax(k,L).gt.z)
     &  gNMDA_supaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of supaxax  
         ENDIF  ! if (mod(O,how_often).eq.0) ...


! Define currents to supaxax   cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for supaxax   cells
       CALL INTEGRATE_supaxaxx (O, time, num_supaxax ,
     &    V_supaxax , curr_supaxax ,
     &    initialize, firstcell, lastcell,
     & gAMPA_supaxax , gNMDA_supaxax , gGABA_A_supaxax ,
     & Mg, 
     & gapcon_supaxax   ,totSDgj_supaxax    ,gjtable_supaxax , dt,
     &  chi_supaxax,mnaf_supaxax,mnap_supaxax,
     &  hnaf_supaxax,mkdr_supaxax,mka_supaxax,
     &  hka_supaxax,mk2_supaxax,hk2_supaxax,
     &  mkm_supaxax,mkc_supaxax,mkahp_supaxax,
     &  mcat_supaxax,hcat_supaxax,mcal_supaxax,
     &  mar_supaxax)

      if (mod(O,5).eq.0) then

!     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
!    &    mpi_double_precision, 0, mpi_comm_world, info)
!     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
!    &    mpi_double_precision, 1, mpi_comm_world, info)
!     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
!    &    mpi_double_precision, 6, mpi_comm_world, info)
!     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
!    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
c           k = table_axgjcompallow_suppyrRS(1)
            k = numcomp_supaxax                  
            gj_axon_supaxax  (L-firstcell+1) = V_supaxax  (k,L)
       end do
  
           call mpi_allgather (    gj_axon_supaxax ,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_supaxax  
       do L = firstcell, lastcell
        distal_axon_supaxax   (L-firstcell+1) = V_supaxax   (59,L)
       end do
  
           call mpi_allgather (distal_axon_supaxax, 
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = 0.d0     
        field_2mm_local(1) = 0.d0     
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
             ENDIF !  if (mod(O,how_often).eq.0) ...

! END thisno for supaxax

c      ELSE IF (THISNO.EQ.4) THEN
       ELSE IF (nodecell(thisno) .eq. 'supLTS   ') THEN
c supLTS

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_supLTS  
          lastcell = firstcell - 1 + ncellspernode_supLTS  

          IF (mod(O,how_often).eq.0) then
c 1st set supLTS  synaptic conductances to 0:

          do i = 1, numcomp_supLTS
          do j = firstcell, lastcell
         gAMPA_supLTS(i,j)      = 0.d0
         gNMDA_supLTS(i,j)      = 0.d0
         gGABA_A_supLTS(i,j)    = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS   -> supLTS
      do i = 1, num_suppyrRS_to_supLTS   
       j = map_suppyrRS_to_supLTS(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_supLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_supLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supLTS(k,L)  = gAMPA_supLTS(k,L) +
     &  gAMPA_suppyrRS_to_supLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_suppyrRS_to_supLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_supLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_suppyrRS_to_supLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_supLTS  
       if (gNMDA_supLTS(k,L).gt.z)
     &  gNMDA_supLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB  -> supLTS
      do i = 1, num_suppyrFRB_to_supLTS   
       j = map_suppyrFRB_to_supLTS(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_supLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_supLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supLTS(k,L)  = gAMPA_supLTS(k,L) +
     &  gAMPA_suppyrFRB_to_supLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_suppyrFRB_to_supLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_supLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_suppyrFRB_to_supLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_supLTS  
       if (gNMDA_supLTS(k,L).gt.z)
     &  gNMDA_supLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supbask    -> supLTS
      do i = 1, num_supbask_to_supLTS  
       j = map_supbask_to_supLTS(i,L) ! j = presynaptic cell
       k = com_supbask_to_supLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_supbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supbask_to_supLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supLTS(k,L)  = gGABA_A_supLTS(k,L) +
     &  gGABA_supbask_to_supLTS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supLTS     -> supLTS
      do i = 1, num_supLTS_to_supLTS  
       j = map_supLTS_to_supLTS(i,L) ! j = presynaptic cell
       k = com_supLTS_to_supLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_supLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supLTS(k,L)  = gGABA_A_supLTS(k,L) +
     &  gGABA_supLTS_to_supLTS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle spinstell  -> supLTS
      do i = 1, num_spinstell_to_supLTS  
       j = map_spinstell_to_supLTS(i,L) ! j = presynaptic cell
       k = com_spinstell_to_supLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_supLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supLTS(k,L)  = gAMPA_supLTS(k,L) +
     &  gAMPA_spinstell_to_supLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_spinstell_to_supLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_supLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_spinstell_to_supLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_supLTS  
       if (gNMDA_supLTS(k,L).gt.z)
     &  gNMDA_supLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB     -> supLTS
      do i = 1, num_tuftIB_to_supLTS  
       j = map_tuftIB_to_supLTS(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_supLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_supLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supLTS(k,L)  = gAMPA_supLTS(k,L) +
     &  gAMPA_tuftIB_to_supLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_tuftIB_to_supLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_supLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_tuftIB_to_supLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_supLTS  
       if (gNMDA_supLTS(k,L).gt.z)
     &  gNMDA_supLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS     -> supLTS
      do i = 1, num_tuftRS_to_supLTS  
       j = map_tuftRS_to_supLTS(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_supLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_supLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supLTS(k,L)  = gAMPA_supLTS(k,L) +
     &  gAMPA_tuftRS_to_supLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_tuftRS_to_supLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_supLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_tuftRS_to_supLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_supLTS  
       if (gNMDA_supLTS(k,L).gt.z)
     &  gNMDA_supLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepLTS    -> supLTS
      do i = 1, num_deepLTS_to_supLTS   
       j = map_deepLTS_to_supLTS(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_supLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_supLTS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supLTS(k,L)  = gGABA_A_supLTS(k,L) +
     &  gGABA_deepLTS_to_supLTS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle nontuftRS  -> supLTS
      do i = 1, num_nontuftRS_to_supLTS  
       j = map_nontuftRS_to_supLTS(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_supLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_supLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supLTS(k,L)  = gAMPA_supLTS(k,L) +
     &  gAMPA_nontuftRS_to_supLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_nontuftRS_to_supLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_supLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supLTS(k,L) = gNMDA_supLTS(k,L) +
     &  gNMDA_nontuftRS_to_supLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_supLTS  
       if (gNMDA_supLTS(k,L).gt.z)
     &  gNMDA_supLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of supLTS   
        ENDIF  ! if (mod(O,how_often).eq.0) ...

! Define currents to supLTS    cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for supLTS    cells

       CALL INTEGRATE_supLTSx  (O, time, num_supLTS  ,
     &    V_supLTS  , curr_supLTS  ,
     &    initialize, firstcell, lastcell,
     & gAMPA_supLTS  , gNMDA_supLTS  , gGABA_A_supLTS  ,
     & Mg, 
     & gapcon_supLTS    ,totSDgj_supLTS     ,gjtable_supLTS  , dt,
     &  chi_supLTS,mnaf_supLTS,mnap_supLTS,
     &  hnaf_supLTS,mkdr_supLTS,mka_supLTS,
     &  hka_supLTS,mk2_supLTS,hk2_supLTS,
     &  mkm_supLTS,mkc_supLTS,mkahp_supLTS,
     &  mcat_supLTS,hcat_supLTS,mcal_supLTS,
     &  mar_supLTS)

      if (mod(O,5).eq.0) then

c     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
c    &    mpi_double_precision, 0, mpi_comm_world, info)
c     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
c    &    mpi_double_precision, 1, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
c    &    mpi_double_precision, 6, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
c    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
c           k = table_axgjcompallow_suppyrRS(1)
            k = numcomp_supLTS                    
            gj_axon_supLTS   (L-firstcell+1) = V_supLTS   (k,L)
       end do
  
           call mpi_allgather (    gj_axon_supLTS  ,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_supLTS   
       do L = firstcell, lastcell
        distal_axon_supLTS    (L-firstcell+1) = V_supLTS    (59,L)
       end do
  
           call mpi_allgather (distal_axon_supLTS,   
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = 0.d0     
        field_2mm_local(1) = 0.d0     
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
         ENDIF  ! if (mod(O,how_often).eq.0) ...

! END thisno for supLTS
 
c      ELSE IF (THISNO.EQ.5) THEN
       ELSE IF (nodecell(thisno) .eq. 'spinstell') THEN
c spinstell

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_spinstell
          lastcell = firstcell - 1 + ncellspernode_spinstell

       IF (mod(O,how_often).eq.0) then
c 1st set spinstell synaptic conductances to 0:

          do i = 1, numcomp_spinstell
          do j = firstcell, lastcell  
         gAMPA_spinstell(i,j)   = 0.d0
         gNMDA_spinstell(i,j)   = 0.d0
         gGABA_A_spinstell(i,j) = 0.d0
         gGABA_B_spinstell(i,j) = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS    -> spinstell
      do i = 1, num_suppyrRS_to_spinstell
       j = map_suppyrRS_to_spinstell(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_spinstell(k,L)  = gAMPA_spinstell(k,L) +
     &  gAMPA_suppyrRS_to_spinstell * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_suppyrRS_to_spinstell * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_spinstell
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_suppyrRS_to_spinstell * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_spinstell
       if (gNMDA_spinstell(k,L).gt.z)
     &  gNMDA_spinstell(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB   -> spinstell
      do i = 1, num_suppyrFRB_to_spinstell
       j = map_suppyrFRB_to_spinstell(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_spinstell(k,L)  = gAMPA_spinstell(k,L) +
     &  gAMPA_suppyrFRB_to_spinstell * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_suppyrFRB_to_spinstell * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_spinstell
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_suppyrFRB_to_spinstell * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_spinstell
       if (gNMDA_spinstell(k,L).gt.z)
     &  gNMDA_spinstell(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supbask     -> spinstell
      do i = 1, num_supbask_to_spinstell
       j = map_supbask_to_spinstell(i,L) ! j = presynaptic cell
       k = com_supbask_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_supbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supbask_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_spinstell(k,L)  = gGABA_A_spinstell(k,L) +
     &  gGABA_supbask_to_spinstell * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supaxax     -> spinstell
      do i = 1, num_supaxax_to_spinstell
       j = map_supaxax_to_spinstell(i,L) ! j = presynaptic cell
       k = com_supaxax_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_supaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supaxax_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_spinstell(k,L)  = gGABA_A_spinstell(k,L) +
     &  gGABA_supaxax_to_spinstell * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supLTS      -> spinstell
      do i = 1, num_supLTS_to_spinstell
       j = map_supLTS_to_spinstell(i,L) ! j = presynaptic cell
       k = com_supLTS_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_spinstell(k,L)  = gGABA_A_spinstell(k,L) +
     &  gGABA_supLTS_to_spinstell * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle spinstell   -> spinstell
      do i = 1, num_spinstell_to_spinstell
       j = map_spinstell_to_spinstell(i,L) ! j = presynaptic cell
       k = com_spinstell_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_spinstell(k,L)  = gAMPA_spinstell(k,L) +
     &  gAMPA_spinstell_to_spinstell * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_spinstell_to_spinstell * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_spinstell
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_spinstell_to_spinstell * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_spinstell
       if (gNMDA_spinstell(k,L).gt.z)
     &  gNMDA_spinstell(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB      -> spinstell
      do i = 1, num_tuftIB_to_spinstell
       j = map_tuftIB_to_spinstell(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_spinstell(k,L)  = gAMPA_spinstell(k,L) +
     &  gAMPA_tuftIB_to_spinstell * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_tuftIB_to_spinstell * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_spinstell
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_tuftIB_to_spinstell * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_spinstell
       if (gNMDA_spinstell(k,L).gt.z)
     &  gNMDA_spinstell(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS      -> spinstell
      do i = 1, num_tuftRS_to_spinstell
       j = map_tuftRS_to_spinstell(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_spinstell(k,L)  = gAMPA_spinstell(k,L) +
     &  gAMPA_tuftRS_to_spinstell * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_tuftRS_to_spinstell * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_spinstell
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_tuftRS_to_spinstell * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_spinstell
       if (gNMDA_spinstell(k,L).gt.z)
     &  gNMDA_spinstell(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask    -> spinstell
      do i = 1, num_deepbask_to_spinstell
       j = map_deepbask_to_spinstell(i,L) ! j = presynaptic cell
       k = com_deepbask_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_spinstell(k,L)  = gGABA_A_spinstell(k,L) +
     &  gGABA_deepbask_to_spinstell * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle deepng     -> spinstell
      do i = 1, num_deepng_to_spinstell
       j = map_deepng_to_spinstell(i,L) ! j = presynaptic cell
       k = com_deepng_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_deepng_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_spinstell(k,L)  = gGABA_A_spinstell(k,L) +
     &  gGABA_deepng_to_spinstell * z      
! end GABA-A part

      gGABA_B_spinstell(k,L) = gGABA_B_spinstell(k,L) +
     &   gGABAB_deepng_to_spinstell * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i



c Handle deepaxax    -> spinstell
      do i = 1, num_deepaxax_to_spinstell
       j = map_deepaxax_to_spinstell(i,L) ! j = presynaptic cell
       k = com_deepaxax_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepaxax_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_spinstell(k,L)  = gGABA_A_spinstell(k,L) +
     &  gGABA_deepaxax_to_spinstell * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepLTS     -> spinstell
      do i = 1, num_deepLTS_to_spinstell
       j = map_deepLTS_to_spinstell(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_spinstell(k,L)  = gGABA_A_spinstell(k,L) +
     &  gGABA_deepLTS_to_spinstell * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle TCR         -> spinstell
      do i = 1, num_TCR_to_spinstell
       j = map_TCR_to_spinstell(i,L) ! j = presynaptic cell
       k = com_TCR_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_spinstell(k,L)  = gAMPA_spinstell(k,L) +
     &  gAMPA_TCR_to_spinstell * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_TCR_to_spinstell * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_spinstell
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_TCR_to_spinstell * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_spinstell 
       if (gNMDA_spinstell(k,L).gt.z)
     &  gNMDA_spinstell(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle nontuftRS   -> spinstell
      do i = 1, num_nontuftRS_to_spinstell
       j = map_nontuftRS_to_spinstell(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_spinstell(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_spinstell
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_spinstell(k,L)  = gAMPA_spinstell(k,L) +
     &  gAMPA_nontuftRS_to_spinstell * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_nontuftRS_to_spinstell * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_spinstell
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_spinstell(k,L) = gNMDA_spinstell(k,L) +
     &  gNMDA_nontuftRS_to_spinstell * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_spinstell
       if (gNMDA_spinstell(k,L).gt.z)
     &  gNMDA_spinstell(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of spinstell
       ENDIF ! if (mod(O,how_often).eq.0) ...

! Define currents to spinstell cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for spinstell cells
       CALL INTEGRATE_spinstelldiegoxB (O, time, num_spinstell,
     &    V_spinstell, curr_spinstell,
     &    initialize, firstcell, lastcell,
     & gAMPA_spinstell, gNMDA_spinstell, gGABA_A_spinstell,
     & gGABA_B_spinstell, Mg, 
     & gapcon_spinstell,totaxgj_spinstell,gjtable_spinstell, dt,
     &  chi_spinstell,mnaf_spinstell,mnap_spinstell,
     &  hnaf_spinstell,mkdr_spinstell,mka_spinstell,
     &  hka_spinstell,mk2_spinstell,hk2_spinstell,
     &  mkm_spinstell,mkc_spinstell,mkahp_spinstell,
     &  mcat_spinstell,hcat_spinstell,mcal_spinstell,
     &  mar_spinstell)

      if (mod(O,5).eq.0) then

c     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
c    &    mpi_double_precision, 0, mpi_comm_world, info)
c     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
c    &    mpi_double_precision, 1, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
c    &    mpi_double_precision, 6, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
c    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
            k = table_axgjcompallow_spinstell(1)
            gj_axon_spinstell(L-firstcell+1) = V_spinstell(k,L)
       end do
  
           call mpi_allgather (    gj_axon_spinstell,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

       IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_spinstell
       do L = firstcell, lastcell
        distal_axon_spinstell (L-firstcell+1) = V_spinstell (57,L)
       end do
  
           call mpi_allgather (distal_axon_spinstell,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = 0.d0     
        field_2mm_local(1) = 0.d0     
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
           ENDIF !  if (mod(O,how_often).eq.0) ...

! END thisno for spinstell

c      ELSE IF (THISNO.EQ.6) THEN
       ELSE IF (nodecell(thisno) .eq. 'tuftIB   ') THEN
c tuftIB

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_tuftIB  
          lastcell = firstcell - 1 + ncellspernode_tuftIB  

         IF (mod(O,how_often).eq.0) then
c 1st set tuftIB    synaptic conductances to 0:

          do i = 1, numcomp_tuftIB
          do j = firstcell, lastcell
         gAMPA_tuftIB(i,j)      = 0.d0
         gNMDA_tuftIB(i,j)      = 0.d0
         gGABA_A_tuftIB(i,j)    = 0.d0
         gGABA_B_tuftIB(i,j)    = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS    -> tuftIB
      do i = 1, num_suppyrRS_to_tuftIB   
       j = map_suppyrRS_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_tuftIB   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftIB(k,L)  = gAMPA_tuftIB(k,L) +
     &  gAMPA_suppyrRS_to_tuftIB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_suppyrRS_to_tuftIB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_tuftIB   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_suppyrRS_to_tuftIB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_tuftIB
       if (gNMDA_tuftIB(k,L).gt.z)
     &  gNMDA_tuftIB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB   -> tuftIB
      do i = 1, num_suppyrFRB_to_tuftIB   
       j = map_suppyrFRB_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_tuftIB   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftIB(k,L)  = gAMPA_tuftIB(k,L) +
     &  gAMPA_suppyrFRB_to_tuftIB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_suppyrFRB_to_tuftIB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_tuftIB   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_suppyrFRB_to_tuftIB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_tuftIB
       if (gNMDA_tuftIB(k,L).gt.z)
     &  gNMDA_tuftIB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supng      -> tuftIB
      do i = 1, num_supng_to_tuftIB
       j = map_supng_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_supng_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_tuftIB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftIB(k,L)  = gGABA_A_tuftIB(k,L) +
     &  gGABA_supng_to_tuftIB * z      
! end GABA-A part

      gGABA_B_tuftIB(k,L) = gGABA_B_tuftIB(k,L) +
     &   gGABAB_supng_to_tuftIB * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle supaxax     -> tuftIB
      do i = 1, num_supaxax_to_tuftIB   
       j = map_supaxax_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_supaxax_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_supaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supaxax_to_tuftIB   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftIB(k,L)  = gGABA_A_tuftIB(k,L) +
     &  gGABA_supaxax_to_tuftIB * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supLTS      -> tuftIB
      do i = 1, num_supLTS_to_tuftIB   
       j = map_supLTS_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_supLTS_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_tuftIB   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftIB(k,L)  = gGABA_A_tuftIB(k,L) +
     &  gGABA_supLTS_to_tuftIB * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle spinstell   -> tuftIB
      do i = 1, num_spinstell_to_tuftIB  
       j = map_spinstell_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_spinstell_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_tuftIB  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftIB(k,L)  = gAMPA_tuftIB(k,L) +
     &  gAMPA_spinstell_to_tuftIB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_spinstell_to_tuftIB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_tuftIB  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_spinstell_to_tuftIB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_tuftIB  
       if (gNMDA_tuftIB(k,L).gt.z)
     &  gNMDA_tuftIB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB      -> tuftIB
      do i = 1, num_tuftIB_to_tuftIB  
       j = map_tuftIB_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_tuftIB  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftIB(k,L)  = gAMPA_tuftIB(k,L) +
     &  gAMPA_tuftIB_to_tuftIB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_tuftIB_to_tuftIB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_tuftIB  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_tuftIB_to_tuftIB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_tuftIB  
       if (gNMDA_tuftIB(k,L).gt.z)
     &  gNMDA_tuftIB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS      -> tuftIB
      do i = 1, num_tuftRS_to_tuftIB   
       j = map_tuftRS_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_tuftIB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftIB(k,L)  = gAMPA_tuftIB(k,L) +
     &  gAMPA_tuftRS_to_tuftIB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_tuftRS_to_tuftIB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_tuftIB
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_tuftRS_to_tuftIB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_tuftIB   
       if (gNMDA_tuftIB(k,L).gt.z)
     &  gNMDA_tuftIB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask    -> tuftIB
      do i = 1, num_deepbask_to_tuftIB   
       j = map_deepbask_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_deepbask_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_tuftIB   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftIB(k,L)  = gGABA_A_tuftIB(k,L) +
     &  gGABA_deepbask_to_tuftIB * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepng      -> tuftIB
      do i = 1, num_deepng_to_tuftIB
       j = map_deepng_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_deepng_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_deepng_to_tuftIB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftIB(k,L)  = gGABA_A_tuftIB(k,L) +
     &  gGABA_deepng_to_tuftIB * z      
! end GABA-A part

      gGABA_B_tuftIB(k,L) = gGABA_B_tuftIB(k,L) +
     &   gGABAB_deepng_to_tuftIB * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle deepaxax    -> tuftIB
      do i = 1, num_deepaxax_to_tuftIB   
       j = map_deepaxax_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_deepaxax_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepaxax_to_tuftIB   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftIB(k,L)  = gGABA_A_tuftIB(k,L) +
     &  gGABA_deepaxax_to_tuftIB * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepLTS     -> tuftIB
      do i = 1, num_deepLTS_to_tuftIB   
       j = map_deepLTS_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_tuftIB   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftIB(k,L)  = gGABA_A_tuftIB(k,L) +
     &  gGABA_deepLTS_to_tuftIB * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle TCR         -> tuftIB
      do i = 1, num_TCR_to_tuftIB
       j = map_TCR_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_TCR_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_tuftIB
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftIB(k,L)  = gAMPA_tuftIB(k,L) +
     &  gAMPA_TCR_to_tuftIB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_TCR_to_tuftIB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_tuftIB
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_TCR_to_tuftIB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_tuftIB 
       if (gNMDA_tuftIB(k,L).gt.z)
     &  gNMDA_tuftIB(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle nontuftRS   -> tuftIB
      do i = 1, num_nontuftRS_to_tuftIB
       j = map_nontuftRS_to_tuftIB(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_tuftIB(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_tuftIB   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftIB(k,L)  = gAMPA_tuftIB(k,L) +
     &  gAMPA_nontuftRS_to_tuftIB * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_nontuftRS_to_tuftIB * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_tuftIB   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftIB(k,L) = gNMDA_tuftIB(k,L) +
     &  gNMDA_nontuftRS_to_tuftIB * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_tuftIB   
       if (gNMDA_tuftIB(k,L).gt.z)
     &  gNMDA_tuftIB(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of tuftIB   
         ENDIF  ! if (mod(O,how_often).eq.0) ....

! Define currents to tuftIB    cells, ectopic spikes,
! tonic synaptic conductances

      if (mod(O,200).eq.0) then
       call durand(seed,num_tuftIB  ,ranvec_tuftIB  ) 
        do L = firstcell, lastcell
         if ((ranvec_tuftIB  (L).gt.0.d0).and.
     &     (ranvec_tuftIB  (L).le.noisepe_tuftIB  )) then
          curr_tuftIB  (60,L) = 0.4d0
          ectr_tuftIB   = ectr_tuftIB   + 1
         else
          curr_tuftIB  (60,L) = 0.d0
         endif 
        end do
      endif

! Call integration routine for tuftIB    cells
c      CALL INTEGRATE_tuftIB (O, time, num_tuftIB,
       CALL INTEGRATE_tuftIBVx3B (O, time, num_tuftIB,
     &    V_tuftIB, curr_tuftIB,
     &  initialize, firstcell, lastcell,
     & gAMPA_tuftIB, gNMDA_tuftIB, gGABA_A_tuftIB,
     & gGABA_B_tuftIB, Mg, 
     & gapcon_tuftIB,totaxgj_tuftIB,gjtable_tuftIB, dt,
     & totaxgj_tuft  , gjtable_tuft  , num_tuftRS   ,
c    & vax_tuftRS   ,
     & v_tuftRS   ,
     &  chi_tuftIB,mnaf_tuftIB,mnap_tuftIB,
     &  hnaf_tuftIB,mkdr_tuftIB,mka_tuftIB,
     &  hka_tuftIB,mk2_tuftIB,hk2_tuftIB,
     &  mkm_tuftIB,mkc_tuftIB,mkahp_tuftIB,
     &  mcat_tuftIB,hcat_tuftIB,mcal_tuftIB,
     &  mar_tuftIB,field_1mm       ,field_2mm       ,
     &  scale_tuftIB_gKAHP, scale_tuftIB_gNaP,
     &  scale_tuftIB_gKM  , scale_tuftIB_gKA, 
     &  scale_tuftIB_gCaL, scale_tuftIB_gKC,
     & rel_axonshift_tuftIB,gCal_tuftIB,Mshift,
     & scale_tuftIB_gAR)

       IF (mod(O,5).eq.0) then
! Set up axonal gj voltage array and broadcast it to node 7
! (tuftRS cells) and receive tuftRS array - for mixed gj.
c      do L = 1, num_tuftIB  
c       vax_tuftIB (L) = V_tuftIB (61,L)
c      end do

c     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
c    &    mpi_double_precision, 0, mpi_comm_world, info)
c     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
c    &    mpi_double_precision, 1, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
c    &    mpi_double_precision, 6, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
c    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
            k = table_axgjcompallow_tuftIB  (1)
            gj_axon_tuftIB   (L-firstcell+1) = V_tuftIB   (k,L)
       end do
  
           call mpi_allgather (    gj_axon_tuftIB  ,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      ENDIF ! vax set-up, broadcasting and receiving
  

        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_tuftIB   
       do L = firstcell, lastcell
        distal_axon_tuftIB    (L-firstcell+1) = V_tuftIB    (60,L)
       end do
  
           call mpi_allgather (distal_axon_tuftIB,  
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = field_1mm
        field_2mm_local(1) = field_2mm
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
           ENDIF  ! if (mod(O,how_often).eq.0) ...

! END thisno for tuftIB

c      ELSE IF (THISNO.EQ.7) THEN
       ELSE IF (nodecell(thisno) .eq. 'tuftRS   ') THEN
c tuftRS

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_tuftRS  
          lastcell = firstcell - 1 + ncellspernode_tuftRS  

         IF (mod(O,how_often).eq.0) then
c 1st set tuftRS    synaptic conductances to 0:

          do i = 1, numcomp_tuftRS
          do j = firstcell, lastcell
         gAMPA_tuftRS(i,j)      = 0.d0 
         gNMDA_tuftRS(i,j)      = 0.d0
         gGABA_A_tuftRS(i,j)    = 0.d0
         gGABA_B_tuftRS(i,j)    = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS    -> tuftRS
      do i = 1, num_suppyrRS_to_tuftRS   
       j = map_suppyrRS_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_tuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftRS(k,L)  = gAMPA_tuftRS(k,L) +
     &  gAMPA_suppyrRS_to_tuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_suppyrRS_to_tuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_tuftRS   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_suppyrRS_to_tuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_tuftRS
       if (gNMDA_tuftRS(k,L).gt.z)
     &  gNMDA_tuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB   -> tuftRS 
      do i = 1, num_suppyrFRB_to_tuftRS   
       j = map_suppyrFRB_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_tuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftRS(k,L)  = gAMPA_tuftRS(k,L) +
     &  gAMPA_suppyrFRB_to_tuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_suppyrFRB_to_tuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_tuftRS   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_suppyrFRB_to_tuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_tuftRS
       if (gNMDA_tuftRS(k,L).gt.z)
     &  gNMDA_tuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supng      -> tuftRS
      do i = 1, num_supng_to_tuftRS
       j = map_supng_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_supng_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_tuftRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftRS(k,L)  = gGABA_A_tuftRS(k,L) +
     &  gGABA_supng_to_tuftRS * z      
! end GABA-A part

      gGABA_B_tuftRS(k,L) = gGABA_B_tuftRS(k,L) +
     &   gGABAB_supng_to_tuftRS * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle supaxax     -> tuftRS
      do i = 1, num_supaxax_to_tuftRS   
       j = map_supaxax_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_supaxax_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_supaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supaxax_to_tuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftRS(k,L)  = gGABA_A_tuftRS(k,L) +
     &  gGABA_supaxax_to_tuftRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supLTS      -> tuftRS
      do i = 1, num_supLTS_to_tuftRS   
       j = map_supLTS_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_supLTS_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_tuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftRS(k,L)  = gGABA_A_tuftRS(k,L) +
     &  gGABA_supLTS_to_tuftRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle spinstell   -> tuftRS
      do i = 1, num_spinstell_to_tuftRS  
       j = map_spinstell_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_spinstell_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_tuftRS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftRS(k,L)  = gAMPA_tuftRS(k,L) +
     &  gAMPA_spinstell_to_tuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_spinstell_to_tuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_tuftRS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_spinstell_to_tuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_tuftRS  
       if (gNMDA_tuftRS(k,L).gt.z)
     &  gNMDA_tuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB      -> tuftRS
      do i = 1, num_tuftIB_to_tuftRS  
       j = map_tuftIB_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_tuftRS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftRS(k,L)  = gAMPA_tuftRS(k,L) +
     &  gAMPA_tuftIB_to_tuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_tuftIB_to_tuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_tuftRS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_tuftIB_to_tuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_tuftRS  
       if (gNMDA_tuftRS(k,L).gt.z)
     &  gNMDA_tuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS      -> tuftRS
      do i = 1, num_tuftRS_to_tuftRS  
       j = map_tuftRS_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_tuftRS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftRS(k,L)  = gAMPA_tuftRS(k,L) +
     &  gAMPA_tuftRS_to_tuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_tuftRS_to_tuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_tuftRS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_tuftRS_to_tuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_tuftRS  
       if (gNMDA_tuftRS(k,L).gt.z)
     &  gNMDA_tuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask    -> tuftRS
      do i = 1, num_deepbask_to_tuftRS   
       j = map_deepbask_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_deepbask_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_tuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftRS(k,L)  = gGABA_A_tuftRS(k,L) +
     &  gGABA_deepbask_to_tuftRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepng      -> tuftRS
      do i = 1, num_deepng_to_tuftRS
       j = map_deepng_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_deepng_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_deepng_to_tuftRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftRS(k,L)  = gGABA_A_tuftRS(k,L) +
     &  gGABA_deepng_to_tuftRS * z      
! end GABA-A part

      gGABA_B_tuftRS(k,L) = gGABA_B_tuftRS(k,L) +
     &   gGABAB_deepng_to_tuftRS * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle deepaxax    -> tuftRS
      do i = 1, num_deepaxax_to_tuftRS   
       j = map_deepaxax_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_deepaxax_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepaxax_to_tuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftRS(k,L)  = gGABA_A_tuftRS(k,L) +
     &  gGABA_deepaxax_to_tuftRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepLTS     -> tuftRS
      do i = 1, num_deepLTS_to_tuftRS   
       j = map_deepLTS_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_tuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_tuftRS(k,L)  = gGABA_A_tuftRS(k,L) +
     &  gGABA_deepLTS_to_tuftRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle TCR         -> tuftRS
      do i = 1, num_TCR_to_tuftRS
       j = map_TCR_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_TCR_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_tuftRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftRS(k,L)  = gAMPA_tuftRS(k,L) +
     &  gAMPA_TCR_to_tuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_TCR_to_tuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_tuftRS
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_TCR_to_tuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_tuftRS 
       if (gNMDA_tuftRS(k,L).gt.z)
     &  gNMDA_tuftRS(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle nontuftRS   -> tuftRS
      do i = 1, num_nontuftRS_to_tuftRS  
       j = map_nontuftRS_to_tuftRS(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_tuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_tuftRS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_tuftRS(k,L)  = gAMPA_tuftRS(k,L) +
     &  gAMPA_nontuftRS_to_tuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_nontuftRS_to_tuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_tuftRS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_tuftRS(k,L) = gNMDA_tuftRS(k,L) +
     &  gNMDA_nontuftRS_to_tuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_tuftRS  
       if (gNMDA_tuftRS(k,L).gt.z)
     &  gNMDA_tuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of tuftRS   
        ENDIF  ! if (mod(O,how_often).eq.0) ...

! Define currents to tuftRS    cells, ectopic spikes,
! tonic synaptic conductances

      if (mod(O,200).eq.0) then
       call durand(seed,num_tuftRS  ,ranvec_tuftRS  ) 
        do L = firstcell, lastcell
         if ((ranvec_tuftRS  (L).gt.0.d0).and.
     &     (ranvec_tuftRS  (L).le.noisepe_tuftRS  )) then
          curr_tuftRS  (60,L) = 0.4d0
          ectr_tuftRS   = ectr_tuftRS   + 1
         else
          curr_tuftRS  (60,L) = 0.d0
         endif 
        end do
      endif

! Call integration routine for tuftRS    cells
       CALL INTEGRATE_tuftRSXXB (O, time, num_tuftRS,
     &    V_tuftRS, curr_tuftRS,
     & initialize, firstcell, lastcell,
     & gAMPA_tuftRS, gNMDA_tuftRS, gGABA_A_tuftRS,
     & gGABA_B_tuftRS, Mg, 
     & gapcon_tuftRS,totaxgj_tuftRS,gjtable_tuftRS, dt,
     & totaxgj_tuft  , gjtable_tuft  , num_tuftIB   ,
c    & vax_tuftIB   ,
     & v_tuftIB   ,
     &  chi_tuftRS,mnaf_tuftRS,mnap_tuftRS,
     &  hnaf_tuftRS,mkdr_tuftRS,mka_tuftRS,
     &  hka_tuftRS,mk2_tuftRS,hk2_tuftRS,
     &  mkm_tuftRS,mkc_tuftRS,mkahp_tuftRS,
     &  mcat_tuftRS,hcat_tuftRS,mcal_tuftRS,
     &  mar_tuftRS,field_1mm       ,field_2mm       )

       IF (mod(O,5).eq.0) then
! Set up axonal gj voltage array and broadcast it to node 6
! (tuftIB cells) and receive tuftIB array - for mixed gj.
c      do L = 1, num_tuftRS  
c       vax_tuftRS (L) = V_tuftRS (61,L)
c      end do

c     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
c    &    mpi_double_precision, 0, mpi_comm_world, info)
c     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
c    &    mpi_double_precision, 1, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
c    &    mpi_double_precision, 6, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
c    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
            k = table_axgjcompallow_tuftRS  (1)
            gj_axon_tuftRS   (L-firstcell+1) = V_tuftRS   (k,L)
       end do
  
           call mpi_allgather (    gj_axon_tuftRS  ,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      ENDIF ! vax set-up, broadcasting and receiving
  

       IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_tuftRS   
       do L = firstcell, lastcell
        distal_axon_tuftRS    (L-firstcell+1) = V_tuftRS    (60,L)
       end do
  
           call mpi_allgather (distal_axon_tuftRS,  
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = field_1mm
        field_2mm_local(1) = field_2mm
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
         ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for tuftRS

c      ELSE IF (THISNO.EQ.8) THEN
       ELSE IF (nodecell(thisno) .eq. 'nontuftRS') THEN
c nontuftRS

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_nontuftRS
          lastcell = firstcell - 1 + ncellspernode_nontuftRS

         IF (mod(O,how_often).eq.0) then
c 1st set nontuftRS synaptic conductances to 0:

          do i = 1, numcomp_nontuftRS
          do j = firstcell, lastcell
         gAMPA_nontuftRS(i,j)   = 0.d0 
         gNMDA_nontuftRS(i,j)   = 0.d0 
         gGABA_A_nontuftRS(i,j) = 0.d0
         gGABA_B_nontuftRS(i,j) = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS   -> nontuftRS
      do i = 1, num_suppyrRS_to_nontuftRS   
       j = map_suppyrRS_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_nontuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_nontuftRS(k,L)  = gAMPA_nontuftRS(k,L) +
     &  gAMPA_suppyrRS_to_nontuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_suppyrRS_to_nontuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_nontuftRS   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_suppyrRS_to_nontuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_nontuftRS
       if (gNMDA_nontuftRS(k,L).gt.z)
     &  gNMDA_nontuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB  -> nontuftRS
      do i = 1, num_suppyrFRB_to_nontuftRS   
       j = map_suppyrFRB_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_nontuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_nontuftRS(k,L)  = gAMPA_nontuftRS(k,L) +
     &  gAMPA_suppyrFRB_to_nontuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_suppyrFRB_to_nontuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_nontuftRS   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_suppyrFRB_to_nontuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_nontuftRS
       if (gNMDA_nontuftRS(k,L).gt.z)
     &  gNMDA_nontuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supng      -> nontuftRS
      do i = 1, num_supng_to_nontuftRS
       j = map_supng_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_supng_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_nontuftRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_nontuftRS(k,L)  = gGABA_A_nontuftRS(k,L) +
     &  gGABA_supng_to_nontuftRS * z      
! end GABA-A part

      gGABA_B_nontuftRS(k,L) = gGABA_B_nontuftRS(k,L) +
     &   gGABAB_supng_to_nontuftRS * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle supaxax    -> nontuftRS
      do i = 1, num_supaxax_to_nontuftRS   
       j = map_supaxax_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_supaxax_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_supaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supaxax_to_nontuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_nontuftRS(k,L)  = gGABA_A_nontuftRS(k,L) +
     &  gGABA_supaxax_to_nontuftRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supLTS     -> nontuftRS
      do i = 1, num_supLTS_to_nontuftRS   
       j = map_supLTS_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_supLTS_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_nontuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_nontuftRS(k,L)  = gGABA_A_nontuftRS(k,L) +
     &  gGABA_supLTS_to_nontuftRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle spinstell  -> nontuftRS
      do i = 1, num_spinstell_to_nontuftRS  
       j = map_spinstell_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_spinstell_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_nontuftRS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_nontuftRS(k,L)  = gAMPA_nontuftRS(k,L) +
     &  gAMPA_spinstell_to_nontuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_spinstell_to_nontuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_nontuftRS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_spinstell_to_nontuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_nontuftRS  
       if (gNMDA_nontuftRS(k,L).gt.z)
     &  gNMDA_nontuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB     -> nontuftRS
      do i = 1, num_tuftIB_to_nontuftRS  
       j = map_tuftIB_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_nontuftRS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_nontuftRS(k,L)  = gAMPA_nontuftRS(k,L) +
     &  gAMPA_tuftIB_to_nontuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_tuftIB_to_nontuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_nontuftRS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_tuftIB_to_nontuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_nontuftRS  
       if (gNMDA_nontuftRS(k,L).gt.z)
     &  gNMDA_nontuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS     -> nontuftRS
      do i = 1, num_tuftRS_to_nontuftRS  
       j = map_tuftRS_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_nontuftRS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_nontuftRS(k,L)  = gAMPA_nontuftRS(k,L) +
     &  gAMPA_tuftRS_to_nontuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_tuftRS_to_nontuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_nontuftRS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_tuftRS_to_nontuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_nontuftRS  
       if (gNMDA_nontuftRS(k,L).gt.z)
     &  gNMDA_nontuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask   -> nontuftRS
      do i = 1, num_deepbask_to_nontuftRS   
       j = map_deepbask_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_deepbask_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_nontuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_nontuftRS(k,L)  = gGABA_A_nontuftRS(k,L) +
     &  gGABA_deepbask_to_nontuftRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepng      -> nontuftRS
      do i = 1, num_deepng_to_nontuftRS
       j = map_deepng_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_deepng_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_deepng_to_nontuftRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_nontuftRS(k,L)  = gGABA_A_nontuftRS(k,L) +
     &  gGABA_deepng_to_nontuftRS * z      
! end GABA-A part

      gGABA_B_nontuftRS(k,L) = gGABA_B_nontuftRS(k,L) +
     &   gGABAB_deepng_to_nontuftRS * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle deepaxax   -> nontuftRS
      do i = 1, num_deepaxax_to_nontuftRS   
       j = map_deepaxax_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_deepaxax_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepaxax(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepaxax(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepaxax_to_nontuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_nontuftRS(k,L)  = gGABA_A_nontuftRS(k,L) +
     &  gGABA_deepaxax_to_nontuftRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepLTS    -> nontuftRS
      do i = 1, num_deepLTS_to_nontuftRS   
       j = map_deepLTS_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_nontuftRS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_nontuftRS(k,L)  = gGABA_A_nontuftRS(k,L) +
     &  gGABA_deepLTS_to_nontuftRS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle TCR        -> nontuftRS
      do i = 1, num_TCR_to_nontuftRS
       j = map_TCR_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_TCR_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_nontuftRS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_nontuftRS(k,L)  = gAMPA_nontuftRS(k,L) +
     &  gAMPA_TCR_to_nontuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_TCR_to_nontuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_nontuftRS
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_TCR_to_nontuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_nontuftRS 
       if (gNMDA_nontuftRS(k,L).gt.z)
     &  gNMDA_nontuftRS(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle nontuftRS  -> nontuftRS
      do i = 1, num_nontuftRS_to_nontuftRS  
       j = map_nontuftRS_to_nontuftRS(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_nontuftRS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_nontuftRS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_nontuftRS(k,L)  = gAMPA_nontuftRS(k,L) +
     &  gAMPA_nontuftRS_to_nontuftRS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_nontuftRS_to_nontuftRS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_nontuftRS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_nontuftRS(k,L) = gNMDA_nontuftRS(k,L) +
     &  gNMDA_nontuftRS_to_nontuftRS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_nontuftRS  
       if (gNMDA_nontuftRS(k,L).gt.z)
     &  gNMDA_nontuftRS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of nontuftRS   
          ENDIF  ! if (mod(O,how_often).eq.0) ...

! Define currents to nontuftRS    cells, ectopic spikes,
! tonic synaptic conductances

      if (mod(O,200).eq.0) then
       call durand(seed,num_nontuftRS  ,ranvec_nontuftRS  ) 
        do L = firstcell, lastcell
         if ((ranvec_nontuftRS  (L).gt.0.d0).and.
     &     (ranvec_nontuftRS  (L).le.noisepe_nontuftRS  )) then
          curr_nontuftRS  (48,L) = 0.4d0
          ectr_nontuftRS   = ectr_nontuftRS   + 1
         else
          curr_nontuftRS  (48,L) = 0.d0
         endif 
        end do
      endif

! Call integration routine for nontuftRS    cells
       CALL INTEGRATE_nontuftRSXXB (O, time, num_nontuftRS,
     &    V_nontuftRS, curr_nontuftRS,
     &  initialize, firstcell, lastcell,
     & gAMPA_nontuftRS, gNMDA_nontuftRS, gGABA_A_nontuftRS,
     & gGABA_B_nontuftRS, Mg, 
     & gapcon_nontuftRS,totaxgj_nontuftRS,gjtable_nontuftRS, dt,
     &  chi_nontuftRS,mnaf_nontuftRS,mnap_nontuftRS,
     &  hnaf_nontuftRS,mkdr_nontuftRS,mka_nontuftRS,
     &  hka_nontuftRS,mk2_nontuftRS,hk2_nontuftRS,
     &  mkm_nontuftRS,mkc_nontuftRS,mkahp_nontuftRS,
     &  mcat_nontuftRS,hcat_nontuftRS,mcal_nontuftRS,
     &  mar_nontuftRS,field_1mm          ,field_2mm          ,
     &  nontuftRS_nFRB)

      if (mod(O,5).eq.0) then

c     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
c    &    mpi_double_precision, 0, mpi_comm_world, info)
c     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
c    &    mpi_double_precision, 1, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
c    &    mpi_double_precision, 6, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
c    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
            k = table_axgjcompallow_nontuftRS(1)
            gj_axon_nontuftRS(L-firstcell+1) = V_nontuftRS(k,L)
       end do
  
           call mpi_allgather (    gj_axon_nontuftRS,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_nontuftRS   
       do L = firstcell, lastcell
        distal_axon_nontuftRS    (L-firstcell+1) = V_nontuftRS    (48,L)
       end do
  
           call mpi_allgather (distal_axon_nontuftRS,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = field_1mm
        field_2mm_local(1) = field_2mm
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
         ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for nontuftRS

c      ELSE IF (THISNO.EQ.9) THEN
       ELSE IF (nodecell(thisno) .eq. 'deepbask ') THEN
c deepbask

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_deepbask
          lastcell = firstcell - 1 + ncellspernode_deepbask

          IF (mod(O,how_often).eq.0) then
c 1st set deepbask  synaptic conductances to 0:

          do i = 1, numcomp_deepbask
          do j = firstcell, lastcell
         gAMPA_deepbask(i,j)    = 0.d0
         gNMDA_deepbask(i,j)    = 0.d0
         gGABA_A_deepbask(i,j)  = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS   -> deepbask
      do i = 1, num_suppyrRS_to_deepbask  
       j = map_suppyrRS_to_deepbask(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_deepbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_suppyrRS_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_suppyrRS_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_deepbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_suppyrRS_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_deepbask  
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB  -> deepbask
      do i = 1, num_suppyrFRB_to_deepbask  
       j = map_suppyrFRB_to_deepbask(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_deepbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_suppyrFRB_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_suppyrFRB_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_deepbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_suppyrFRB_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_deepbask  
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supLTS     -> deepbask
      do i = 1, num_supLTS_to_deepbask    
       j = map_supLTS_to_deepbask(i,L) ! j = presynaptic cell
       k = com_supLTS_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_deepbask    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepbask(k,L)  = gGABA_A_deepbask(k,L) +
     &  gGABA_supLTS_to_deepbask * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle spinstell  -> deepbask
      do i = 1, num_spinstell_to_deepbask   
       j = map_spinstell_to_deepbask(i,L) ! j = presynaptic cell
       k = com_spinstell_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_deepbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_spinstell_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_spinstell_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_deepbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_spinstell_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_deepbask  
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB     -> deepbask
      do i = 1, num_tuftIB_to_deepbask   
       j = map_tuftIB_to_deepbask(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_deepbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_tuftIB_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_tuftIB_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_deepbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_tuftIB_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_deepbask  
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS     -> deepbask
      do i = 1, num_tuftRS_to_deepbask   
       j = map_tuftRS_to_deepbask(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_deepbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_tuftRS_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_tuftRS_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_deepbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_tuftRS_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_deepbask  
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask   -> deepbask
      do i = 1, num_deepbask_to_deepbask    
       j = map_deepbask_to_deepbask(i,L) ! j = presynaptic cell
       k = com_deepbask_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_deepbask    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepbask(k,L)  = gGABA_A_deepbask(k,L) +
     &  gGABA_deepbask_to_deepbask * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle deepng     -> deepbask
      do i = 1, num_deepng_to_deepbask    
       j = map_deepng_to_deepbask(i,L) ! j = presynaptic cell
       k = com_deepng_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepng_to_deepbask    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepbask(k,L)  = gGABA_A_deepbask(k,L) +
     &  gGABA_deepng_to_deepbask * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepLTS    -> deepbask
      do i = 1, num_deepLTS_to_deepbask    
       j = map_deepLTS_to_deepbask(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_deepbask    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepbask(k,L)  = gGABA_A_deepbask(k,L) +
     &  gGABA_deepLTS_to_deepbask * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle TCR        -> deepbask
      do i = 1, num_TCR_to_deepbask 
       j = map_TCR_to_deepbask(i,L) ! j = presynaptic cell
       k = com_TCR_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_deepbask 
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_TCR_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_TCR_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_deepbask 
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_TCR_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_deepbask  
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle nontuftRS  -> deepbask
      do i = 1, num_nontuftRS_to_deepbask
       j = map_nontuftRS_to_deepbask(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_deepbask
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_nontuftRS_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_nontuftRS_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_deepbask
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_nontuftRS_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_deepbask
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of deepbask    
         ENDIF ! if (mod(O,how_often).eq.0) ...

! Define currents to deepbask     cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for deepbask     cells
       CALL INTEGRATE_deepbaskx  (O, time, num_deepbask ,
     &    V_deepbask , curr_deepbask ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepbask, gNMDA_deepbask, gGABA_A_deepbask,
     & Mg, 
     & gapcon_deepbask  ,totSDgj_deepbask   ,gjtable_deepbask, dt,
     &  chi_deepbask,mnaf_deepbask,mnap_deepbask,
     &  hnaf_deepbask,mkdr_deepbask,mka_deepbask,
     &  hka_deepbask,mk2_deepbask,hk2_deepbask,
     &  mkm_deepbask,mkc_deepbask,mkahp_deepbask,
     &  mcat_deepbask,hcat_deepbask,mcal_deepbask,
     &  mar_deepbask)

      if (mod(O,5).eq.0) then

c     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
c    &    mpi_double_precision, 0, mpi_comm_world, info)
c     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
c    &    mpi_double_precision, 1, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
c    &    mpi_double_precision, 6, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
c    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
c           k = table_axgjcompallow_suppyrRS(1)
            k = numcomp_deepbask                  
            gj_axon_deepbask (L-firstcell+1) = V_deepbask (k,L)
       end do
  
           call mpi_allgather (    gj_axon_deepbask,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_deepbask    
       do L = firstcell, lastcell
        distal_axon_deepbask     (L-firstcell+1) = V_deepbask     (59,L)
       end do
  
           call mpi_allgather (distal_axon_deepbask,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = 0.d0     
        field_2mm_local(1) = 0.d0     
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
           ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for deepbask

       ELSE IF (nodecell(thisno) .eq. 'deepng   ') THEN
c deepng  

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_deepng  
          lastcell = firstcell - 1 + ncellspernode_deepng  

          IF (mod(O,how_often).eq.0) then
c 1st set deepng    synaptic conductances to 0:

          do i = 1, numcomp_deepng  
          do j = firstcell, lastcell
         gAMPA_deepng  (i,j)    = 0.d0
         gNMDA_deepng  (i,j)    = 0.d0
         gGABA_A_deepng  (i,j)  = 0.d0
          end do
          end do

         do L = firstcell, lastcell

c Handle spinstell  -> deepng
      do i = 1, num_spinstell_to_deepng   
       j = map_spinstell_to_deepng(i,L) ! j = presynaptic cell
       k = com_spinstell_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_deepng  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepng(k,L)  = gAMPA_deepng(k,L) +
     &  gAMPA_spinstell_to_deepng * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_spinstell_to_deepng * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_deepng  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_spinstell_to_deepng * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_deepng  
       if (gNMDA_deepng(k,L).gt.z)
     &  gNMDA_deepng(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB     -> deepng
      do i = 1, num_tuftIB_to_deepng   
       j = map_tuftIB_to_deepng(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_deepng  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepng(k,L)  = gAMPA_deepng(k,L) +
     &  gAMPA_tuftIB_to_deepng * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_tuftIB_to_deepng * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_deepng  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_tuftIB_to_deepng * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_deepng  
       if (gNMDA_deepng(k,L).gt.z)
     &  gNMDA_deepng(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS     -> deepng
      do i = 1, num_tuftRS_to_deepng   
       j = map_tuftRS_to_deepng(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_deepng  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepng(k,L)  = gAMPA_deepng(k,L) +
     &  gAMPA_tuftRS_to_deepng * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_tuftRS_to_deepng * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_deepng  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_tuftRS_to_deepng * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_deepng  
       if (gNMDA_deepng(k,L).gt.z)
     &  gNMDA_deepng(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask   -> deepng  
      do i = 1, num_deepbask_to_deepng      
       j = map_deepbask_to_deepng  (i,L) ! j = presynaptic cell
       k = com_deepbask_to_deepng  (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_deepng      
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepng  (k,L)  = gGABA_A_deepng  (k,L) +
     &  gGABA_deepbask_to_deepng   * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle deepng     -> deepng
      do i = 1, num_deepng_to_deepng    
       j = map_deepng_to_deepng(i,L) ! j = presynaptic cell
       k = com_deepng_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepng_to_deepng    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepng(k,L)  = gGABA_A_deepng(k,L) +
     &  gGABA_deepng_to_deepng * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle TCR        -> deepng
      do i = 1, num_TCR_to_deepng 
       j = map_TCR_to_deepng(i,L) ! j = presynaptic cell
       k = com_TCR_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_deepng 
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepng(k,L)  = gAMPA_deepng(k,L) +
     &  gAMPA_TCR_to_deepng * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_TCR_to_deepng * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_deepng 
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_TCR_to_deepng * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_deepng  
       if (gNMDA_deepng(k,L).gt.z)
     &  gNMDA_deepng(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle nontuftRS  -> deepng
      do i = 1, num_nontuftRS_to_deepng
       j = map_nontuftRS_to_deepng(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_deepng
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepng(k,L)  = gAMPA_deepng(k,L) +
     &  gAMPA_nontuftRS_to_deepng * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_nontuftRS_to_deepng * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_deepng
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_nontuftRS_to_deepng * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_deepng
       if (gNMDA_deepng(k,L).gt.z)
     &  gNMDA_deepng  (k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of deepng      
         ENDIF ! if (mod(O,how_often).eq.0) ...

! Define currents to deepng       cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for deepng     cells
       CALL INTEGRATE_deepng  (O, time, num_deepng ,
     &    V_deepng , curr_deepng ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepng, gNMDA_deepng, gGABA_A_deepng,
     & Mg, 
     & gapcon_deepng  ,totSDgj_deepng   ,gjtable_deepng, dt,
     &  chi_deepng,mnaf_deepng,mnap_deepng,
     &  hnaf_deepng,mkdr_deepng,mka_deepng,
     &  hka_deepng,mk2_deepng,hk2_deepng,
     &  mkm_deepng,mkc_deepng,mkahp_deepng,
     &  mcat_deepng,hcat_deepng,mcal_deepng,
     &  mar_deepng)

      if (mod(O,5).eq.0) then


       do L = firstcell, lastcell
c           k = table_axgjcompallow_suppyrRS(1)
            k = numcomp_deepng                  
            gj_axon_deepng (L-firstcell+1) = V_deepng (k,L)
       end do
  
           call mpi_allgather (    gj_axon_deepng,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_deepng    
       do L = firstcell, lastcell
        distal_axon_deepng     (L-firstcell+1) = V_deepng     (59,L)
       end do
  
           call mpi_allgather (distal_axon_deepng,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = 0.d0     
        field_2mm_local(1) = 0.d0     
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
           ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for deepng  

c      ELSE IF (THISNO.EQ.10) THEN
       ELSE IF (nodecell(thisno) .eq. 'deepaxax ') THEN
c deepaxax

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_deepaxax
          lastcell = firstcell - 1 + ncellspernode_deepaxax

        IF (mod(O,how_often).eq.0) then
c 1st set deepaxax  synaptic conductances to 0:

          do i = 1, numcomp_deepaxax
          do j = firstcell, lastcell
         gAMPA_deepaxax(i,j)    = 0.d0
         gNMDA_deepaxax(i,j)    = 0.d0
         gGABA_A_deepaxax(i,j)  = 0.d0 
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS   -> deepaxax
      do i = 1, num_suppyrRS_to_deepaxax  
       j = map_suppyrRS_to_deepaxax(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_deepaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_deepaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepaxax(k,L)  = gAMPA_deepaxax(k,L) +
     &  gAMPA_suppyrRS_to_deepaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_suppyrRS_to_deepaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_deepaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_suppyrRS_to_deepaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_deepaxax  
       if (gNMDA_deepaxax(k,L).gt.z)
     &  gNMDA_deepaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB  -> deepaxax
      do i = 1, num_suppyrFRB_to_deepaxax  
       j = map_suppyrFRB_to_deepaxax(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_deepaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_deepaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepaxax(k,L)  = gAMPA_deepaxax(k,L) +
     &  gAMPA_suppyrFRB_to_deepaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_suppyrFRB_to_deepaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_deepaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_suppyrFRB_to_deepaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_deepaxax  
       if (gNMDA_deepaxax(k,L).gt.z)
     &  gNMDA_deepaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supLTS     -> deepaxax
      do i = 1, num_supLTS_to_deepaxax    
       j = map_supLTS_to_deepaxax(i,L) ! j = presynaptic cell
       k = com_supLTS_to_deepaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_deepaxax    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepaxax(k,L)  = gGABA_A_deepaxax(k,L) +
     &  gGABA_supLTS_to_deepaxax * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle spinstell  -> deepaxax
      do i = 1, num_spinstell_to_deepaxax   
       j = map_spinstell_to_deepaxax(i,L) ! j = presynaptic cell
       k = com_spinstell_to_deepaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_deepaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepaxax(k,L)  = gAMPA_deepaxax(k,L) +
     &  gAMPA_spinstell_to_deepaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_spinstell_to_deepaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_deepaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_spinstell_to_deepaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_deepaxax  
       if (gNMDA_deepaxax(k,L).gt.z)
     &  gNMDA_deepaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB     -> deepaxax
      do i = 1, num_tuftIB_to_deepaxax   
       j = map_tuftIB_to_deepaxax(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_deepaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_deepaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepaxax(k,L)  = gAMPA_deepaxax(k,L) +
     &  gAMPA_tuftIB_to_deepaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_tuftIB_to_deepaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_deepaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_tuftIB_to_deepaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_deepaxax  
       if (gNMDA_deepaxax(k,L).gt.z)
     &  gNMDA_deepaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS     -> deepaxax
      do i = 1, num_tuftRS_to_deepaxax   
       j = map_tuftRS_to_deepaxax(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_deepaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_deepaxax  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepaxax(k,L)  = gAMPA_deepaxax(k,L) +
     &  gAMPA_tuftRS_to_deepaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_tuftRS_to_deepaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_deepaxax  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_tuftRS_to_deepaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_deepaxax  
       if (gNMDA_deepaxax(k,L).gt.z)
     &  gNMDA_deepaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask   -> deepaxax
      do i = 1, num_deepbask_to_deepaxax    
       j = map_deepbask_to_deepaxax(i,L) ! j = presynaptic cell
       k = com_deepbask_to_deepaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_deepaxax    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepaxax(k,L)  = gGABA_A_deepaxax(k,L) +
     &  gGABA_deepbask_to_deepaxax * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepLTS    -> deepaxax
      do i = 1, num_deepLTS_to_deepaxax    
       j = map_deepLTS_to_deepaxax(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_deepaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_deepaxax    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepaxax(k,L)  = gGABA_A_deepaxax(k,L) +
     &  gGABA_deepLTS_to_deepaxax * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle TCR        -> deepaxax
      do i = 1, num_TCR_to_deepaxax 
       j = map_TCR_to_deepaxax(i,L) ! j = presynaptic cell
       k = com_TCR_to_deepaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_TCR_to_deepaxax 
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepaxax(k,L)  = gAMPA_deepaxax(k,L) +
     &  gAMPA_TCR_to_deepaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_TCR_to_deepaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_deepaxax 
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_TCR_to_deepaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_deepaxax  
       if (gNMDA_deepaxax(k,L).gt.z)
     &  gNMDA_deepaxax(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle nontuftRS  -> deepaxax
      do i = 1, num_nontuftRS_to_deepaxax
       j = map_nontuftRS_to_deepaxax(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_deepaxax(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_deepaxax
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepaxax(k,L)  = gAMPA_deepaxax(k,L) +
     &  gAMPA_nontuftRS_to_deepaxax * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_nontuftRS_to_deepaxax * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_deepaxax
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepaxax(k,L) = gNMDA_deepaxax(k,L) +
     &  gNMDA_nontuftRS_to_deepaxax * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_deepaxax
       if (gNMDA_deepaxax(k,L).gt.z)
     &  gNMDA_deepaxax(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of deepaxax    
        ENDIF  !  if (mod(O,how_often).eq.0) ...

! Define currents to deepaxax     cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for deepaxax     cells
       CALL INTEGRATE_deepaxaxx (O, time, num_deepaxax ,
     &    V_deepaxax , curr_deepaxax ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepaxax, gNMDA_deepaxax, gGABA_A_deepaxax,
     & Mg, 
     & gapcon_deepaxax  ,totSDgj_deepaxax   ,gjtable_deepaxax, dt,
     &  chi_deepaxax,mnaf_deepaxax,mnap_deepaxax,
     &  hnaf_deepaxax,mkdr_deepaxax,mka_deepaxax,
     &  hka_deepaxax,mk2_deepaxax,hk2_deepaxax,
     &  mkm_deepaxax,mkc_deepaxax,mkahp_deepaxax,
     &  mcat_deepaxax,hcat_deepaxax,mcal_deepaxax,
     &  mar_deepaxax)

      if (mod(O,5).eq.0) then

c     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
c    &    mpi_double_precision, 0, mpi_comm_world, info)
c     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
c    &    mpi_double_precision, 1, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
c    &    mpi_double_precision, 6, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
c    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
c           k = table_axgjcompallow_suppyrRS(1)
            k = numcomp_deepaxax                  
            gj_axon_deepaxax (L-firstcell+1) = V_deepaxax (k,L)
       end do
  
           call mpi_allgather (    gj_axon_deepaxax,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_deepaxax    
       do L = firstcell, lastcell
        distal_axon_deepaxax     (L-firstcell+1) = V_deepaxax     (59,L)
       end do
  
           call mpi_allgather (distal_axon_deepaxax,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = 0.d0     
        field_2mm_local(1) = 0.d0     
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
        ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for deepaxax

c      ELSE IF (THISNO.EQ.11) THEN
       ELSE IF (nodecell(thisno) .eq. 'deepLTS  ') THEN
c deepLTS

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_deepLTS 
          lastcell = firstcell - 1 + ncellspernode_deepLTS 

       IF (mod(O,how_often).eq.0) then
c 1st set deepLTS   synaptic conductances to 0:

          do i = 1, numcomp_deepLTS
          do j = firstcell, lastcell
         gAMPA_deepLTS(i,j)     = 0.d0
         gNMDA_deepLTS(i,j)     = 0.d0
         gGABA_A_deepLTS(i,j)   = 0.d0 
          end do
          end do

         do L = firstcell, lastcell
c Handle suppyrRS   -> deepLTS
      do i = 1, num_suppyrRS_to_deepLTS   
       j = map_suppyrRS_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_suppyrRS_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrRS_to_deepLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_suppyrRS_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_suppyrRS_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrRS_to_deepLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_suppyrRS_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrRS_to_deepLTS  
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle suppyrFRB  -> deepLTS
      do i = 1, num_suppyrFRB_to_deepLTS   
       j = map_suppyrFRB_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_suppyrFRB_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_suppyrFRB(j)  ! enumerate presyn. spikes
        presyntime = outtime_suppyrFRB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_suppyrFRB_to_deepLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_suppyrFRB_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_suppyrFRB_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_suppyrFRB_to_deepLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_suppyrFRB_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_suppyrFRB_to_deepLTS  
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supLTS     -> deepLTS
      do i = 1, num_supLTS_to_deepLTS     
       j = map_supLTS_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_supLTS_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_supLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supLTS_to_deepLTS     
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepLTS(k,L)  = gGABA_A_deepLTS(k,L) +
     &  gGABA_supLTS_to_deepLTS * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle spinstell  -> deepLTS
      do i = 1, num_spinstell_to_deepLTS    
       j = map_spinstell_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_spinstell_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_spinstell(j)  ! enumerate presyn. spikes
        presyntime = outtime_spinstell(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_spinstell_to_deepLTS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_spinstell_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_spinstell_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_spinstell_to_deepLTS   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_spinstell_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_spinstell_to_deepLTS  
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftIB     -> deepLTS
      do i = 1, num_tuftIB_to_deepLTS    
       j = map_tuftIB_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_tuftIB_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftIB(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftIB(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftIB_to_deepLTS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_tuftIB_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_tuftIB_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftIB_to_deepLTS   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_tuftIB_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftIB_to_deepLTS   
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle tuftRS     -> deepLTS
      do i = 1, num_tuftRS_to_deepLTS    
       j = map_tuftRS_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_tuftRS_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_tuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_tuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_tuftRS_to_deepLTS   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_tuftRS_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_tuftRS_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_tuftRS_to_deepLTS   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_tuftRS_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_tuftRS_to_deepLTS   
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask   -> deepLTS
      do i = 1, num_deepbask_to_deepLTS     
       j = map_deepbask_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_deepbask_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_deepLTS     
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepLTS(k,L)  = gGABA_A_deepLTS(k,L) +
     &  gGABA_deepbask_to_deepLTS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepLTS    -> deepLTS
      do i = 1, num_deepLTS_to_deepLTS     
       j = map_deepLTS_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_deepLTS     
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepLTS(k,L)  = gGABA_A_deepLTS(k,L) +
     &  gGABA_deepLTS_to_deepLTS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle nontuftRS  -> deepLTS
      do i = 1, num_nontuftRS_to_deepLTS
       j = map_nontuftRS_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_deepLTS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_nontuftRS_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_nontuftRS_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_deepLTS 
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_nontuftRS_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_deepLTS
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of deepLTS     
         ENDIF  !  if (mod(O,how_often).eq.0) ...

! Define currents to deepLTS      cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for deepLTS      cells
       CALL INTEGRATE_deepLTSx  (O, time, num_deepLTS  ,
     &    V_deepLTS  , curr_deepLTS  ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepLTS  , gNMDA_deepLTS  , gGABA_A_deepLTS  ,
     & Mg, 
     & gapcon_deepLTS  ,totSDgj_deepLTS  ,gjtable_deepLTS  , dt,
     &  chi_deepLTS,mnaf_deepLTS,mnap_deepLTS,
     &  hnaf_deepLTS,mkdr_deepLTS,mka_deepLTS,
     &  hka_deepLTS,mk2_deepLTS,hk2_deepLTS,
     &  mkm_deepLTS,mkc_deepLTS,mkahp_deepLTS,
     &  mcat_deepLTS,hcat_deepLTS,mcal_deepLTS,
     &  mar_deepLTS)

      if (mod(O,5).eq.0) then

c     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
c    &    mpi_double_precision, 0, mpi_comm_world, info)
c     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
c    &    mpi_double_precision, 1, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
c    &    mpi_double_precision, 6, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
c    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
c           k = table_axgjcompallow_suppyrRS(1)
            k = numcomp_deepLTS                   
            gj_axon_deepLTS  (L-firstcell+1) = V_deepLTS  (k,L)
       end do
  
           call mpi_allgather (    gj_axon_deepLTS ,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_deepLTS     
       do L = firstcell, lastcell
        distal_axon_deepLTS      (L-firstcell+1) = V_deepLTS      (59,L)
       end do
  
           call mpi_allgather (distal_axon_deepLTS, 
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = 0.d0     
        field_2mm_local(1) = 0.d0     
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
        ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for deepLTS

c      ELSE IF (THISNO.EQ.12) THEN
       ELSE IF (nodecell(thisno) .eq. 'TCR      ') THEN
c TCR

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_TCR     
          lastcell = firstcell - 1 + ncellspernode_TCR     

        IF (mod(O,how_often).eq.0) then
c 1st set TCR synaptic conductances to 0:

          do i = 1, numcomp_TCR
          do j = firstcell, lastcell
         gAMPA_TCR(i,j)         = 0.d0 
         gNMDA_TCR(i,j)         = 0.d0
         gGABA_A_TCR(i,j)       = 0.d0 
         gGABA_B_TCR(i,j)       = 0.d0 
          end do
          end do

         do L = firstcell, lastcell
c Handle nRT       -> TCR
      do i = 1, num_nRT_to_TCR     
       j = map_nRT_to_TCR(i,L) ! j = presynaptic cell
       k = com_nRT_to_TCR(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nRT(j)  ! enumerate presyn. spikes
        presyntime = outtime_nRT(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part
        dexparg1 = delta / tauGABA1_nRT_to_TCR     
c note that dexparg1 = MINUS the actual arg. to dexp
         if (dexparg1.le.5.d0) then
          z1 = dexptablesmall (int(dexparg1*1000.d0))
         else if (dexparg1.le.100.d0) then
          z1 = dexptablebig (int(dexparg1*10.d0))
         else
          z1 = 0.d0
         endif

        dexparg2 = delta / tauGABA2_nRT_to_TCR     
c note that dexparg2 = MINUS the actual arg. to dexp
         if (dexparg2.le.5.d0) then
          z2 = dexptablesmall (int(dexparg2*1000.d0))
         else if (dexparg2.le.100.d0) then
          z2 = dexptablebig (int(dexparg2*10.d0))
         else
          z2 = 0.d0
         endif

      gGABA_A_TCR(k,L)  = gGABA_A_TCR(k,L) +
     &  gGABA_nRT_to_TCR(j) * (0.625d0 * z1 + 0.375d0 * z2) 
! end GABA-A part


      gGABA_B_TCR(k,L) = gGABA_B_TCR(k,L) +
     &   gGABAB_nRT_to_TCR * otis_table(k0)
! end GABA-B part
       end do ! m
      end do ! i


c Handle nontuftRS -> TCR
      do i = 1, num_nontuftRS_to_TCR
       j = map_nontuftRS_to_TCR(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_TCR(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime - cort_thal_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_TCR
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_TCR(k,L)  = gAMPA_TCR(k,L) +
     &  gAMPA_nontuftRS_to_TCR * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_TCR(k,L) = gNMDA_TCR(k,L) +
     &  gNMDA_nontuftRS_to_TCR * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_TCR
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_TCR(k,L) = gNMDA_TCR(k,L) +
     &  gNMDA_nontuftRS_to_TCR * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_TCR 
       if (gNMDA_TCR(k,L).gt.z)
     &  gNMDA_TCR(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


         end do
c End enumeration of TCR         
          ENDIF  !  if (mod(O,how_often).eq.0) ...

! Define currents to TCR          cells, ectopic spikes,
! tonic synaptic conductances

      if (mod(O,200).eq.0) then
       call durand(seed,num_TCR     ,ranvec_TCR     ) 
        do L = firstcell, lastcell
         if ((ranvec_TCR     (L).gt.0.d0).and.
     &     (ranvec_TCR     (L).le.noisepe_TCR     )) then
          curr_TCR     (135,L) = 0.4d0
          ectr_TCR      = ectr_TCR      + 1
         else
          curr_TCR     (135,L) = 0.d0
         endif 
        end do
      endif

c       GOTO 9144 ! SKIP TCR INTEGRATION IN ISOLATED CTX
! Call integration routine for TCR          cells
       CALL INTEGRATE_tcrxB     (O, time, num_tcr      ,
     &    V_tcr      , curr_tcr      ,
     & initialize, firstcell, lastcell,
     & gAMPA_tcr      , gNMDA_tcr      , gGABA_A_tcr      ,
     & gGABA_B_tcr, Mg, 
     & gapcon_tcr      ,totaxgj_tcr      ,gjtable_tcr      , dt,
     &  chi_tcr,mnaf_tcr,mnap_tcr,
     &  hnaf_tcr,mkdr_tcr,mka_tcr,
     &  hka_tcr,mk2_tcr,hk2_tcr,
     &  mkm_tcr,mkc_tcr,mkahp_tcr,
     &  mcat_tcr,hcat_tcr,mcal_tcr,
     &  mar_tcr)
9144    CONTINUE

      if (mod(O,5).eq.0) then

c     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
c    &    mpi_double_precision, 0, mpi_comm_world, info)
c     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
c    &    mpi_double_precision, 1, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
c    &    mpi_double_precision, 6, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
c    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
            k = table_axgjcompallow_TCR(1)
            gj_axon_TCR      (L-firstcell+1) = V_TCR      (k,L)
       end do
  
           call mpi_allgather (    gj_axon_TCR     ,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

         IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_TCR         
       do L = firstcell, lastcell
        distal_axon_TCR  (L-firstcell+1) = V_TCR          (135,L)
       end do
  
           call mpi_allgather (distal_axon_TCR,     
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = 0.d0       
        field_2mm_local(1) = 0.d0     
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
        ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for TCR

c      ELSE IF (THISNO.EQ.13) THEN
       ELSE IF (nodecell(thisno) .eq. 'nRT      ') THEN
c nRT

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_nRT     
          lastcell = firstcell - 1 + ncellspernode_nRT     

        IF (mod(O,how_often).eq.0) then
c 1st set nRT synaptic conductances to 0:

          do i = 1, numcomp_nRT
          do j = firstcell, lastcell
         gAMPA_nRT(i,j)         = 0.d0 
         gNMDA_nRT(i,j)         = 0.d0
         gGABA_A_nRT(i,j)       = 0.d0
         gGABA_B_nRT(i,j)       = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle TCR        -> nRT
      do i = 1, num_TCR_to_nRT
       j = map_TCR_to_nRT(i,L) ! j = presynaptic cell
       k = com_TCR_to_nRT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_TCR(j)  ! enumerate presyn. spikes
        presyntime = outtime_TCR(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_TCR_to_nRT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_nRT(k,L)  = gAMPA_nRT(k,L) +
     &  gAMPA_TCR_to_nRT * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_nRT(k,L) = gNMDA_nRT(k,L) +
     &  gNMDA_TCR_to_nRT * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_TCR_to_nRT 
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_nRT(k,L) = gNMDA_nRT(k,L) +
     &  gNMDA_TCR_to_nRT * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_TCR_to_nRT
       if (gNMDA_nRT(k,L).gt.z)
     &  gNMDA_nRT(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle nRT        -> nRT
      do i = 1, num_nRT_to_nRT     
       j = map_nRT_to_nRT(i,L) ! j = presynaptic cell
       k = com_nRT_to_nRT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nRT(j)  ! enumerate presyn. spikes
        presyntime = outtime_nRT(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part
        dexparg1 = delta / tauGABA1_nRT_to_nRT     
c note that dexparg1 = MINUS the actual arg. to dexp
         if (dexparg1.le.5.d0) then
          z1 = dexptablesmall (int(dexparg1*1000.d0))
         else if (dexparg1.le.100.d0) then
          z1 = dexptablebig (int(dexparg1*10.d0))
         else
          z1 = 0.d0
         endif

        dexparg2 = delta / tauGABA2_nRT_to_nRT     
c note that dexparg2 = MINUS the actual arg. to dexp
         if (dexparg2.le.5.d0) then
          z2 = dexptablesmall (int(dexparg2*1000.d0))
         else if (dexparg2.le.100.d0) then
          z2 = dexptablebig (int(dexparg2*10.d0))
         else
          z2 = 0.d0
         endif

      gGABA_A_nRT(k,L)  = gGABA_A_nRT(k,L) +
     &  gGABA_nRT_to_nRT * (0.56d0 * z1 + 0.44d0 * z2) 
! end GABA-A part

      gGABA_B_nRT(k,L) = gGABA_B_nRT(k,L) +
     &   gGABAB_nRT_to_nRT * otis_table(k0)

! end GABA-B part
       end do ! m
      end do ! i


c Handle nontuftRS  -> nRT
      do i = 1, num_nontuftRS_to_nRT
       j = map_nontuftRS_to_nRT(i,L) ! j = presynaptic cell
       k = com_nontuftRS_to_nRT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_nontuftRS(j)  ! enumerate presyn. spikes
        presyntime = outtime_nontuftRS(m,j)
        delta = time - presyntime - cort_thal_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_nontuftRS_to_nRT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_nRT(k,L)  = gAMPA_nRT(k,L) +
     &  gAMPA_nontuftRS_to_nRT * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_nRT(k,L) = gNMDA_nRT(k,L) +
     &  gNMDA_nontuftRS_to_nRT * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_nontuftRS_to_nRT
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_nRT(k,L) = gNMDA_nRT(k,L) +
     &  gNMDA_nontuftRS_to_nRT * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_nontuftRS_to_nRT 
       if (gNMDA_nRT(k,L).gt.z)
     &  gNMDA_nRT(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


         end do
c End enumeration of nRT         
        ENDIF  !  if (mod(O,how_often).eq.0) ...

! Define currents to nRT          cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for nRT          cells
       CALL INTEGRATE_nrtxB      (O, time, num_nRT      ,
     &    V_nRT      , curr_nRT      ,
     & initialize, firstcell, lastcell,
     & gAMPA_nRT      , gNMDA_nRT      , gGABA_A_nRT      ,
     & gGABA_B_nRT, Mg, 
     & gapcon_nRT      ,totSDgj_nRT      ,gjtable_nRT      , dt,
     &  chi_nRT,mnaf_nRT,mnap_nRT,
     &  hnaf_nRT,mkdr_nRT,mka_nRT,
     &  hka_nRT,mk2_nRT,hk2_nRT,
     &  mkm_nRT,mkc_nRT,mkahp_nRT,
     &  mcat_nRT,hcat_nRT,mcal_nRT,
     &  mar_nRT)

      if (mod(O,5).eq.0) then

c     call mpi_bcast (vax_suppyrRS  ,num_suppyrRS,
c    &    mpi_double_precision, 0, mpi_comm_world, info)
c     call mpi_bcast (vax_suppyrFRB ,num_suppyrFRB,
c    &    mpi_double_precision, 1, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftIB    ,num_tuftIB  ,
c    &    mpi_double_precision, 6, mpi_comm_world, info)
c     call mpi_bcast (vax_tuftRS    ,num_tuftRS   ,
c    &    mpi_double_precision, 7, mpi_comm_world, info)

       do L = firstcell, lastcell
c           k = table_axgjcompallow_suppyrRS(1)
            k = numcomp_nRT                      
            gj_axon_nRT      (L-firstcell+1) = V_nRT      (k,L)
       end do
  
           call mpi_allgather (    gj_axon_nRT     ,
     &  maxcellspernode, mpi_double_precision,
     &      gj_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

      endif

         IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_nRT         
       do L = firstcell, lastcell
        distal_axon_nRT  (L-firstcell+1) = V_nRT          (59,L)
       end do
  
           call mpi_allgather (distal_axon_nRT,      
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_1mm_local(1) = 0.d0      
        field_2mm_local(1) = 0.d0       
           call mpi_allgather (field_1mm_local,     
     &  1              , mpi_double_precision,
     &  field_1mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_2mm_local,     
     &  1              , mpi_double_precision,
     &  field_2mm_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
         ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for nRT

       ENDIF  ! if (mod(O,how_often).eq.0) then ...

! Update distal axon vectors, then outctr's and outtime tables.
! This code is common to all the nodes.
! Some of this section adapted from supergj.f
c     IF (mod(O,how_often).eq.0) then
      IF (mod(O,  5      ).eq.0) then ! Necessary because gj data also
!  being updated, not just synaptic
c Construct distal axon vectors, taking into account the structure of
c distal_axon_global: let m = maxcellspernode;
c then nodesfor_suppyrRS segments, each m entries long;
c following this, nodesfor_suppyrFRB segments, each m entries long; etc....
! Do the same for voltages at sites of possible axonal gj.

            ictr = 0 ! will keep track of which segment in distal_axon_global

            j = table_axgjcompallow_suppyrRS (1)
            do k = 1, nodesfor_suppyrRS
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_suppyrRS  ! offset in ldistal_axon_suppyrRS
        do i = m, m + ncellspernode_suppyrRS - 1 
         ldistal_axon_suppyrRS(i   ) = distal_axon_global(ioffset+i-m)
         v_suppyrRS  (j,i)          =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = table_axgjcompallow_suppyrFRB(1)
            do k = 1, nodesfor_suppyrFRB
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_suppyrFRB ! offset in ldistal_axon_suppyrFRB
        do i = m, m + ncellspernode_suppyrFRB - 1 
         ldistal_axon_suppyrFRB(i   ) = distal_axon_global(ioffset+i-m)
         v_suppyrFRB (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = numcomp_supbask
            do k = 1, nodesfor_supbask  
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_supbask   ! offset in ldistal_axon_supbask  
        do i = m, m + ncellspernode_supbask   - 1 
         ldistal_axon_supbask  (i   ) = distal_axon_global(ioffset+i-m)
         v_supbask   (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = numcomp_supng
            do k = 1, nodesfor_supng  
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_supng   ! offset in ldistal_axon_supng  
        do i = m, m + ncellspernode_supng   - 1 
         ldistal_axon_supng  (i   ) = distal_axon_global(ioffset+i-m)
         v_supng   (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = numcomp_supaxax
            do k = 1, nodesfor_supaxax  
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_supaxax   ! offset in ldistal_axon_supaxax  
        do i = m, m + ncellspernode_supaxax   - 1 
         ldistal_axon_supaxax  (i   ) = distal_axon_global(ioffset+i-m)
         v_supaxax   (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = numcomp_supLTS
            do k = 1, nodesfor_supLTS   
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_supLTS    ! offset in ldistal_axon_supLTS   
        do i = m, m + ncellspernode_supLTS    - 1 
         ldistal_axon_supLTS   (i   ) = distal_axon_global(ioffset+i-m)
         v_supLTS    (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = table_axgjcompallow_spinstell(1)
            do k = 1, nodesfor_spinstell
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_spinstell ! offset in ldistal_axon_spinstell
        do i = m, m + ncellspernode_spinstell - 1 
         ldistal_axon_spinstell(i   ) = distal_axon_global(ioffset+i-m)
         v_spinstell (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = table_axgjcompallow_tuftIB   (1)
            do k = 1, nodesfor_tuftIB   
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_tuftIB    ! offset in ldistal_axon_tuftIB   
        do i = m, m + ncellspernode_tuftIB    - 1 
         ldistal_axon_tuftIB   (i   ) = distal_axon_global(ioffset+i-m)
         v_tuftIB    (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = table_axgjcompallow_tuftRS   (1)
            do k = 1, nodesfor_tuftRS   
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_tuftRS    ! offset in ldistal_axon_tuftRS   
        do i = m, m + ncellspernode_tuftRS    - 1 
         ldistal_axon_tuftRS   (i   ) = distal_axon_global(ioffset+i-m)
         v_tuftRS    (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = table_axgjcompallow_nontuftRS(1)
            do k = 1, nodesfor_nontuftRS   
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_nontuftRS    ! offset in ldistal_axon_nontuftRS   
        do i = m, m + ncellspernode_nontuftRS    - 1 
         ldistal_axon_nontuftRS  (i ) = distal_axon_global(ioffset+i-m)
         v_nontuftRS (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = numcomp_deepbask
            do k = 1, nodesfor_deepbask    
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_deepbask     ! offset in ldistal_axon_deepbask    
        do i = m, m + ncellspernode_deepbask     - 1 
         ldistal_axon_deepbask   (i ) = distal_axon_global(ioffset+i-m)
         v_deepbask  (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = numcomp_deepng
            do k = 1, nodesfor_deepng    
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_deepng     ! offset in ldistal_axon_deepng    
        do i = m, m + ncellspernode_deepng     - 1 
         ldistal_axon_deepng   (i ) = distal_axon_global(ioffset+i-m)
         v_deepng  (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = numcomp_deepaxax
            do k = 1, nodesfor_deepaxax    
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_deepaxax     ! offset in ldistal_axon_deepaxax    
        do i = m, m + ncellspernode_deepaxax     - 1 
         ldistal_axon_deepaxax   (i ) = distal_axon_global(ioffset+i-m)
         v_deepaxax  (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = numcomp_deepLTS
            do k = 1, nodesfor_deepLTS     
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_deepLTS      ! offset in ldistal_axon_deepLTS     
        do i = m, m + ncellspernode_deepLTS      - 1 
         ldistal_axon_deepLTS    (i ) = distal_axon_global(ioffset+i-m)
         v_deepLTS   (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = table_axgjcompallow_TCR      (1)
            do k = 1, nodesfor_TCR         
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_TCR          ! offset in ldistal_axon_TCR         
        do i = m, m + ncellspernode_TCR          - 1 
         ldistal_axon_TCR        (i ) = distal_axon_global(ioffset+i-m)
         v_TCR       (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do

            j = numcomp_nRT
            do k = 1, nodesfor_nRT         
              ictr = ictr + 1
              ioffset = (ictr - 1) * maxcellspernode + 1 ! offset in distal_axon_global
              m = 1 + (k-1) * ncellspernode_nRT          ! offset in ldistal_axon_nRT         
        do i = m, m + ncellspernode_nRT          - 1 
         ldistal_axon_nRT        (i ) = distal_axon_global(ioffset+i-m)
         v_nRT       (j,i)           =     gj_axon_global(ioffset+i-m)
        end do
             end do
c End updating of distal axon vectors.

       do L = 1, num_suppyrRS
	 if (ldistal_axon_suppyrRS(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_suppyrRS(L).eq.0) then
	    outctr_suppyrRS(L) = 1
	    outtime_suppyrRS(1,L) = time
          else
      if ((time-outtime_suppyrRS(outctr_suppyrRS(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_suppyrRS(L) = outctr_suppyrRS(L) + 1
	     outtime_suppyrRS (outctr_suppyrRS(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_suppyrRS

       do L = 1, num_suppyrFRB
	 if (ldistal_axon_suppyrFRB(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_suppyrFRB(L).eq.0) then
	    outctr_suppyrFRB(L) = 1
	    outtime_suppyrFRB(1,L) = time
          else
      if ((time-outtime_suppyrFRB(outctr_suppyrFRB(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_suppyrFRB(L) = outctr_suppyrFRB(L) + 1
	     outtime_suppyrFRB (outctr_suppyrFRB(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_suppyrFRB

       do L = 1, num_supbask   
	 if (ldistal_axon_supbask(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_supbask(L).eq.0) then
	    outctr_supbask(L) = 1
	    outtime_supbask(1,L) = time
          else
      if ((time-outtime_supbask(outctr_supbask(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_supbask(L) = outctr_supbask(L) + 1
	     outtime_supbask (outctr_supbask(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_supbask  

       do L = 1, num_supng   
	 if (ldistal_axon_supng(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_supng(L).eq.0) then
	    outctr_supng(L) = 1
	    outtime_supng(1,L) = time
          else
      if ((time-outtime_supng(outctr_supng(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_supng(L) = outctr_supng(L) + 1
	     outtime_supng (outctr_supng(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_supng  

       do L = 1, num_supaxax   
	 if (ldistal_axon_supaxax(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_supaxax(L).eq.0) then
	    outctr_supaxax(L) = 1
	    outtime_supaxax(1,L) = time
          else
      if ((time-outtime_supaxax(outctr_supaxax(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_supaxax(L) = outctr_supaxax(L) + 1
	     outtime_supaxax (outctr_supaxax(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_supaxax  

       do L = 1, num_supLTS    
	 if (ldistal_axon_supLTS(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_supLTS(L).eq.0) then
	    outctr_supLTS(L) = 1
	    outtime_supLTS(1,L) = time
          else
      if ((time-outtime_supLTS(outctr_supLTS(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_supLTS(L) = outctr_supLTS(L) + 1
	     outtime_supLTS (outctr_supLTS(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_supLTS  

       do L = 1, num_spinstell 
	 if (ldistal_axon_spinstell(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_spinstell(L).eq.0) then
	    outctr_spinstell(L) = 1
	    outtime_spinstell(1,L) = time
          else
      if ((time-outtime_spinstell(outctr_spinstell(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_spinstell(L) = outctr_spinstell(L) + 1
	     outtime_spinstell (outctr_spinstell(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_spinstell

       do L = 1, num_tuftIB    
	 if (ldistal_axon_tuftIB(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_tuftIB(L).eq.0) then
	    outctr_tuftIB(L) = 1
	    outtime_tuftIB(1,L) = time
          else
      if ((time-outtime_tuftIB(outctr_tuftIB(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_tuftIB(L) = outctr_tuftIB(L) + 1
	     outtime_tuftIB (outctr_tuftIB(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_tuftIB   

       do L = 1, num_tuftRS    
	 if (ldistal_axon_tuftRS(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_tuftRS(L).eq.0) then
	    outctr_tuftRS(L) = 1
	    outtime_tuftRS(1,L) = time
          else
      if ((time-outtime_tuftRS(outctr_tuftRS(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_tuftRS(L) = outctr_tuftRS(L) + 1
	     outtime_tuftRS (outctr_tuftRS(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_tuftRS   

       do L = 1, num_nontuftRS    
	 if (ldistal_axon_nontuftRS(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_nontuftRS(L).eq.0) then
	    outctr_nontuftRS(L) = 1
	    outtime_nontuftRS(1,L) = time
          else
      if ((time-outtime_nontuftRS(outctr_nontuftRS(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_nontuftRS(L) = outctr_nontuftRS(L) + 1
	     outtime_nontuftRS (outctr_nontuftRS(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_nontuftRS   

       do L = 1, num_deepbask     
	 if (ldistal_axon_deepbask(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_deepbask(L).eq.0) then
	    outctr_deepbask(L) = 1
	    outtime_deepbask(1,L) = time
          else
      if ((time-outtime_deepbask(outctr_deepbask(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_deepbask(L) = outctr_deepbask(L) + 1
	     outtime_deepbask (outctr_deepbask(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_deepbask   

       do L = 1, num_deepng     
	 if (ldistal_axon_deepng(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_deepng(L).eq.0) then
	    outctr_deepng(L) = 1
	    outtime_deepng(1,L) = time
          else
      if ((time-outtime_deepng(outctr_deepng(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_deepng(L) = outctr_deepng(L) + 1
	     outtime_deepng (outctr_deepng(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_deepng   

       do L = 1, num_deepaxax     
	 if (ldistal_axon_deepaxax(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_deepaxax(L).eq.0) then
	    outctr_deepaxax(L) = 1
	    outtime_deepaxax(1,L) = time
          else
      if ((time-outtime_deepaxax(outctr_deepaxax(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_deepaxax(L) = outctr_deepaxax(L) + 1
	     outtime_deepaxax (outctr_deepaxax(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_deepaxax   

       do L = 1, num_deepLTS      
	 if (ldistal_axon_deepLTS(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_deepLTS(L).eq.0) then
	    outctr_deepLTS(L) = 1
	    outtime_deepLTS(1,L) = time
          else
      if ((time-outtime_deepLTS(outctr_deepLTS(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_deepLTS(L) = outctr_deepLTS(L) + 1
	     outtime_deepLTS (outctr_deepLTS(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_deepLTS   

       do L = 1, num_TCR      
	 if (ldistal_axon_TCR(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_TCR(L).eq.0) then
	    outctr_TCR(L) = 1
	    outtime_TCR(1,L) = time
          else
      if ((time-outtime_TCR(outctr_TCR(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_TCR(L) = outctr_TCR(L) + 1
	     outtime_TCR (outctr_TCR(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_TCR   

       do L = 1, num_nRT      
	 if (ldistal_axon_nRT(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
	  if (outctr_nRT(L).eq.0) then
	    outctr_nRT(L) = 1
	    outtime_nRT(1,L) = time
          else
      if ((time-outtime_nRT(outctr_nRT(L),L))
     &   .gt. axon_refrac_time) then
	     outctr_nRT(L) = outctr_nRT(L) + 1
	     outtime_nRT (outctr_nRT(L),L) = time
            endif
          endif
	 endif
       end do  ! do L = 1, num_nRT   

       field_1mm_tot = 0.d0
       field_2mm_tot = 0.d0
        do i = 1, numnodes
         field_1mm_tot = field_1mm_tot + field_1mm_global(i)
         field_2mm_tot = field_2mm_tot + field_2mm_global(i)
        end do

      ENDIF  ! if (mod(O,how_often).eq.0) ...
       ! CHANGED to if (mod(O,5).eq.0)...
! End updating outctr's and outtime tables, and computing fields

! Set up output data to be written
       if (mod(O, 50) == 0) then
c      if (thisno.eq.0) then
       IF (nodecell(thisno) .eq. 'suppyrRS ') THEN
c suppyrRS
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_suppyrRS
          lastcell = firstcell - 1 + ncellspernode_suppyrRS

        outrcd( 1) = time
        outrcd( 2) = v_suppyrRS(1,firstcell+1)
        outrcd( 3) = v_suppyrRS(numcomp_suppyrRS,firstcell+1)
        outrcd( 4) = v_suppyrRS(43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_suppyrRS(1,i)
          end do
        outrcd( 5) = z / dble(lastcell - firstcell + 1) ! - av. cell somata 
         z = 0.d0
          do i = 1, numcomp_suppyrRS
           z = z + gAMPA_suppyrRS(i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_suppyrRS
           z = z + gNMDA_suppyrRS(i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_suppyrRS
           z = z + gGABA_A_suppyrRS(i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_suppyrRS
           z = z + gGABA_B_suppyrRS(i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_suppyrRS(1,firstcell+2)
        outrcd(11) = v_suppyrRS( 1,firstcell+3)
          z = 0.d0
          do i = firstcell, lastcell
           if(v_suppyrRS(numcomp_suppyrRS,i) .gt. 0.d0) z = z + 1.d0
          end do
        outrcd(12) = z   
        outrcd(13) = field_1mm_tot     
        outrcd(14) = field_2mm_tot     

        outrcd(15) = v_suppyrRS( 1,firstcell+4)
        outrcd(16) = v_suppyrRS( 1,firstcell+5)
c       outrcd(16) = mkahp_suppyrRS (53,firstcell+1)
c       outrcd(16) = mkDR_suppyrRS (numcomp,firstcell+1)
        outrcd(17) = mnap_suppyrRS  (69,firstcell+1)
        outrcd(18) = hka_suppyrRS (numcomp,firstcell+1)
        outrcd(19) = mkm_suppyrRS (numcomp,firstcell+1)
        outrcd(20) = chi_suppyrRS (1, firstcell+1)

            if (place(thisno).eq.1) then
      OPEN(11,FILE='alphaY67.suppyrRS')
      WRITE (11,FMT='(20F10.4)') (OUTRCD(I),I=1,20)
            end if

         outrcd( 1) = time
         outrcd( 2) = v_suppyrRS( 1,firstcell + 10)
         outrcd( 3) = v_suppyrRS(numcomp, firstcell+10)
         outrcd( 4) = v_suppyrRS(53,firstcell + 10)
         outrcd( 5) = v_suppyrRS( 1,firstcell + 11)
         outrcd( 6) = v_suppyrRS(53,firstcell + 11)
         outrcd( 7) = v_suppyrRS( 1,firstcell + 12)
         outrcd( 8) = v_suppyrRS(53,firstcell + 12)
         outrcd( 9) = v_suppyrRS( 1,firstcell + 12)
         outrcd(10) = v_suppyrRS(53,firstcell + 13)
            if (place(thisno).eq.1) then
      OPEN(81,FILE='alphaY67.suppyrRSA')
      WRITE (81,FMT='(10F10.4)') (OUTRCD(I),I=1,10)
            end if

c      else if (thisno.eq.1) then
       else IF (nodecell(thisno) .eq. 'suppyrFRB') THEN
c suppyrFRB
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_suppyrFRB
          lastcell = firstcell - 1 + ncellspernode_suppyrFRB

        outrcd( 1) = time
        outrcd( 2) = v_suppyrFRB(1,firstcell+1)
        outrcd( 3) = v_suppyrFRB(numcomp_suppyrFRB,firstcell+1)
        outrcd( 4) = v_suppyrFRB(43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_suppyrFRB(1,i)
          end do
        outrcd( 5) = z / dble(lastcell - firstcell + 1) ! - av. cell somata 
         z = 0.d0
          do i = 1, numcomp_suppyrFRB
           z = z + gAMPA_suppyrFRB(i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_suppyrFRB
           z = z + gNMDA_suppyrFRB(i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_suppyrFRB
           z = z + gGABA_A_suppyrFRB(i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_suppyrFRB
           z = z + gGABA_B_suppyrFRB(i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd(10) = v_suppyrFRB(1,firstcell+2)
        outrcd(11) = v_suppyrFRB(1,firstcell+3)
        outrcd(12) = field_1mm_tot      
        outrcd(13) = field_2mm_tot      

        if (place(thisno).eq.1) then
      OPEN(12,FILE='alphaY67.suppyrFRB')
c     WRITE (12,FMT='(13F10.4)') (OUTRCD(I),I=1,13)
        endif

c      else if (thisno.eq.2) then
       else IF (nodecell(thisno) .eq. 'supbask  ') THEN
c supbask 
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_supbask 
          lastcell = firstcell - 1 + ncellspernode_supbask 

        outrcd( 1) = time
        outrcd( 2) = v_supbask  (1,firstcell+1)
        outrcd( 3) = v_supbask  (numcomp_supbask,firstcell+1)
        outrcd( 4) = v_supbask  (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_supbask(1,i)
          end do
        outrcd( 5) = z / dble(lastcell - firstcell + 1  ) ! - av. cell somata 
         z = 0.d0
          do i = 1, numcomp_supbask   
           z = z + gAMPA_supbask  (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_supbask   
           z = z + gNMDA_supbask  (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_supbask  
           z = z + gGABA_A_supbask  (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_supbask  (1,firstcell+2)
        outrcd(10) = v_supbask  (1,firstcell+3)
        outrcd(11) = field_1mm_tot      
        outrcd(12) = field_2mm_tot      
         if (place(thisno).eq.1) then
      OPEN(13,FILE='alphaY67.supbask')
c     WRITE (13,FMT='(12F10.4)') (OUTRCD(I),I=1,12)
         endif

       else IF (nodecell(thisno) .eq. 'supng  ') THEN
c supng 
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_supng 
          lastcell = firstcell - 1 + ncellspernode_supng 

        outrcd( 1) = time
        outrcd( 2) = v_supng  (1,firstcell+1)
        outrcd( 3) = v_supng  (numcomp_supng,firstcell+1)
        outrcd( 4) = v_supng  (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_supng(1,i)
          end do
        outrcd( 5) = z / dble(lastcell - firstcell + 1  ) ! - av. cell somata 
         z = 0.d0
          do i = 1, numcomp_supng   
           z = z + gAMPA_supng  (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_supng   
           z = z + gNMDA_supng  (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_supng  
           z = z + gGABA_A_supng  (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_supng  (1,firstcell+2)
        outrcd(10) = v_supng  (1,firstcell+3)
        outrcd(11) = field_1mm_tot      
        outrcd(12) = field_2mm_tot      
         if (place(thisno).eq.1) then
      OPEN(33,FILE='alphaY67.supng')
c     WRITE (33,FMT='(12F10.4)') (OUTRCD(I),I=1,12)
         endif

c      else if (thisno.eq.3) then
       else IF (nodecell(thisno) .eq. 'supaxax  ') THEN
c supaxax 
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_supaxax 
          lastcell = firstcell - 1 + ncellspernode_supaxax 

        outrcd( 1) = time
        outrcd( 2) = v_supaxax  (1,firstcell+1)
        outrcd( 3) = v_supaxax  (numcomp_supaxax  ,firstcell+1)
        outrcd( 4) = v_supaxax  (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_supaxax(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_supaxax  
           z = z + gAMPA_supaxax  (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_supaxax  
           z = z + gNMDA_supaxax  (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_supaxax  
           z = z + gGABA_A_supaxax  (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_supaxax  (1,firstcell+2)
        outrcd(10) = v_supaxax  (1,firstcell+3)
        outrcd(11) = field_1mm_tot
        outrcd(12) = field_2mm_tot
          if (place(thisno).eq.1) then
      OPEN(14,FILE='alphaY67.supaxax')
      WRITE (14,FMT='(12F10.4)') (OUTRCD(I),I=1,12)
          endif

c      else if (thisno.eq.4) then
       else IF (nodecell(thisno) .eq. 'supLTS   ') THEN
c supLTS  
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_supLTS  
          lastcell = firstcell - 1 + ncellspernode_supLTS  

        outrcd( 1) = time
        outrcd( 2) = v_supLTS   (1,firstcell+1)
        outrcd( 3) = v_supLTS   (numcomp_supLTS   ,firstcell+1)
        outrcd( 4) = v_supLTS   (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_supLTS(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_supLTS   
           z = z + gAMPA_supLTS   (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_supLTS   
           z = z + gNMDA_supLTS   (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_supLTS   
           z = z + gGABA_A_supLTS   (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_supLTS   (1,firstcell+2)
        outrcd(10) = v_supLTS   (1,firstcell+3)
        outrcd(11) = field_1mm_tot
        outrcd(12) = field_2mm_tot
         if (place(thisno).eq.1) then
      OPEN(15,FILE='alphaY67.supLTS')
      WRITE (15,FMT='(12F10.4)') (OUTRCD(I),I=1,12)
         endif

c      else if (thisno.eq.5) then
       else IF (nodecell(thisno) .eq. 'spinstell') THEN
c spinstell
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_spinstell
          lastcell = firstcell - 1 + ncellspernode_spinstell

        outrcd( 1) = time
        outrcd( 2) = v_spinstell(1,firstcell+1)
        outrcd( 3) = v_spinstell(numcomp_spinstell,firstcell+1)
        outrcd( 4) = v_spinstell(43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_spinstell(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_spinstell
           z = z + gAMPA_spinstell(i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 1 
         z = 0.d0
          do i = 1, numcomp_spinstell
           z = z + gNMDA_spinstell(i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 1 
         z = 0.d0
          do i = 1, numcomp_spinstell
           z = z + gGABA_A_spinstell(i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 1 
         z = 0.d0
          do i = 1, numcomp_spinstell
           z = z + gGABA_B_spinstell(i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 1 
        outrcd(10) = v_spinstell(1,firstcell+2)
        outrcd(11) = v_spinstell(1,firstcell+3)
        outrcd(12) = v_spinstell(1,firstcell+4)
        outrcd(13) = v_spinstell(1,firstcell+5)

c       outrcd(12) = field_1mm_tot    
c       outrcd(13) = field_2mm_tot

          if (place(thisno).eq.1) then
      OPEN(16,FILE='alphaY67.spinstell')
      WRITE (16,FMT='(13F10.4)') (OUTRCD(I),I=1,13)
          endif

c      do i = 1, num_spinstell
c       if ((time.ge.300.d0).and.(v_spinstell(1,i).gt.10.d0)) then
c          OPEN(30,FILE='alphaY67.spinstell.raster')
c          WRITE(30,FMT='(2F10.4)') time, dble(i)
c       endif
c      end do

c      else if (thisno.eq.6) then
       else IF (nodecell(thisno) .eq. 'tuftIB   ') THEN
c tuftIB  
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_tuftIB  
          lastcell = firstcell - 1 + ncellspernode_tuftIB  

        outrcd( 1) = time
        outrcd( 2) = v_tuftIB   (1,firstcell+1)
        outrcd( 3) = v_tuftIB   (numcomp_tuftIB   ,firstcell+1)
        outrcd( 4) = v_tuftIB   (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_tuftIB(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_tuftIB   
           z = z + gAMPA_tuftIB   (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_tuftIB   
           z = z + gNMDA_tuftIB   (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_tuftIB   
           z = z + gGABA_A_tuftIB   (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_tuftIB   
           z = z + gGABA_B_tuftIB   (i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_tuftIB   (1,firstcell+2)
        outrcd(11) = v_tuftIB   (1,firstcell+3)
        outrcd(12) = field_1mm_tot   
        outrcd(13) = field_2mm_tot   
        outrcd(14) = v_tuftIB   (1,firstcell+4)
        outrcd(15) = v_tuftIB   (numcomp_tuftIB   ,firstcell+4)
          z = 0.d0
        do L = 1, num_tuftIB
          if (ldistal_axon_tuftIB(L).ge.0.d0) z = z + 1.d0
        end do
        outrcd(16) = z ! should be number of distal tuftIB axons overshooting
        outrcd(17) = v_tuftIB (1,firstcell + 5)
        outrcd(18) = v_tuftIB (1,firstcell + 6)
        outrcd(19) = v_tuftIB (1,firstcell + 7)
c       outrcd(20) = v_tuftIB (1,firstcell + 8)
        outrcd(20) = 1000.d0 * noisepe_tuftIB 
          if (place(thisno).eq.1) then
c     OPEN(17,FILE='alphaY67.tuftIB')
c     WRITE (17,FMT='(20F10.4)') (OUTRCD(I),I=1,20)
          else
c           write(6,9091) 'tuftIB', thisno, time, v_tuftIB(1,firstcell),
c    &            v_tuftIB(1,lastcell)
9091        format(a6,i4,3f10.4)
          endif

c      else if (thisno.eq.7) then
       else IF (nodecell(thisno) .eq. 'tuftRS   ') THEN
c tuftRS  
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_tuftRS  
          lastcell = firstcell - 1 + ncellspernode_tuftRS  

        outrcd( 1) = time
        outrcd( 2) = v_tuftRS   (1,firstcell+1)
        outrcd( 3) = v_tuftRS   (numcomp_tuftRS   ,firstcell+1)
        outrcd( 4) = v_tuftRS   (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_tuftRS(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_tuftRS   
           z = z + gAMPA_tuftRS   (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_tuftRS   
           z = z + gNMDA_tuftRS   (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_tuftRS   
           z = z + gGABA_A_tuftRS   (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_tuftRS   
           z = z + gGABA_B_tuftRS   (i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_tuftRS   (1,firstcell+2)
        outrcd(11) = v_tuftRS   (1,firstcell+3)
        outrcd(12) = field_1mm_tot   
        outrcd(13) = field_2mm_tot    
          if (place(thisno).eq.1) then
c     OPEN(18,FILE='alphaY67.tuftRS')
c     WRITE (18,FMT='(13F10.4)') (OUTRCD(I),I=1,13)
          end if

c      else if (thisno.eq.8) then
       else IF (nodecell(thisno) .eq. 'nontuftRS') THEN
c nontuftRS
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_nontuftRS
          lastcell = firstcell - 1 + ncellspernode_nontuftRS

        outrcd( 1) = time
        outrcd( 2) = v_nontuftRS(1,firstcell+1)
        outrcd( 3) = v_nontuftRS(numcomp_nontuftRS,firstcell+1)
        outrcd( 4) = v_nontuftRS(43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_nontuftRS(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_nontuftRS
           z = z + gAMPA_nontuftRS(i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_nontuftRS
           z = z + gNMDA_nontuftRS(i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_nontuftRS
           z = z + gGABA_A_nontuftRS(i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_nontuftRS
           z = z + gGABA_B_nontuftRS(i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_nontuftRS(1,firstcell+2)
        outrcd(11) = v_nontuftRS(1,firstcell+3)
        outrcd(12) = field_1mm_tot      
        outrcd(13) = field_2mm_tot       
        outrcd(14) = v_nontuftRS(1,lastcell     -2)
        outrcd(15) = v_nontuftRS(1,lastcell     -1)
        outrcd(16) = v_nontuftRS(1,lastcell     )
          if (place(thisno).eq.1) then
c     OPEN(19,FILE='alphaY67.nontuftRS')
c     WRITE (19,FMT='(16F10.4)') (OUTRCD(I),I=1,16)
          else
c      write(6,9092) 'nontuftRS',thisno,time,v_nontuftRS(1,firstcell),
c    &            v_nontuftRS(1,lastcell)
9092        format(a9,i4,3f10.4)
          endif

c      else if (thisno.eq.9) then
       else IF (nodecell(thisno) .eq. 'deepbask ') THEN
c deepbask 
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_deepbask
          lastcell = firstcell - 1 + ncellspernode_deepbask

        outrcd( 1) = time
        outrcd( 2) = v_deepbask (1,firstcell+1)
        outrcd( 3) = v_deepbask (numcomp_deepbask ,firstcell+1)
        outrcd( 4) = v_deepbask (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_deepbask(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_deepbask 
           z = z + gAMPA_deepbask (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepbask 
           z = z + gNMDA_deepbask (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepbask 
           z = z + gGABA_A_deepbask (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_deepbask (1,firstcell+2)
        outrcd(10) = v_deepbask (1,firstcell+3)
        outrcd(11) = field_1mm_tot
        outrcd(12) = field_2mm_tot
           if (place(thisno).eq.1) then
c     OPEN(20,FILE='alphaY67.deepbask')
c     WRITE (20,FMT='(12F10.4)') (OUTRCD(I),I=1,12)
           endif

c      do i = 1, num_deepbask 
c       if ((time.ge.300.d0).and.(v_deepbask (1,i).gt.10.d0)) then
c          OPEN(31,FILE='alphaY67.deepbask.raster')
c          WRITE(31,FMT='(2F10.4)') time, dble(i)
c       endif
c      end do

       else IF (nodecell(thisno) .eq. 'deepng ') THEN
c deepng 
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_deepng
          lastcell = firstcell - 1 + ncellspernode_deepng

        outrcd( 1) = time
        outrcd( 2) = v_deepng (1,firstcell+1)
        outrcd( 3) = v_deepng (numcomp_deepng ,firstcell+1)
        outrcd( 4) = v_deepng (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_deepng(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_deepng 
           z = z + gAMPA_deepng (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepng 
           z = z + gNMDA_deepng (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepng 
           z = z + gGABA_A_deepng (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_deepng (1,firstcell+2)
        outrcd(10) = v_deepng (1,firstcell+3)
        outrcd(11) = field_1mm_tot
        outrcd(12) = field_2mm_tot
           if (place(thisno).eq.1) then
c     OPEN(34,FILE='alphaY67.deepng')
c     WRITE (34,FMT='(12F10.4)') (OUTRCD(I),I=1,12)
           endif


c      else if (thisno.eq.10) then
       else IF (nodecell(thisno) .eq. 'deepaxax ') THEN
c deepaxax
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_deepaxax
          lastcell = firstcell - 1 + ncellspernode_deepaxax

        outrcd( 1) = time
        outrcd( 2) = v_deepaxax (1,firstcell+1)
        outrcd( 3) = v_deepaxax (numcomp_deepaxax ,firstcell+1)
        outrcd( 4) = v_deepaxax (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_deepaxax(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_deepaxax 
           z = z + gAMPA_deepaxax (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepaxax 
           z = z + gNMDA_deepaxax (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepaxax 
           z = z + gGABA_A_deepaxax (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_deepaxax (1,firstcell+2)
        outrcd(10) = v_deepaxax (1,firstcell+3)
        outrcd(11) = field_1mm_tot
        outrcd(12) = field_2mm_tot
           if (place(thisno).eq.1) then
c     OPEN(21,FILE='alphaY67.deepaxax')
c     WRITE (21,FMT='(12F10.4)') (OUTRCD(I),I=1,12)
           endif

c      else if (thisno.eq.11) then
       else IF (nodecell(thisno) .eq. 'deepLTS  ') THEN
c deepLTS 
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_deepLTS 
          lastcell = firstcell - 1 + ncellspernode_deepLTS 

        outrcd( 1) = time
        outrcd( 2) = v_deepLTS  (1,firstcell+1)
        outrcd( 3) = v_deepLTS  (numcomp_deepLTS  ,firstcell+1)
        outrcd( 4) = v_deepLTS  (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_deepLTS(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_deepLTS  
           z = z + gAMPA_deepLTS  (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepLTS   
           z = z + gNMDA_deepLTS  (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepLTS  
           z = z + gGABA_A_deepLTS  (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_deepLTS  (1,firstcell+2)
        outrcd(10) = v_deepLTS  (1,firstcell+3)
        outrcd(11) = v_deepLTS  (1,firstcell+4)
        outrcd(12) = v_deepLTS  (1,firstcell+5)
c       outrcd(11) = field_1mm_tot
c       outrcd(12) = field_2mm_tot
         if (place(thisno).eq.1) then
      OPEN(22,FILE='alphaY67.deepLTS')
      WRITE (22,FMT='(12F10.4)') (OUTRCD(I),I=1,12)
         endif

c      else if (thisno.eq.12) then
       else IF (nodecell(thisno) .eq. 'TCR      ') THEN
c TCR      
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_TCR     
          lastcell = firstcell - 1 + ncellspernode_TCR      

        outrcd( 1) = time
        outrcd( 2) = v_TCR      (1,firstcell+1)
        outrcd( 3) = v_TCR      (numcomp_TCR ,firstcell+1)
        outrcd( 4) = v_TCR      (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_TCR(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_TCR      
           z = z + gAMPA_TCR      (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_TCR      
           z = z + gNMDA_TCR      (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_TCR      
           z = z + gGABA_A_TCR      (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_TCR      
           z = z + gGABA_B_TCR      (i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_TCR      (1,firstcell+2)
        outrcd(11) = v_TCR      (1,firstcell+3)

          z = 0.d0
          do i = firstcell, lastcell
           if(v_TCR     (numcomp_TCR     ,i) .gt. 0.d0) z = z + 1.d0
          end do
        outrcd(12) = z   

        outrcd(13) = v_TCR (1,33)
        outrcd(14) = field_1mm_tot
        outrcd(15) = field_2mm_tot

         if (place(thisno).eq.1) then
c     OPEN(23,FILE='alphaY67.TCR')
c     WRITE (23,FMT='(15F10.4)') (OUTRCD(I),I=1,15)
         endif

c      else if (thisno.eq.13) then
       else IF (nodecell(thisno) .eq. 'nRT      ') THEN
c nRT       
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode_nRT        
          lastcell = firstcell - 1 + ncellspernode_nRT       

        outrcd( 1) = time
        outrcd( 2) = v_nRT      (1,firstcell+1)
        outrcd( 3) = v_nRT      (numcomp_nRT ,firstcell+1)
        outrcd( 4) = v_nRT      (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_nRT(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_nRT       
           z = z + gAMPA_nRT      (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_nRT         
           z = z + gNMDA_nRT      (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_nRT        
           z = z + gGABA_A_nRT      (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_nRT        
           z = z + gGABA_B_nRT      (i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_nRT      (1,firstcell+2)
        outrcd(11) = v_nRT      (1,firstcell+3)

          z = 0.d0
          do i = firstcell, lastcell
           if(v_nRT     (numcomp_nRT     ,i) .gt. 0.d0) z = z + 1.d0
          end do
        outrcd(12) = z   
        outrcd(13) = field_1mm_tot
        outrcd(14) = field_2mm_tot
         if (place(thisno).eq.1) then
c     OPEN(24,FILE='alphaY67.nRT')
c     WRITE (24,FMT='(14F10.4)') (OUTRCD(I),I=1,14)
         endif
       endif ! checking thisno

       endif ! mod(O, 50) = 0

        goto 1000
c END guts of main program

2000    CONTINUE
        time2 = gettime()
         if (thisno.eq.0) then
        write(6,3434) time2 - time1
         endif
3434    format(' elapsed time = ',f8.0,' secs')

        call mpi_finalize (info)
             END
c end main routine


c 11 Sept 2006, start with /interact/integrate_suppyrRSXP.f & add GABA-B
! 7 Nov. 2005: modify integrate_suppyrRSX.f to allow for Colbert-Pan axon.
!29 July 2005: modify groucho/integrate_suppyrRS.f, for a separate
! call for initialization, and to integrate only selected cells.
! Integration routine for suppyrRS cells
! Routine adapted from scortn in supergj.f
       SUBROUTINE INTEGRATE_suppyrRSXPB (O, time, numcell,     
     &    V, curr, initialize, firstcell, lastcell,
     & gAMPA, gNMDA, gGABA_A, gGABA_B,
     & Mg, 
     & gapcon  ,totaxgj   ,gjtable, dt,
     & totaxgj_mix, gjtable_mix, num_other, vax,
     &  chi,mnaf,mnap,
     &  hnaf,mkdr,mka,
     &  hka,mk2,hk2,
     &  mkm,mkc,mkahp,
     &  mcat,hcat,mcal,
     &  mar,field_1mm,field_2mm,rel_axonshift)
! num_other here = num suppyrFRB, vax = vax_suppyrFRB

       SAVE

       INTEGER, PARAMETER:: numcomp = 74
! numcomp here must be compatible with numcomp_suppyrRS in calling prog.
       INTEGER  numcell, num_other
       INTEGER initialize, firstcell, lastcell
       INTEGER J1, I, J, K, K1, K2, K3, L, L1, O
       REAL*8 c(numcomp), curr(numcomp,numcell)
       REAL*8  Z, Z1, Z2, Z3, Z4, Z5, Z6, Z7, DT, time
       REAL*8 Z8, Z9, Z7a, Z4a, Z1a
       integer totaxgj, gjtable(totaxgj,4)
       integer totaxgj_mix, gjtable_mix(totaxgj_mix,4)
       real*8 gapcon, gAMPA(numcomp,numcell),
     &        gNMDA(numcomp,numcell), gGABA_A(numcomp,numcell),
     &        gGABA_B(numcomp,numcell)
       real*8 Mg, V(numcomp,numcell)
!      real*8 vax (num_other), rel_axonshift
       real*8 vax (numcomp,num_other) ! "vax" now is voltages of another cell type,
! assumed to have same number of compartments as this cell type

c CINV is 1/C, i.e. inverse capacitance
       real*8 chi(numcomp,numcell),
     & mnaf(numcomp,numcell),mnap(numcomp,numcell),
     x hnaf(numcomp,numcell), mkdr(numcomp,numcell),
     x mka(numcomp,numcell),hka(numcomp,numcell),
     x mk2(numcomp,numcell), cinv(numcomp),
     x hk2(numcomp,numcell),mkm(numcomp,numcell),
     x mkc(numcomp,numcell),mkahp(numcomp,numcell),
     x mcat(numcomp,numcell),hcat(numcomp,numcell),
     x mcal(numcomp,numcell), betchi(numcomp),
     x mar(numcomp,numcell),jacob(numcomp,numcomp),
     x gam(0: numcomp,0: numcomp),gL(numcomp),gnaf(numcomp),
     x gnap(numcomp),gkdr(numcomp),gka(numcomp),
     x gk2(numcomp),gkm(numcomp),
     x gkc(numcomp),gkahp(numcomp),
     x gcat(numcomp),gcaL(numcomp),gar(numcomp),
     x cafor(numcomp)
       real*8
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
       real*8 vL(numcomp),vk(numcomp),vna,var,vca,vgaba_a
       real*8 depth(12), membcurr(12), field_1mm, field_2mm
       integer level(numcomp)

        INTEGER NEIGH(numcomp,10), NNUM(numcomp)
        INTEGER igap1, igap2
c the f's are the functions giving 1st derivatives for evolution of
c the differential equations for the voltages (v), calcium (chi), and
c other state variables.
       real*8 fv(numcomp), fchi(numcomp),
     x fmnaf(numcomp),fhnaf(numcomp),fmkdr(numcomp),
     x fmka(numcomp),fhka(numcomp),fmk2(numcomp),
     x fhk2(numcomp),fmnap(numcomp),
     x fmkm(numcomp),fmkc(numcomp),fmkahp(numcomp),
     x fmcat(numcomp),fhcat(numcomp),fmcal(numcomp),
     x fmar(numcomp)

c below are for calculating the partial derivatives
       real*8 dfv_dv(numcomp,numcomp), dfv_dchi(numcomp),
     x  dfv_dmnaf(numcomp),  dfv_dmnap(numcomp),
     x  dfv_dhnaf(numcomp),dfv_dmkdr(numcomp),
     x  dfv_dmka(numcomp),dfv_dhka(numcomp),
     x  dfv_dmk2(numcomp),dfv_dhk2(numcomp),
     x  dfv_dmkm(numcomp),dfv_dmkc(numcomp),
     x  dfv_dmkahp(numcomp),dfv_dmcat(numcomp),
     x  dfv_dhcat(numcomp),dfv_dmcal(numcomp),
     x  dfv_dmar(numcomp)

        real*8 dfchi_dv(numcomp), dfchi_dchi(numcomp),
     x dfmnaf_dmnaf(numcomp), dfmnaf_dv(numcomp),
     x dfhnaf_dhnaf(numcomp),
     x dfmnap_dmnap(numcomp), dfmnap_dv(numcomp),
     x dfhnaf_dv(numcomp),dfmkdr_dmkdr(numcomp),
     x dfmkdr_dv(numcomp),
     x dfmka_dmka(numcomp),dfmka_dv(numcomp),
     x dfhka_dhka(numcomp),dfhka_dv(numcomp),
     x dfmk2_dmk2(numcomp),dfmk2_dv(numcomp),
     x dfhk2_dhk2(numcomp),dfhk2_dv(numcomp),
     x dfmkm_dmkm(numcomp),dfmkm_dv(numcomp),
     x dfmkc_dmkc(numcomp),dfmkc_dv(numcomp),
     x dfmcat_dmcat(numcomp),dfmcat_dv(numcomp),dfhcat_dhcat(numcomp),
     x dfhcat_dv(numcomp),dfmcal_dmcal(numcomp),dfmcal_dv(numcomp),
     x dfmar_dmar(numcomp),dfmar_dv(numcomp),dfmkahp_dchi(numcomp),
     x dfmkahp_dmkahp(numcomp), dt2

       REAL*8 OPEN(numcomp),gamma(numcomp),gamma_prime(numcomp)
c gamma is function of chi used in calculating KC conductance
       REAL*8 alpham_ahp(numcomp), alpham_ahp_prime(numcomp)
       REAL*8 gna_tot(numcomp),gk_tot(numcomp),gca_tot(numcomp)
       REAL*8 gca_high(numcomp), gar_tot(numcomp)
c this will be gCa conductance corresponding to high-thresh channels

       real*8 persistentNa_shift, fastNa_shift_SD,
     x   fastNa_shift_axon

       REAL*8 A, BB1, BB2  ! params. for FNMDA.f


c          if (O.eq.1) then
           if (initialize.eq.0) then
c do initialization

c Program fnmda assumes A, BB1, BB2 defined in calling program
c as follows:
         A = DEXP(-2.847d0)
         BB1 = DEXP(-.693d0)
         BB2 = DEXP(-3.101d0)

c       goto 4000
       CALL   SCORT_SETUP_suppyrRS
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)

        CALL SCORTMAJ_suppyrRS
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM,depth,level)

          do i = 1, numcomp
             cinv(i) = 1.d0 / c(i)
          end do
4000      CONTINUE

           do i = 1, numcomp
          vL(i) = -70.d0
          vK(i) = -95.d0
           end do

        VNA = 50.d0
        VCA = 125.d0
        VAR = -43.d0
        VAR = -35.d0
c -43 mV from Huguenard & McCormick
        VGABA_A = -81.d0
c       write(6,901) VNa, VCa, VK(1), O
901     format('VNa =',f6.2,' VCa =',f6.2,' VK =',f6.2,
     &   ' O = ',i3)

c ? initialize membrane state variables?
         do L = 1, numcell  
         do i = 1, numcomp
        v(i,L) = VL(i)
	chi(i,L) = 0.d0
	mnaf(i,L) = 0.d0
	mkdr(i,L) = 0.d0
	mk2(i,L) = 0.d0
	mkm(i,L) = 0.d0
	mkc(i,L) = 0.d0
	mkahp(i,L) = 0.d0
	mcat(i,L) = 0.d0
	mcal(i,L) = 0.d0
         end do
         end do

          do L = 1, numcell
        k1 = idnint (4.d0 * (v(1,L) + 120.d0))

            do i = 1, numcomp
      hnaf(i,L) = alphah_naf(k1)/(alphah_naf(k1)
     &       +betah_naf(k1))
      hka(i,L) = alphah_ka(k1)/(alphah_ka(k1)
     &                               +betah_ka(k1))
      hk2(i,L) = alphah_k2(k1)/(alphah_k2(k1)
     &                                +betah_k2(k1))
      hcat(i,L)=alphah_cat(k1)/(alphah_cat(k1)
     &                                +betah_cat(k1))
c     mar=alpham_ar(k1)/(alpham_ar(k1)+betam_ar(k1))
      mar(i,L) = .25d0
             end do
           end do


             do i = 1, numcomp
	    open(i) = 0.d0
            gkm(i) = 2.d0 * gkm(i)
             end do

         do i = 1, 68
c          gnaf(i) = 0.8d0 * 1.25d0 * gnaf(i) ! factor of 0.8 added 19 Nov. 2005
c          gnaf(i) = 0.9d0 * 1.25d0 * gnaf(i) ! Back to 0.9, 29 Nov. 2005
           gnaf(i) = 0.6d0 * 1.25d0 * gnaf(i) ! 
! NOTE THAT THERE IS QUESTION OF HOW TO COMPARE BEHAVIOR OF PYRAMID IN NETWORK WITH
! SIMULATIONS OF SINGLE CELL.  IN FORMER CASE, THERE IS LARGE AXONAL SHUNT THROUGH
! gj(s), NOT PRESENT IN SINGLE CELL MODEL.  THEREFORE, HIGHER AXONAL gNa MIGHT BE
! NECESSARY FOR SPIKE PROPAGATION.
c          gnaf(i) = 0.9d0 * 1.25d0 * gnaf(i) ! factor of 0.9 added 20 Nov. 2005
           gkdr(i) = 1.25d0 * gkdr(i)
         end do
 
c Perhaps reduce fast gNa on IS
          gnaf(69) = 1.00d0 * gnaf(69)
c         gnaf(69) = 0.25d0 * gnaf(69)
          gnaf(70) = 1.00d0 * gnaf(70)
c         gnaf(70) = 0.25d0 * gnaf(70)

c Perhaps reduce coupling between soma and IS
c         gam(1,69) = 0.15d0 * gam(1,69)
c         gam(69,1) = 0.15d0 * gam(69,1)

c Determine which cells are FRB and which RS.
c           if (L.le.(nL2-nFRB)) then
c              z1 = 0.60d0 ! for gNaP, axon
               z1 = 0.00d0 ! for gNaP, axon
               z1a= 0.4d0  ! for soma, dendrites 
c              z1a= 0.3d0  ! for soma, dendrites 
c              z1a= 0.5d0  ! for soma, dendrites 
c              z2 = 2.6d0
c              z2 = 2.1d0
!              z2 = 2.3d0
               z2 = 2.00d0
c              z2 = 1.2d0 ! value 1.2 tried Feb. 21, 2013
c              z3 = 1.0d0
               z3 = 0.80d0
c              z3 = 0.0d0 ! Note reduction from 0.4, to prevent
c slow hyperpolarization that seems to mess up gamma.
c              z4 = 1.0d0
c              z4 = 0.250d0
               z4 = 0.65d0 ! for axon
!              z4 = 0.80d0 ! for axon
               z4a= 0.35d0 ! for soma/dend
c              z5 = 0.35d0
               z5 = 0.30d0
c              z6 =  5.d0 ! NOTE possibly increased gCaT & other changes - see scortalphaY674
               z6 =  1.d0 ! NOTE possibly increased gCaT & other changes - see scortalphaY674
c              z6 = 10.d0 ! NOTE possibly increased gCaT & other changes - see scortalphaY674
               z7 = 2.d0 ! axon gKA
c              z7 = 4.d0 ! axon gKA
               z7a= 2.d0 ! SD gKA
c              z7 = 0.d0
               z8 = 1.d0
               z9 = 0.d0
c              z9 = 1.d0

c RS cell
c           else
c              z1 = 1.50d0
c              z2 = 0.6d0
c              z3 = 1.d0
c              z4 = 1.d0
c FRB cell
c           endif
             do i = 1, numcomp
           if (i.ge.69) then ! axon  
              gnap(i) = z1 * gnap(i)
           else  ! soma/dendrites
              gnap(i) = z1a* gnap(i)
           endif
              gkc (i) = z2 * gkc (i)
              gkahp(i) = z3 * gkahp(i)
         if (i.ge.69) then ! Axon
              gkm (i) = z4 * gkm (i)
         else
              gkm (i) = z4a * gkm (i) ! soma/dendrites
         endif
              gCaL(i) = z5 * gCaL(i)
              gCaT(i) = z6 * gCaT(i)
              gK2(i)  = z8 * gK2(i)
              gAR(i)  = z9 * gAR(i) 
             end do

            do i = 69, 74
             gKA(i) = z7 * gKA(i)
            end do
            do i = 1, 68
             gKA(i) = z7a* gKA(i)
            end do

              goto 6000

          endif
c End initialization

          do i = 1, 12
           membcurr(i) = 0.d0
          end do

c                  goto 2001


c             do L = 1, numcell
              do L = firstcell, lastcell

	  do i = 1, numcomp
	  do j = 1, nnum(i)
	   if (neigh(i,j).gt.numcomp) then
          write(6,433) i, j, L
433       format(' ls ',3x,3i5)
           endif
	end do
	end do

       DO I = 1, numcomp
          FV(I) = -GL(I) * (V(I,L) - VL(i)) * cinv(i)
          DO J = 1, NNUM(I)
             K = NEIGH(I,J)
302     FV(I) = FV(I) + GAM(I,K) * (V(K,L) - V(I,L)) * cinv(i)
           END DO
       END DO
301    CONTINUE


       CALL FNMDA (V, OPEN, numcell, numcomp, MG, L,
     &                 A, BB1, BB2)

      DO I = 1, numcomp
       FV(I) = FV(I) + ( CURR(I,L)
     X   - (gampa(I,L) + open(i) * gnmda(I,L))*V(I,L)
     X   - ggaba_a(I,L)*(V(I,L)-Vgaba_a) 
     X   - ggaba_b(I,L)*(V(I,L)-VK(i)  ) ) * cinv(i)
c above assumes equil. potential for AMPA & NMDA = 0 mV
      END DO
421      continue

       do m = 1, totaxgj
        if (gjtable(m,1).eq.L) then
         L1 = gjtable(m,3)
         igap1 = gjtable(m,2)
         igap2 = gjtable(m,4)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        else if (gjtable(m,3).eq.L) then
         L1 = gjtable(m,1)
         igap1 = gjtable(m,4)
         igap2 = gjtable(m,2)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        endif
       end do ! do m

       do m = 1, totaxgj_mix
        if (gjtable_mix(m,1).eq.L) then
         L1 = gjtable_mix(m,3)
         igap1 = gjtable_mix(m,2)
 	fv(igap1) = fv(igap1) + gapcon *
c    &   (vax(L1) - v(igap1,L)) * cinv(igap1)
     &   (vax(igap1,L1) - v(igap1,L)) * cinv(igap1)
c Note that compartments assumed to match
        endif
       end do ! do m

c      do i = 1, ngap(L) ! OBSOLETE CODE FROM SUPERGJ.F
c      L1 = list_gap(L,i)
c	fv(axsite) = fv(axsite) + gapcon *
c    &   (vaxgap_global(L1) - v(axsite,L)) * cinv(axsite)
c      end do  ! OBSOLETE CODE FROM SUPERGJ.F

       do i = 1, numcomp
        gamma(i) = dmin1 (1.d0, .004d0 * chi(i,L))
        if (chi(i,L).le.250.d0) then
          gamma_prime(i) = .004d0
        else
          gamma_prime(i) = 0.d0
        endif
c         endif
       end do

      DO I = 1, numcomp
       gna_tot(i) = gnaf(i) * (mnaf(i,L)**3) * hnaf(i,L) +
     x     gnap(i) * mnap(i,L)
       gk_tot(i) = gkdr(i) * (mkdr(i,L)**4) +
     x             gka(i)  * (mka(i,L)**4) * hka(i,L) +
     x             gk2(i)  * mk2(i,L) * hk2(i,L) +
     x             gkm(i)  * mkm(i,L) +
     x             gkc(i)  * mkc(i,L) * gamma(i) +
     x             gkahp(i)* mkahp(i,L)
       gca_tot(i) = gcat(i) * (mcat(i,L)**2) * hcat(i,L) +
     x              gcaL(i) * (mcaL(i,L)**2)
       gca_high(i) =
     x              gcaL(i) * (mcaL(i,L)**2)
       gar_tot(i) = gar(i) * mar(i,L)


       FV(I) = FV(I) - ( gna_tot(i) * (v(i,L) - vna)
     X  + gk_tot(i) * (v(i,L) - vK(i))
     X  + gca_tot(i) * (v(i,L) - vCa)
     X  + gar_tot(i) * (v(i,L) - var) ) * cinv(i)
c        endif
       END DO
88           continue

         do i = 1, numcomp
         do j = 1, numcomp
          if (i.ne.j) then
            dfv_dv(i,j) = jacob(i,j)
          else
            dfv_dv(i,j) = jacob(i,i) - cinv(i) *
     X  (gna_tot(i) + gk_tot(i) + gca_tot(i) + gar_tot(i)
     X   + ggaba_a(i,L) + ggaba_b(i,L) + gampa(i,L)
     X   + open(i) * gnmda(I,L) )
          endif
         end do
         end do

           do i = 1, numcomp
        dfv_dchi(i)  = - cinv(i) * gkc(i) * mkc(i,L) *
     x                     gamma_prime(i) * (v(i,L)-vK(i))
        dfv_dmnaf(i) = -3.d0 * cinv(i) * (mnaf(i,L)**2) *
     X    (gnaf(i) * hnaf(i,L)          ) * (v(i,L) - vna)
        dfv_dmnap(i) = - cinv(i) *
     X    (               gnap(i)) * (v(i,L) - vna)
        dfv_dhnaf(i) = - cinv(i) * gnaf(i) * (mnaf(i,L)**3) *
     X                    (v(i,L) - vna)
        dfv_dmkdr(i) = -4.d0 * cinv(i) * gkdr(i) * (mkdr(i,L)**3)
     X                   * (v(i,L) - vK(i))
        dfv_dmka(i)  = -4.d0 * cinv(i) * gka(i) * (mka(i,L)**3) *
     X                   hka(i,L) * (v(i,L) - vK(i))
        dfv_dhka(i)  = - cinv(i) * gka(i) * (mka(i,L)**4) *
     X                    (v(i,L) - vK(i))
      dfv_dmk2(i) = - cinv(i) * gk2(i) * hk2(i,L) * (v(i,L)-vK(i))
      dfv_dhk2(i) = - cinv(i) * gk2(i) * mk2(i,L) * (v(i,L)-vK(i))
      dfv_dmkm(i) = - cinv(i) * gkm(i) * (v(i,L) - vK(i))
      dfv_dmkc(i) = - cinv(i)*gkc(i) * gamma(i) * (v(i,L)-vK(i))
        dfv_dmkahp(i)= - cinv(i) * gkahp(i) * (v(i,L) - vK(i))
        dfv_dmcat(i)  = -2.d0 * cinv(i) * gcat(i) * mcat(i,L) *
     X                    hcat(i,L) * (v(i,L) - vCa)
        dfv_dhcat(i) = - cinv(i) * gcat(i) * (mcat(i,L)**2) *
     X                  (v(i,L) - vCa)
        dfv_dmcal(i) = -2.d0 * cinv(i) * gcal(i) * mcal(i,L) *
     X                      (v(i,L) - vCa)
        dfv_dmar(i) = - cinv(i) * gar(i) * (v(i,L) - var)
            end do

         do i = 1, numcomp
          fchi(i) = - cafor(i) * gca_high(i) * (v(i,L) - vca)
     x       - betchi(i) * chi(i,L)
          dfchi_dv(i) = - cafor(i) * gca_high(i)
          dfchi_dchi(i) = - betchi(i)
         end do

       do i = 1, numcomp
c Note possible increase in rate at which AHP current develops
c       alpham_ahp(i) = dmin1(0.2d-4 * chi(i,L),0.01d0)
        alpham_ahp(i) = dmin1(1.0d-4 * chi(i,L),0.01d0)
        if (chi(i,L).le.500.d0) then
c         alpham_ahp_prime(i) = 0.2d-4
          alpham_ahp_prime(i) = 1.0d-4
        else
          alpham_ahp_prime(i) = 0.d0
        endif
       end do

       do i = 1, numcomp
        fmkahp(i) = alpham_ahp(i) * (1.d0 - mkahp(i,L))
c    x                  -.001d0 * mkahp(i,L)
c    x                  -.010d0 * mkahp(i,L)
     x                  -.007d0 * mkahp(i,L)
c       dfmkahp_dmkahp(i) = - alpham_ahp(i) - .001d0
        dfmkahp_dmkahp(i) = - alpham_ahp(i) - .010d0
        dfmkahp_dchi(i) = alpham_ahp_prime(i) *
     x                     (1.d0 - mkahp(i,L))
       end do

          do i = 1, numcomp

       K1 = IDNINT ( 4.d0 * (V(I,L) + 120.d0) )
       IF (K1.GT.640) K1 = 640
       IF (K1.LT.  0) K1 =   0

c      persistentNa_shift =  0.d0
c      persistentNa_shift =  8.d0
       persistentNa_shift = 10.d0
       K2 = IDNINT ( 4.d0 * (V(I,L)+persistentNa_shift+ 120.d0) )
       IF (K2.GT.640) K2 = 640
       IF (K2.LT.  0) K2 =   0

c            fastNa_shift = -2.0d0
c            fastNa_shift = -2.5d0
             fastNa_shift_SD = -3.5d0
             fastNa_shift_axon = fastNa_shift_SD + rel_axonshift 
       K0 = IDNINT ( 4.d0 * (V(I,L)+  fastNa_shift_SD+ 120.d0) )
       IF (K0.GT.640) K0 = 640
       IF (K0.LT.  0) K0 =   0
       K3 = IDNINT ( 4.d0 * (V(I,L)+  fastNa_shift_axon+ 120.d0) )
       IF (K3.GT.640) K3 = 640
       IF (K3.LT.  0) K3 =   0

         if (i.le.68) then   ! FOR SD
        fmnaf(i) = alpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k0) * mnaf(i,L)
        fhnaf(i) = alphah_naf(k0) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k0) * hnaf(i,L)
         else  ! for axon
        fmnaf(i) = alpham_naf(k3) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k3) * mnaf(i,L)
        fhnaf(i) = alphah_naf(k3) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k3) * hnaf(i,L)
         endif
        fmnap(i) = alpham_naf(k2) * (1.d0 - mnap(i,L)) -
     X              betam_naf(k2) * mnap(i,L)
        fmkdr(i) = alpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X              betam_kdr(k1) * mkdr(i,L)
        fmka(i)  = alpham_ka (k1) * (1.d0 - mka(i,L)) -
     X              betam_ka (k1) * mka(i,L)
        fhka(i)  = alphah_ka (k1) * (1.d0 - hka(i,L)) -
     X              betah_ka (k1) * hka(i,L)
        fmk2(i)  = alpham_k2 (k1) * (1.d0 - mk2(i,L)) -
     X              betam_k2 (k1) * mk2(i,L)
        fhk2(i)  = alphah_k2 (k1) * (1.d0 - hk2(i,L)) -
     X              betah_k2 (k1) * hk2(i,L)
        fmkm(i)  = alpham_km (k1) * (1.d0 - mkm(i,L)) -
     X              betam_km (k1) * mkm(i,L)
        fmkc(i)  = alpham_kc (k1) * (1.d0 - mkc(i,L)) -
     X              betam_kc (k1) * mkc(i,L)
        fmcat(i) = alpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X              betam_cat(k1) * mcat(i,L)
        fhcat(i) = alphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X              betah_cat(k1) * hcat(i,L)
        fmcaL(i) = alpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X              betam_caL(k1) * mcaL(i,L)
        fmar(i)  = alpham_ar (k1) * (1.d0 - mar(i,L)) -
     X              betam_ar (k1) * mar(i,L)

       dfmnaf_dv(i) = dalpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X                  dbetam_naf(k0) * mnaf(i,L)
       dfmnap_dv(i) = dalpham_naf(k2) * (1.d0 - mnap(i,L)) -
     X                  dbetam_naf(k2) * mnap(i,L)
       dfhnaf_dv(i) = dalphah_naf(k1) * (1.d0 - hnaf(i,L)) -
     X                  dbetah_naf(k1) * hnaf(i,L)
       dfmkdr_dv(i) = dalpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X                  dbetam_kdr(k1) * mkdr(i,L)
       dfmka_dv(i)  = dalpham_ka(k1) * (1.d0 - mka(i,L)) -
     X                  dbetam_ka(k1) * mka(i,L)
       dfhka_dv(i)  = dalphah_ka(k1) * (1.d0 - hka(i,L)) -
     X                  dbetah_ka(k1) * hka(i,L)
       dfmk2_dv(i)  = dalpham_k2(k1) * (1.d0 - mk2(i,L)) -
     X                  dbetam_k2(k1) * mk2(i,L)
       dfhk2_dv(i)  = dalphah_k2(k1) * (1.d0 - hk2(i,L)) -
     X                  dbetah_k2(k1) * hk2(i,L)
       dfmkm_dv(i)  = dalpham_km(k1) * (1.d0 - mkm(i,L)) -
     X                  dbetam_km(k1) * mkm(i,L)
       dfmkc_dv(i)  = dalpham_kc(k1) * (1.d0 - mkc(i,L)) -
     X                  dbetam_kc(k1) * mkc(i,L)
       dfmcat_dv(i) = dalpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X                  dbetam_cat(k1) * mcat(i,L)
       dfhcat_dv(i) = dalphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X                  dbetah_cat(k1) * hcat(i,L)
       dfmcaL_dv(i) = dalpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X                  dbetam_caL(k1) * mcaL(i,L)
       dfmar_dv(i)  = dalpham_ar(k1) * (1.d0 - mar(i,L)) -
     X                  dbetam_ar(k1) * mar(i,L)

       dfmnaf_dmnaf(i) =  - alpham_naf(k0) - betam_naf(k0)
       dfmnap_dmnap(i) =  - alpham_naf(k2) - betam_naf(k2)
       dfhnaf_dhnaf(i) =  - alphah_naf(k1) - betah_naf(k1)
       dfmkdr_dmkdr(i) =  - alpham_kdr(k1) - betam_kdr(k1)
       dfmka_dmka(i)  =   - alpham_ka (k1) - betam_ka (k1)
       dfhka_dhka(i)  =   - alphah_ka (k1) - betah_ka (k1)
       dfmk2_dmk2(i)  =   - alpham_k2 (k1) - betam_k2 (k1)
       dfhk2_dhk2(i)  =   - alphah_k2 (k1) - betah_k2 (k1)
       dfmkm_dmkm(i)  =   - alpham_km (k1) - betam_km (k1)
       dfmkc_dmkc(i)  =   - alpham_kc (k1) - betam_kc (k1)
       dfmcat_dmcat(i) =  - alpham_cat(k1) - betam_cat(k1)
       dfhcat_dhcat(i) =  - alphah_cat(k1) - betah_cat(k1)
       dfmcaL_dmcaL(i) =  - alpham_caL(k1) - betam_caL(k1)
       dfmar_dmar(i)  =   - alpham_ar (k1) - betam_ar (k1)

          end do

       dt2 = 0.5d0 * dt * dt

        do i = 1, numcomp
          v(i,L) = v(i,L) + dt * fv(i)
           do j = 1, numcomp
        v(i,L) = v(i,L) + dt2 * dfv_dv(i,j) * fv(j)
           end do
        v(i,L) = v(i,L) + dt2 * ( dfv_dchi(i) * fchi(i)
     X          + dfv_dmnaf(i) * fmnaf(i)
     X          + dfv_dmnap(i) * fmnap(i)
     X          + dfv_dhnaf(i) * fhnaf(i)
     X          + dfv_dmkdr(i) * fmkdr(i)
     X          + dfv_dmka(i)  * fmka(i)
     X          + dfv_dhka(i)  * fhka(i)
     X          + dfv_dmk2(i)  * fmk2(i)
     X          + dfv_dhk2(i)  * fhk2(i)
     X          + dfv_dmkm(i)  * fmkm(i)
     X          + dfv_dmkc(i)  * fmkc(i)
     X          + dfv_dmkahp(i)* fmkahp(i)
     X          + dfv_dmcat(i)  * fmcat(i)
     X          + dfv_dhcat(i) * fhcat(i)
     X          + dfv_dmcaL(i) * fmcaL(i)
     X          + dfv_dmar(i)  * fmar(i) )

        chi(i,L) = chi(i,L) + dt * fchi(i) + dt2 *
     X   (dfchi_dchi(i) * fchi(i) + dfchi_dv(i) * fv(i))
        mnaf(i,L) = mnaf(i,L) + dt * fmnaf(i) + dt2 *
     X   (dfmnaf_dmnaf(i) * fmnaf(i) + dfmnaf_dv(i)*fv(i))
        mnap(i,L) = mnap(i,L) + dt * fmnap(i) + dt2 *
     X   (dfmnap_dmnap(i) * fmnap(i) + dfmnap_dv(i)*fv(i))
        hnaf(i,L) = hnaf(i,L) + dt * fhnaf(i) + dt2 *
     X   (dfhnaf_dhnaf(i) * fhnaf(i) + dfhnaf_dv(i)*fv(i))
        mkdr(i,L) = mkdr(i,L) + dt * fmkdr(i) + dt2 *
     X   (dfmkdr_dmkdr(i) * fmkdr(i) + dfmkdr_dv(i)*fv(i))
        mka(i,L) =  mka(i,L) + dt * fmka(i) + dt2 *
     X   (dfmka_dmka(i) * fmka(i) + dfmka_dv(i) * fv(i))
        hka(i,L) =  hka(i,L) + dt * fhka(i) + dt2 *
     X   (dfhka_dhka(i) * fhka(i) + dfhka_dv(i) * fv(i))
        mk2(i,L) =  mk2(i,L) + dt * fmk2(i) + dt2 *
     X   (dfmk2_dmk2(i) * fmk2(i) + dfmk2_dv(i) * fv(i))
        hk2(i,L) =  hk2(i,L) + dt * fhk2(i) + dt2 *
     X   (dfhk2_dhk2(i) * fhk2(i) + dfhk2_dv(i) * fv(i))
        mkm(i,L) =  mkm(i,L) + dt * fmkm(i) + dt2 *
     X   (dfmkm_dmkm(i) * fmkm(i) + dfmkm_dv(i) * fv(i))
        mkc(i,L) =  mkc(i,L) + dt * fmkc(i) + dt2 *
     X   (dfmkc_dmkc(i) * fmkc(i) + dfmkc_dv(i) * fv(i))
        mkahp(i,L) = mkahp(i,L) + dt * fmkahp(i) + dt2 *
     X (dfmkahp_dmkahp(i)*fmkahp(i) + dfmkahp_dchi(i)*fchi(i))
        mcat(i,L) =  mcat(i,L) + dt * fmcat(i) + dt2 *
     X   (dfmcat_dmcat(i) * fmcat(i) + dfmcat_dv(i) * fv(i))
        hcat(i,L) =  hcat(i,L) + dt * fhcat(i) + dt2 *
     X   (dfhcat_dhcat(i) * fhcat(i) + dfhcat_dv(i) * fv(i))
        mcaL(i,L) =  mcaL(i,L) + dt * fmcaL(i) + dt2 *
     X   (dfmcaL_dmcaL(i) * fmcaL(i) + dfmcaL_dv(i) * fv(i))
        mar(i,L) =   mar(i,L) + dt * fmar(i) + dt2 *
     X   (dfmar_dmar(i) * fmar(i) + dfmar_dv(i) * fv(i))
c            endif
         end do

! Add membrane currents into membcurr for appropriate compartments
          do i = 1, 9
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 14, 21
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 26, 33
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 39, 68
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do

            end do
c Finish loop L = 1 to numcell

         field_1mm = 0.d0
         field_2mm = 0.d0

         do i = 1, 12
          field_1mm = field_1mm + membcurr(i) / dabs(1000.d0 - depth(i))
          field_2mm = field_2mm + membcurr(i) / dabs(2000.d0 - depth(i))
         end do

2001          CONTINUE

6000    END



C  SETS UP TABLES FOR RATE FUNCTIONS
       SUBROUTINE SCORT_SETUP_suppyrRS
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)
      INTEGER I,J,K
      real*8 minf, hinf, taum, tauh, V, Z, shift_hnaf,
     X  shift_mkdr,
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
C FOR VOLTAGE, RANGE IS -120 TO +40 MV (absol.), 0.25 MV RESOLUTION


       DO 1, I = 0, 640
          V = dble(I)
          V = (V / 4.d0) - 120.d0

c gNa
           minf = 1.d0/(1.d0 + dexp((-V-38.d0)/10.d0))
           if (v.le.-30.d0) then
            taum = .025d0 + .14d0*dexp((v+30.d0)/10.d0)
           else
            taum = .02d0 + .145d0*dexp((-v-30.d0)/10.d0)
           endif
c from principal c. data, Martina & Jonas 1997, tau x 0.5
c Note that minf about the same for interneuron & princ. cell.
           alpham_naf(i) = minf / taum
           betam_naf(i) = 1.d0/taum - alpham_naf(i)

            shift_hnaf =  0.d0
        hinf = 1.d0/(1.d0 +
     x     dexp((v + shift_hnaf + 62.9d0)/10.7d0))
        tauh = 0.15d0 + 1.15d0/(1.d0+dexp((v+37.d0)/15.d0))
c from princ. cell data, Martina & Jonas 1997, tau x 0.5
            alphah_naf(i) = hinf / tauh
            betah_naf(i) = 1.d0/tauh - alphah_naf(i)

          shift_mkdr = 0.d0
c delayed rectifier, non-inactivating
       minf = 1.d0/(1.d0+dexp((-v-shift_mkdr-29.5d0)/10.0d0))
            if (v.le.-10.d0) then
             taum = .25d0 + 4.35d0*dexp((v+10.d0)/10.d0)
            else
             taum = .25d0 + 4.35d0*dexp((-v-10.d0)/10.d0)
            endif
              alpham_kdr(i) = minf / taum
              betam_kdr(i) = 1.d0 /taum - alpham_kdr(i)
c from Martina, Schultz et al., 1998. See espec. Table 1.

c A current: Huguenard & McCormick 1992, J Neurophysiol (TCR)
            minf = 1.d0/(1.d0 + dexp((-v-60.d0)/8.5d0))
            hinf = 1.d0/(1.d0 + dexp((v+78.d0)/6.d0))
        taum = .185d0 + .5d0/(dexp((v+35.8d0)/19.7d0) +
     x                            dexp((-v-79.7d0)/12.7d0))
        if (v.le.-63.d0) then
         tauh = .5d0/(dexp((v+46.d0)/5.d0) +
     x                  dexp((-v-238.d0)/37.5d0))
        else
         tauh = 9.5d0
        endif
           alpham_ka(i) = minf/taum
           betam_ka(i) = 1.d0 / taum - alpham_ka(i)
           alphah_ka(i) = hinf / tauh
           betah_ka(i) = 1.d0 / tauh - alphah_ka(i)

c h-current (anomalous rectifier), Huguenard & McCormick, 1992
           minf = 1.d0/(1.d0 + dexp((v+75.d0)/5.5d0))
           taum = 1.d0/(dexp(-14.6d0 -0.086d0*v) +
     x                   dexp(-1.87 + 0.07d0*v))
           alpham_ar(i) = minf / taum
           betam_ar(i) = 1.d0 / taum - alpham_ar(i)

c K2 K-current, McCormick & Huguenard
             minf = 1.d0/(1.d0 + dexp((-v-10.d0)/17.d0))
             hinf = 1.d0/(1.d0 + dexp((v+58.d0)/10.6d0))
            taum = 4.95d0 + 0.5d0/(dexp((v-81.d0)/25.6d0) +
     x                  dexp((-v-132.d0)/18.d0))
            tauh = 60.d0 + 0.5d0/(dexp((v-1.33d0)/200.d0) +
     x                  dexp((-v-130.d0)/7.1d0))
             alpham_k2(i) = minf / taum
             betam_k2(i) = 1.d0/taum - alpham_k2(i)
             alphah_k2(i) = hinf / tauh
             betah_k2(i) = 1.d0 / tauh - alphah_k2(i)

c voltage part of C-current, using 1994 kinetics, shift 60 mV
              if (v.le.-10.d0) then
       alpham_kc(i) = (2.d0/37.95d0)*dexp((v+50.d0)/11.d0 -
     x                                     (v+53.5)/27.d0)
       betam_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)-alpham_kc(i)
               else
       alpham_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)
       betam_kc(i) = 0.d0
               endif

c high-threshold gCa, from 1994, with 60 mV shift & no inactivn.
            alpham_cal(i) = 1.6d0/(1.d0+dexp(-.072d0*(v-5.d0)))
            betam_cal(i) = 0.1d0 * ((v+8.9d0)/5.d0) /
     x          (dexp((v+8.9d0)/5.d0) - 1.d0)

c M-current, from plast.f, with 60 mV shift
        alpham_km(i) = .02d0/(1.d0+dexp((-v-20.d0)/5.d0))
        betam_km(i) = .01d0 * dexp((-v-43.d0)/18.d0)

c T-current, from Destexhe, Neubig et al., 1998
         minf = 1.d0/(1.d0 + dexp((-v-56.d0)/6.2d0))
         hinf = 1.d0/(1.d0 + dexp((v+80.d0)/4.d0))
         taum = 0.204d0 + .333d0/(dexp((v+15.8d0)/18.2d0) +
     x                  dexp((-v-131.d0)/16.7d0))
          if (v.le.-81.d0) then
         tauh = 0.333 * dexp((v+466.d0)/66.6d0)
          else
         tauh = 9.32d0 + 0.333d0*dexp((-v-21.d0)/10.5d0)
          endif
              alpham_cat(i) = minf / taum
              betam_cat(i) = 1.d0/taum - alpham_cat(i)
              alphah_cat(i) = hinf / tauh
              betah_cat(i) = 1.d0 / tauh - alphah_cat(i)

1        CONTINUE

         do  i = 0, 639

      dalpham_naf(i) = (alpham_naf(i+1)-alpham_naf(i))/.25d0
      dbetam_naf(i) = (betam_naf(i+1)-betam_naf(i))/.25d0
      dalphah_naf(i) = (alphah_naf(i+1)-alphah_naf(i))/.25d0
      dbetah_naf(i) = (betah_naf(i+1)-betah_naf(i))/.25d0
      dalpham_kdr(i) = (alpham_kdr(i+1)-alpham_kdr(i))/.25d0
      dbetam_kdr(i) = (betam_kdr(i+1)-betam_kdr(i))/.25d0
      dalpham_ka(i) = (alpham_ka(i+1)-alpham_ka(i))/.25d0
      dbetam_ka(i) = (betam_ka(i+1)-betam_ka(i))/.25d0
      dalphah_ka(i) = (alphah_ka(i+1)-alphah_ka(i))/.25d0
      dbetah_ka(i) = (betah_ka(i+1)-betah_ka(i))/.25d0
      dalpham_k2(i) = (alpham_k2(i+1)-alpham_k2(i))/.25d0
      dbetam_k2(i) = (betam_k2(i+1)-betam_k2(i))/.25d0
      dalphah_k2(i) = (alphah_k2(i+1)-alphah_k2(i))/.25d0
      dbetah_k2(i) = (betah_k2(i+1)-betah_k2(i))/.25d0
      dalpham_km(i) = (alpham_km(i+1)-alpham_km(i))/.25d0
      dbetam_km(i) = (betam_km(i+1)-betam_km(i))/.25d0
      dalpham_kc(i) = (alpham_kc(i+1)-alpham_kc(i))/.25d0
      dbetam_kc(i) = (betam_kc(i+1)-betam_kc(i))/.25d0
      dalpham_cat(i) = (alpham_cat(i+1)-alpham_cat(i))/.25d0
      dbetam_cat(i) = (betam_cat(i+1)-betam_cat(i))/.25d0
      dalphah_cat(i) = (alphah_cat(i+1)-alphah_cat(i))/.25d0
      dbetah_cat(i) = (betah_cat(i+1)-betah_cat(i))/.25d0
      dalpham_caL(i) = (alpham_cal(i+1)-alpham_cal(i))/.25d0
      dbetam_caL(i) = (betam_cal(i+1)-betam_cal(i))/.25d0
      dalpham_ar(i) = (alpham_ar(i+1)-alpham_ar(i))/.25d0
      dbetam_ar(i) = (betam_ar(i+1)-betam_ar(i))/.25d0
       end do
2      CONTINUE

         do i = 640, 640
      dalpham_naf(i) =  dalpham_naf(i-1)
      dbetam_naf(i) =  dbetam_naf(i-1)
      dalphah_naf(i) = dalphah_naf(i-1)
      dbetah_naf(i) = dbetah_naf(i-1)
      dalpham_kdr(i) =  dalpham_kdr(i-1)
      dbetam_kdr(i) =  dbetam_kdr(i-1)
      dalpham_ka(i) =  dalpham_ka(i-1)
      dbetam_ka(i) =  dbetam_ka(i-1)
      dalphah_ka(i) =  dalphah_ka(i-1)
      dbetah_ka(i) =  dbetah_ka(i-1)
      dalpham_k2(i) =  dalpham_k2(i-1)
      dbetam_k2(i) =  dbetam_k2(i-1)
      dalphah_k2(i) =  dalphah_k2(i-1)
      dbetah_k2(i) =  dbetah_k2(i-1)
      dalpham_km(i) =  dalpham_km(i-1)
      dbetam_km(i) =  dbetam_km(i-1)
      dalpham_kc(i) =  dalpham_kc(i-1)
      dbetam_kc(i) =  dbetam_kc(i-1)
      dalpham_cat(i) =  dalpham_cat(i-1)
      dbetam_cat(i) =  dbetam_cat(i-1)
      dalphah_cat(i) =  dalphah_cat(i-1)
      dbetah_cat(i) =  dbetah_cat(i-1)
      dalpham_caL(i) =  dalpham_caL(i-1)
      dbetam_caL(i) =  dbetam_caL(i-1)
      dalpham_ar(i) =  dalpham_ar(i-1)
      dbetam_ar(i) =  dbetam_ar(i-1)
       end do   

4000   END

        SUBROUTINE SCORTMAJ_suppyrRS
C BRANCHED ACTIVE DENDRITES
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM,depth,level)
c Conductances: leak gL, coupling g, delayed rectifier gKDR, A gKA,
c C gKC, AHP gKAHP, K2 gK2, M gKM, low thresh Ca gCAT, high thresh
c gCAL, fast Na gNAF, persistent Na gNAP, h or anom. rectif. gAR.
c Note VAR = equil. potential for anomalous rectifier.
c Soma = comp. 1; 10 dendrites each with 13 compartments, 6-comp. axon
c Drop "glc"-like terms, just using "gl"-like
c CAFOR corresponds to "phi" in Traub et al., 1994
c Consistent set of units: nF, mV, ms, nA, microS

       INTEGER, PARAMETER:: numcomp = 74
! numcomp here must be compatible with numcomp_suppyrRS in calling prog.
        REAL*8 C(numcomp),GL(numcomp), GAM(0:numcomp, 0:numcomp)
        REAL*8 GNAF(numcomp),GCAT(numcomp), GKAHP(numcomp)
        REAL*8 GKDR(numcomp),GKA(numcomp),GKC(numcomp)
        REAL*8 GK2(numcomp),GNAP(numcomp),GAR(numcomp)
        REAL*8 GKM(numcomp), gcal(numcomp), CDENS
        REAL*8 JACOB(numcomp,numcomp),RI_SD,RI_AXON,RM_SD,RM_AXON
        INTEGER LEVEL(numcomp)
        REAL*8 GNAF_DENS(0:12), GCAT_DENS(0:12), GKDR_DENS(0:12)
        REAL*8 GKA_DENS(0:12), GKC_DENS(0:12), GKAHP_DENS(0:12)
        REAL*8 GCAL_DENS(0:12), GK2_DENS(0:12), GKM_DENS(0:12)
        REAL*8 GNAP_DENS(0:12), GAR_DENS(0:12)
        REAL*8 RES, RINPUT, Z, ELEN(numcomp)
        REAL*8 RSOMA, PI, BETCHI(numcomp), CAFOR(numcomp)
        REAL*8 RAD(numcomp), LEN(numcomp), GAM1, GAM2
        REAL*8 RIN, D(numcomp), AREA(numcomp), RI
        INTEGER NEIGH(numcomp,10), NNUM(numcomp), i, j, k, it
C FOR ESTABLISHING TOPOLOGY OF COMPARTMENTS
        real*8 depth(12) ! depth in microns of levels 1-12, assuming soma in middle
! of layer 2/3 at depth 850 mu

        depth(1) = 850.d0
        depth(2) = 885.d0
        depth(3) = 920.d0
        depth(4) = 955.d0
        depth(5) = 825.d0
        depth(6) = 775.d0
        depth(7) = 725.d0
        depth(8) = 690.d0
        depth(9) = 655.d0
        depth(10) = 620.d0
        depth(11) = 585.d0
        depth(12) = 550.d0

        RI_SD = 250.d0
        RM_SD = 50000.d0
        RI_AXON = 100.d0
        RM_AXON = 1000.d0
        CDENS = 0.9d0

        PI = 3.14159d0

       do i = 0, 12
        gnaf_dens(i) = 10.d0
       end do
c       gnaf_dens(0) = 400.d0
        gnaf_dens(0) = 120.d0
c       gnaf_dens(0) = 100.d0
c       gnaf_dens(1) = 120.d0
        gnaf_dens(1) = 110.d0
        gnaf_dens(2) =  75.d0
        gnaf_dens(5) = 100.d0
        gnaf_dens(6) =  75.d0

       do i = 0, 12
        gkdr_dens(i) = 0.d0
       end do
c       gkdr_dens(0) = 400.d0
c       gkdr_dens(0) = 100.d0
c       gkdr_dens(0) = 150.d0
!       gkdr_dens(0) = 200.d0
        gkdr_dens(0) = 300.d0
c       gkdr_dens(1) = 100.d0
        gkdr_dens(1) = 150.d0
        gkdr_dens(2) =  75.d0
        gkdr_dens(5) = 100.d0
        gkdr_dens(6) =  75.d0

        gnap_dens(0) = 0.d0
        do i = 1, 12
          gnap_dens(i) = 0.0040d0 * gnaf_dens(i)
c         gnap_dens(i) = 0.002d0 * gnaf_dens(i)
c         gnap_dens(i) = 0.0030d0 * gnaf_dens(i)
        end do

        gcat_dens(0) = 0.d0
        do i = 1, 12
c         gcat_dens(i) = 0.5d0
!         gcat_dens(i) = 0.1d0
          gcat_dens(i) = 0.0d0
        end do

        gcaL_dens(0) = 0.d0
        do i = 1, 6
!         gcaL_dens(i) = 0.5d0
          gcaL_dens(i) = 0.80d0
        end do
        gCaL_dens(1) = 0.2d0
        gCaL_dens(2) = 0.2d0
        gCaL_dens(3) = 0.2d0 
        do i = 7, 12
!         gcaL_dens(i) = 0.5d0
          gcaL_dens(i) = 2.0d0 
        end do

       do i = 0, 12
c       gka_dens(i) = 3.d0
        gka_dens(i) = 4.d0
       end do
        gka_dens(1) = 30.d0
        gka_dens(5) = 30.d0

      do i = 0, 12
c        gkc_dens(i)  = 12.00d0
         gkc_dens(i)  =  0.00d0
c        gkc_dens(i)  =  2.00d0
c        gkc_dens(i)  =  7.00d0
      end do
         gkc_dens(0) =  0.00d0
         gkc_dens(1) = 7.5d0
c        gkc_dens(1) =12.5d0
         gkc_dens(2) = 7.5d0
         gkc_dens(5) = 7.5d0
         gkc_dens(6) = 7.5d0

c       gkm_dens(0) = 2.d0 ! 9 Nov. 2005, see scort-pan.f of today
        gkm_dens(0) = 4.d0 ! 9 Nov. 2005, see scort-pan.f of today
! Above suppresses doublets, but still allows FRB with appropriate
! gNaP, gKC, and rel_axonshift (e.g. 6 mV)
        do i = 1, 12
         gkm_dens(i) = 2.5d0 * 1.50d0
        end do

        do i = 0, 12
c       gk2_dens(i) = 1.d0
!       gk2_dens(i) = 0.1d0
        gk2_dens(i) = 0.0d0
        end do
        gk2_dens(0) = 0.d0

        gkahp_dens(0) = 0.d0
        do i = 1, 12
c        gkahp_dens(i) = 0.200d0
!        gkahp_dens(i) = 0.100d0
c        gkahp_dens(i) = 0.030d0
         gkahp_dens(i) = 0.000d0
        end do

        gar_dens(0) = 0.d0
        do i = 1, 12
c        gar_dens(i) = 0.25d0
         gar_dens(i) = 0.00d0
        end do

c       WRITE   (6,9988)
9988    FORMAT(2X,'I',4X,'NADENS',' CADENS(T)',' KDRDEN',' KAHPDE',
     X     ' KCDENS',' KADENS')
        DO 9989, I = 0, 12
c         WRITE (6,9990) I, gnaf_dens(i), gcat_dens(i), gkdr_dens(i),
c    X  gkahp_dens(i), gkc_dens(i), gka_dens(i)
9990    FORMAT(2X,I2,2X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2)
9989    CONTINUE


        level(1) = 1
        do i = 2, 13
         level(i) = 2
        end do
        do i = 14, 25
           level(i) = 3
        end do
        do i = 26, 37
           level(i) = 4
        end do
        level(38) = 5
        level(39) = 6
        level(40) = 7
        level(41) = 8
        level(42) = 8
        level(43) = 9
        level(44) = 9
        do i = 45, 52
           level(i) = 10
        end do
        do i = 53, 60
           level(i) = 11
        end do
        do i = 61, 68
           level(i) = 12
        end do

        do i =  69, 74
         level(i) = 0
        end do

c connectivity of axon
        nnum( 69) = 2
        nnum( 70) = 3
        nnum( 71) = 3
        nnum( 73) = 3
        nnum( 72) = 1
        nnum( 74) = 1
         neigh(69,1) =  1
         neigh(69,2) = 70
         neigh(70,1) = 69
         neigh(70,2) = 71
         neigh(70,3) = 73
         neigh(71,1) = 70
         neigh(71,2) = 72
         neigh(71,3) = 73
         neigh(73,1) = 70
         neigh(73,2) = 71
         neigh(73,3) = 74
         neigh(72,1) = 71
         neigh(74,1) = 73

c connectivity of SD part
          nnum(1) = 10
          neigh(1,1) = 69
          neigh(1,2) =  2
          neigh(1,3) =  3
          neigh(1,4) =  4
          neigh(1,5) =  5
          neigh(1,6) =  6
          neigh(1,7) =  7
          neigh(1,8) =  8
          neigh(1,9) =  9
          neigh(1,10) = 38

          do i = 2, 9
           nnum(i) = 2
           neigh(i,1) = 1
           neigh(i,2) = i + 12
          end do

          do i = 14, 21
            nnum(i) = 2
            neigh(i,1) = i - 12
            neigh(i,2) = i + 12
          end do

          do i = 26, 33
            nnum(i) = 1
            neigh(i,1) = i - 12
          end do

          do i = 10, 13
            nnum(i) = 2
            neigh(i,1) = 38
            neigh(i,2) = i + 12
          end do

          do i = 22, 25
            nnum(i) = 2
            neigh(i,1) = i - 12
            neigh(i,2) = i + 12
          end do

          do i = 34, 37
            nnum(i) = 1
            neigh(i,1) = i - 12
          end do

          nnum(38) = 6
          neigh(38,1) = 1
          neigh(38,2) = 39
          neigh(38,3) = 10
          neigh(38,4) = 11
          neigh(38,5) = 12
          neigh(38,6) = 13

          nnum(39) = 2
          neigh(39,1) = 38
          neigh(39,2) = 40

          nnum(40) = 3
          neigh(40,1) = 39
          neigh(40,2) = 41
          neigh(40,3) = 42

          nnum(41) = 3
          neigh(41,1) = 40
          neigh(41,2) = 42
          neigh(41,3) = 43

          nnum(42) = 3
          neigh(42,1) = 40
          neigh(42,2) = 41
          neigh(42,3) = 44

           nnum(43) = 5
           neigh(43,1) = 41
           neigh(43,2) = 45
           neigh(43,3) = 46
           neigh(43,4) = 47
           neigh(43,5) = 48

           nnum(44) = 5
           neigh(44,1) = 42
           neigh(44,2) = 49
           neigh(44,3) = 50
           neigh(44,4) = 51
           neigh(44,5) = 52

           nnum(45) = 5
           neigh(45,1) = 43
           neigh(45,2) = 53
           neigh(45,3) = 46
           neigh(45,4) = 47
           neigh(45,5) = 48

           nnum(46) = 5
           neigh(46,1) = 43
           neigh(46,2) = 54
           neigh(46,3) = 45
           neigh(46,4) = 47
           neigh(46,5) = 48

           nnum(47) = 5
           neigh(47,1) = 43
           neigh(47,2) = 55
           neigh(47,3) = 45
           neigh(47,4) = 46
           neigh(47,5) = 48

           nnum(48) = 5
           neigh(48,1) = 43
           neigh(48,2) = 56
           neigh(48,3) = 45
           neigh(48,4) = 46
           neigh(48,5) = 47

           nnum(49) = 5
           neigh(49,1) = 44
           neigh(49,2) = 57
           neigh(49,3) = 50
           neigh(49,4) = 51
           neigh(49,5) = 52

           nnum(50) = 5
           neigh(50,1) = 44
           neigh(50,2) = 58
           neigh(50,3) = 49
           neigh(50,4) = 51
           neigh(50,5) = 52

           nnum(51) = 5
           neigh(51,1) = 44
           neigh(51,2) = 59
           neigh(51,3) = 49
           neigh(51,4) = 50
           neigh(51,5) = 52

           nnum(52) = 5
           neigh(52,1) = 44
           neigh(52,2) = 60
           neigh(52,3) = 49
           neigh(52,4) = 51
           neigh(52,5) = 50

          do i = 53, 60
           nnum(i) = 2
           neigh(i,1) = i - 8
           neigh(i,2) = i + 8
          end do

          do i = 61, 68
           nnum(i) = 1
           neigh(i,1) = i - 8
          end do

c        DO 332, I = 1, 74
         DO I = 1, 74
c          WRITE(6,3330) I, NEIGH(I,1),NEIGH(I,2),NEIGH(I,3),NEIGH(I,4),
c    X NEIGH(I,5),NEIGH(I,6),NEIGH(I,7),NEIGH(I,8),NEIGH(I,9),
c    X NEIGH(I,10)
3330     FORMAT(2X,11I5)
         END DO
332      CONTINUE
c         DO 858, I = 1, 74
          DO I = 1, 74
c          DO 858, J = 1, NNUM(I)
           DO J = 1, NNUM(I)
            K = NEIGH(I,J)
            IT = 0
c           DO 859, L = 1, NNUM(K)
            DO  L = 1, NNUM(K)
             IF (NEIGH(K,L).EQ.I) IT = 1
            END DO
859         CONTINUE
             IF (IT.EQ.0) THEN
c             WRITE(6,8591) I, K
8591          FORMAT(' ASYMMETRY IN NEIGH MATRIX ',I4,I4)
              STOP
             ENDIF
          END DO
          END DO
858       CONTINUE

c length and radius of axonal compartments
c Note shortened "initial segment"
          len(69) = 25.d0
          do i = 70, 74
            len(i) = 50.d0
          end do
          rad( 69) = 0.90d0
c         rad( 69) = 0.80d0
          rad( 70) = 0.7d0
          do i = 71, 74
           rad(i) = 0.5d0
          end do

c  length and radius of SD compartments
          len(1) = 15.d0
          rad(1) =  8.d0

          do i = 2, 68
           len(i) = 50.d0
          end do

          do i = 2, 37
            rad(i) = 0.5d0
          end do

          z = 4.0d0
          rad(38) = z
          rad(39) = 0.9d0 * z
          rad(40) = 0.8d0 * z
          rad(41) = 0.5d0 * z
          rad(42) = 0.5d0 * z
          rad(43) = 0.5d0 * z
          rad(44) = 0.5d0 * z
          do i = 45, 68
           rad(i) = 0.2d0 * z
          end do


c       WRITE(6,919)
919     FORMAT('COMPART.',' LEVEL ',' RADIUS ',' LENGTH(MU)')
c       DO 920, I = 1, 74
c920      WRITE(6,921) I, LEVEL(I), RAD(I), LEN(I)
921     FORMAT(I3,5X,I2,3X,F6.2,1X,F6.1,2X,F4.3)

        DO 120, I = 1, 74
          AREA(I) = 2.d0 * PI * RAD(I) * LEN(I)
      if((i.gt.1).and.(i.le.68)) area(i) = 2.d0 * area(i)
C    CORRECTION FOR CONTRIBUTION OF SPINES TO AREA
          K = LEVEL(I)
          C(I) = CDENS * AREA(I) * (1.D-8)

           if (k.ge.1) then
          GL(I) = (1.D-2) * AREA(I) / RM_SD
           else
          GL(I) = (1.D-2) * AREA(I) / RM_AXON
           endif

          GNAF(I) = GNAF_DENS(K) * AREA(I) * (1.D-5)
          GNAP(I) = GNAP_DENS(K) * AREA(I) * (1.D-5)
          GCAT(I) = GCAT_DENS(K) * AREA(I) * (1.D-5)
          GKDR(I) = GKDR_DENS(K) * AREA(I) * (1.D-5)
          GKA(I) = GKA_DENS(K) * AREA(I) * (1.D-5)
          GKC(I) = GKC_DENS(K) * AREA(I) * (1.D-5)
          GKAHP(I) = GKAHP_DENS(K) * AREA(I) * (1.D-5)
          GKM(I) = GKM_DENS(K) * AREA(I) * (1.D-5)
          GCAL(I) = GCAL_DENS(K) * AREA(I) * (1.D-5)
          GK2(I) = GK2_DENS(K) * AREA(I) * (1.D-5)
          GAR(I) = GAR_DENS(K) * AREA(I) * (1.D-5)
c above conductances should be in microS
120           continue

         Z = 0.d0
c        DO 1019, I = 2, 68
         DO I = 2, 68
           Z = Z + AREA(I)
         END DO
1019     CONTINUE
c        WRITE(6,1020) Z
1020     FORMAT(2X,' TOTAL DENDRITIC AREA ',F7.0)

c       DO 140, I = 1, 74
        DO I = 1, 74
c       DO 140, K = 1, NNUM(I)
        DO K = 1, NNUM(I)
         J = NEIGH(I,K)
           if (level(i).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM1 =100.d0 * PI * RAD(I) * RAD(I) / ( RI * LEN(I) )

           if (level(j).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM2 =100.d0 * PI * RAD(J) * RAD(J) / ( RI * LEN(J) )
         GAM(I,J) = 2.d0/( (1.d0/GAM1) + (1.d0/GAM2) )
	 END DO
	 END DO

140     CONTINUE
c gam computed in microS

c       DO 299, I = 1, 74
        DO I = 1, 74
299       BETCHI(I) = .05d0
        END DO
        BETCHI( 1) =  .01d0

c       DO 300, I = 1, 74
        DO I = 1, 74
c300     D(I) = 2.D-4
300     D(I) = 5.D-4
        END DO
c       DO 301, I = 1, 74
        DO I = 1, 74
         IF (LEVEL(I).EQ.1) D(I) = 2.D-3
        END DO
301     CONTINUE
C  NOTE NOTE NOTE  (DIFFERENT FROM SWONG)


c      DO 160, I = 1, 74
       DO I = 1, 74
160     CAFOR(I) = 5200.d0 / (AREA(I) * D(I))
       END DO
C     NOTE CORRECTION

c       do 200, i = 1, 74
        do i = 1, numcomp
200     C(I) = 1000.d0 * C(I)
        end do
C     TO GO FROM MICROF TO NF.

c     DO 909, I = 1, 74
      DO I = 1, numcomp
       JACOB(I,I) = - GL(I)
c     DO 909, J = 1, NNUM(I)
      DO J = 1, NNUM(I)
         K = NEIGH(I,J)
         IF (I.EQ.K) THEN
c            WRITE(6,510) I
510          FORMAT(' UNEXPECTED SYMMETRY IN NEIGH ',I4)
         ENDIF
         JACOB(I,K) = GAM(I,K)
         JACOB(I,I) = JACOB(I,I) - GAM(I,K)
       END DO
       END DO
909   CONTINUE

c 15 Jan. 2001: make correction for c(i)
          do i = 1, numcomp
          do j = 1, numcomp
             jacob(i,j) = jacob(i,j) / c(i)
          end do
          end do

c      DO 500, I = 1, 74
       DO I = 1, 74
c       WRITE (6,501) I,C(I)
501     FORMAT(1X,I3,' C(I) = ',F7.4)
       END DO
500     CONTINUE
        END

! 9 Sept. 2006: start with integrate_nontuftRSXX.f from isoldeepVFOK, and
! add GABA-B
c 31 Mar. 2005, modify with lower axonal fast gNa, shift axonal gNaF
c 11 Oct. 2004, allow for some of the cells (highest nontuftRS_nFRB
c of them) to be FRB, with altered gNaP, gKC & gCaL - see
c layVtup.f.3Feb04

c 3 Nov. 2003, modify layVrsp.f (layer 5 nontufted. pyr cell with thin
c apical dendrite), for use in groucho.f

c 14 May 2003.  Copy program from Rose and modify for mpi.
c 19 June 2001. Taken from scortpd.f.  Parallel.
c 19 June 2001: layer V RS cell with "thin" dendrite, no apical tuft.
c See Kim & Connors, Mason & Larkman.
c 44 SD compartments, 6 axonal, total 50.
c 5 basal and 6 apical oblique dendrites, each with 3 compartments.
c 10-compartment apical dendrite, no branches (apart from obliques)

c 13 April 2001, version of scortp.f, for looking at dendritic activities.
c  7 April 2001, parallel version of scort.f
c 30 March 2001: layer 2/3 pyramidal cell, with geometry (as much as
c possible) from Guy Major thesis; start with tcr.f.
c Total 74 compartments: 6 axon.  8 basal and 3 oblique dendrites, each
c with 3 compartments: apical shaft and branch; 8 3-segment pieces in
c the "apical tuft".

c Revised tcr.f, using modifications developed in short.f
c 22 Feb. 2001: alter persistent gNaP to have lower threshold and
c 1st power activation; in addition, try increasing activation
c threshold of fast gNa, as per Parri & Crunelli 1998.
c 25 Jan. 2001, single TCR cell, modification of nrt.f
c TCR cell has 10 short dendrites, each with 13 compartments.
c Soma is compartment 1; axon is 132-137, with structure as in
c  nRT cell model.  Each dendrite has 2 layers of trifurcations.

c 28 Dec. 2000, begin converting interneuron program to nRT cell.
c Soma will be comp. 1.  4 equivalent dendrites, each with 13 comps.
c (so 53 SD compartments).  Branching axon with 6 compartments - 59
c compartments in all.  Try one integration program for whole structure.
c Currents: leak, fast Na (naf), persistent Na (nap), fast DR (kdr),
c A-current (ka), K2 current, M-current (km), C current (kc), AHP
c (kahp), T-current (cat), high-thresh. Ca (CAL), h-current = anomalous
c rectifier (ar).

         SUBROUTINE integrate_nontuftRSXXB (O, time, numcell, V, curr,
     &    initialize, firstcell, lastcell,
     &    gAMPA, gNMDA, gGABA_A, gGABA_B,
     &    Mg, gapcon, totaxgj, gjtable, dt,
     &  chi,mnaf,mnap,
     &  hnaf,mkdr,mka,
     &  hka,mk2,hk2,
     &  mkm,mkc,mkahp,
     &  mcat,hcat,mcal,
     &  mar,field_1mm,field_2mm,
     &  nontuftRS_nFRB)

         SAVE

         integer, parameter:: numcomp = 50
c numcomp = number of compartments, including 6 in the axon.

       integer numcell, totaxgj, gjtable(totaxgj,4)
       integer initialize, firstcell, lastcell
       INTEGER J1, I, J, K, L, O, k0, K1, k2, nontuftRS_nFRB
       REAL*8  Z, Z1, Z2, Z3, curr(numcomp,numcell)
       REAL*8  mg, time, gapcon, dt, c(numcomp)

c CINV is 1/C, i.e. inverse capacitance

       real*8 v(numcomp,numcell),chi(numcomp,numcell),
     x  mnaf(numcomp,numcell),mnap(numcomp,numcell),
     x hnaf(numcomp,numcell),mkdr(numcomp,numcell),
     x mka(numcomp,numcell),hka(numcomp,numcell),
     x mk2(numcomp,numcell), cinv(numcomp),
     x hk2(numcomp,numcell),mkm(numcomp,numcell),
     x mkc(numcomp,numcell),mkahp(numcomp,numcell),
     x mcat(numcomp,numcell),hcat(numcomp,numcell),
     x mcal(numcomp,numcell),
     x mar(numcomp,numcell),
     x jacob(numcomp,numcomp),betchi(numcomp),
     x gam(0: numcomp,0: numcomp),gL(numcomp),gnaf(numcomp),
     x gnap(numcomp),gkdr(numcomp),gka(numcomp),
     x gk2(numcomp),gkm(numcomp),gkc(numcomp),gkahp(numcomp),
     x gcat(numcomp),gcaL(numcomp),gar(numcomp),
     x gampa(numcomp,numcell),gnmda(numcomp,numcell),
     x ggaba_a(numcomp,numcell),cafor(numcomp),
     x ggaba_b(numcomp,numcell)

! define separate conductance parameters for RS vs FRB behaviors
        real*8 gnap_FRB(numcomp), gkc_FRB(numcomp), gcal_FRB(numcomp)
        real*8 gnap_RS(numcomp), gkc_RS(numcomp), gcal_RS(numcomp)

       real*8
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
       real*8 vL(numcomp),vk(numcomp),vna,var,vca,vgaba_a
       real*8 outrcd(20), A, BB1, BB2

        INTEGER NEIGH(numcomp, 7), NNUM(numcomp)

c the f's are the functions giving 1st derivatives for evolution of
c the differential equations for the voltages (v), calcium (chi), and
c other state variables.
       real*8 fv(numcomp), fchi(numcomp),
     x fmnaf(numcomp),fhnaf(numcomp),fmkdr(numcomp),
     x fmka(numcomp),fhka(numcomp),fmar(numcomp),
     x fmk2(numcomp),fhk2(numcomp),fmnap(numcomp),
     x fmkm(numcomp),fmkc(numcomp),fmkahp(numcomp),
     x fmcat(numcomp),fhcat(numcomp),fmcal(numcomp)

c below are for calculating the partial derivatives
       real*8 dfv_dv(numcomp,numcomp), dfv_dchi(numcomp),
     x    dfv_dmnaf(numcomp),
     x    dfv_dmnap(numcomp),
     x  dfv_dhnaf(numcomp),dfv_dmkdr(numcomp),
     x  dfv_dmka(numcomp),dfv_dhka(numcomp),
     x  dfv_dmk2(numcomp),dfv_dhk2(numcomp),
     x  dfv_dmkm(numcomp),dfv_dmkc(numcomp),
     x  dfv_dmkahp(numcomp),dfv_dmcat(numcomp),
     x  dfv_dhcat(numcomp),dfv_dmcal(numcomp),
     x  dfv_dmar(numcomp)

        real*8 dfchi_dv(numcomp), dfchi_dchi(numcomp),
     x dfmnaf_dmnaf(numcomp), dfmnaf_dv(numcomp),
     x dfhnaf_dhnaf(numcomp),
     x dfmnap_dmnap(numcomp), dfmnap_dv(numcomp),
     x dfhnaf_dv(numcomp),dfmkdr_dmkdr(numcomp),
     x dfmkdr_dv(numcomp),
     x dfmka_dmka(numcomp),dfmka_dv(numcomp),
     x dfhka_dhka(numcomp),dfhka_dv(numcomp),
     x dfmk2_dmk2(numcomp),dfmk2_dv(numcomp),
     x dfhk2_dhk2(numcomp),dfhk2_dv(numcomp),
     x dfmkm_dmkm(numcomp),dfmkm_dv(numcomp),
     x dfmkc_dmkc(numcomp),dfmkc_dv(numcomp),
     x dfmcat_dmcat(numcomp),dfmcat_dv(numcomp),
     x dfhcat_dhcat(numcomp),
     x dfhcat_dv(numcomp),dfmcal_dmcal(numcomp),
     x dfmcal_dv(numcomp),
     x dfmar_dmar(numcomp),dfmar_dv(numcomp),
     x dfmkahp_dchi(numcomp),
     x dfmkahp_dmkahp(numcomp), dt2

       REAL*8 OPEN(numcomp),gamma(numcomp),gamma_prime(numcomp)
c gamma is function of chi used in calculating KC conductance
       REAL*8 alpham_ahp(numcomp), alpham_ahp_prime(numcomp)
       REAL*8 gna_tot(numcomp),gk_tot(numcomp)
       REAL*8 gca_tot(numcomp),gar_tot(numcomp)
       REAL*8 gca_high(numcomp)
c this will be gCa conductance corresponding to high-thresh channels
       real*8 depth(14), membcurr(14), field_1mm, field_2mm
       integer level(numcomp)

       double precision:: persistentNa_shift, fastNa_shift_SD
       double precision::                     fastNa_shift_axon

c Do initialization on 1st time step
        if (initialize.eq.0) then
c       if (O.eq.1) then

c Program assumes A, BB1, BB2 defined in calling program
c as follows:
         A = DEXP(-2.847d0)
         BB1 = DEXP(-.693d0)
         BB2 = DEXP(-3.101d0)

       CALL  nontuftRS_SETUP
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)

        CALL nontuftRSMAJ (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM,depth,level)

          do i = 1, numcomp
             cinv(i) = 1.d0 / c(i)
          end do

             vL = -70.d0

             vK = -95.d0

        VNA = 50.d0
        VCA = 125.d0
        VAR = -43.d0
        VAR = -35.d0
c -43 mV from Huguenard & McCormick
c       VGABA_A = -81.d0
        VGABA_A = -75.d0


c ? initialize membrane state variables?
        v = VL(1)

        k1 = idnint (4.d0 * (v(1,1) + 120.d0))

      hnaf = alphah_naf(k1)/(alphah_naf(k1)+betah_naf(k1))
      hka = alphah_ka(k1)/(alphah_ka(k1)+betah_ka(k1))
      hk2 = alphah_k2(k1)/(alphah_k2(k1)+betah_k2(k1))
      hcat=alphah_cat(k1)/(alphah_cat(k1)+betah_cat(k1))
c     mar=alpham_ar(k1)/(alpham_ar(k1)+betam_ar(k1))
      mar= .25d0
      mnaf = 0.d0
      mkdr = 0.d0
      mka = 0.d0
      mk2 = 0.d0
      mkm = 0.d0
      mkc = 0.d0
      mkahp = 0.d0
      mcat = 0.d0
      mcal = 0.d0
      mnap = 0.d0
      chi = 0.d0

           z1 = 0.10d0
c          z2 = 2.0d0
           z2 = 1.0d0
c This should give an RS cell?

! 1st scale for FRB cells
             do i = 1, numcomp
          gnap_FRB(i) = 1.2d0 * gnap(i)
          gkc_FRB(i)  = 1.d0  * gkc(i)
          gcal_FRB(i) = 0.4d0 * gcal(i)
             end do

! now scale for RS cells
             do i = 1, numcomp
         gnap_RS(i) = z1 * gnap(i)
         gkc_RS (i) = z2 * gkc (i)
         gcal_RS(i) = gcal(i)
             end do

           goto 4000

      endif !  End initialization

             do i = 1, 14
              membcurr(i) = 0.d0
             end do

c         do L = 1, numcell
          do L = firstcell, lastcell

           if (L.le.(numcell - nontuftRS_nFRB)) then
              gnap = gnap_RS
              gkc = gkc_RS
              gcal = gcal_RS
           else
              gnap = gnap_FRB
              gkc  = gkc_FRB
              gcal = gcal_FRB
           endif  

       DO 301, I = 1, numcomp
          FV(I) = -GL(I) * (V(I,L) - VL(i)) * cinv(i)
          DO 302, J = 1, NNUM(I)
             K = NEIGH(I,J)
302     FV(I) = FV(I) + GAM(I,K) * (V(K,L) - V(I,L)) * cinv(i)
301    CONTINUE


        CALL FNMDA (V, OPEN, numcell, numcomp, MG, L, 
     &    A, BB1, BB2)

      DO 421, I = 1, numcomp
       FV(I) = FV(I) + ( CURR(I,L)
     X   - (gampa(I,L) + open(i) * gnmda(I,L))*V(I,L)
     X   - ggaba_a(I,L)*(V(I,L)-Vgaba_a) 
     X   - ggaba_b(I,L)*(V(I,L)-VK(i)  ) ) * cinv(i)
c above assumes equil. potential for AMPA & NMDA = 0 mV
421      continue

! gj code here

       do m = 1, totaxgj
        if (gjtable(m,1).eq.L) then
         L1 = gjtable(m,3)
         igap1 = gjtable(m,2)
         igap2 = gjtable(m,4)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        else if (gjtable(m,3).eq.L) then
         L1 = gjtable(m,1)
         igap1 = gjtable(m,4)
         igap2 = gjtable(m,2)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        endif
       end do ! do m

       do i = 1, numcomp
        gamma(i) = dmin1 (1.d0, .004d0 * chi(i,L))
        if (chi(i,L).le.250.d0) then
          gamma_prime(i) = .004d0
        else
          gamma_prime(i) = 0.d0
        endif
c         endif
       end do

      DO 88, I = 1, numcomp
       gna_tot(i) = gnaf(i) * (mnaf(i,L)**3) * hnaf(i,L) +
     x     gnap(i) * mnap(i,L)
       gk_tot(i) = gkdr(i) * (mkdr(i,L)**4) +
     x             gka(i)  * (mka(i,L)**4) * hka(i,L) +
     x             gk2(i)  * mk2(i,L) * hk2(i,L) +
     x             gkm(i)  * mkm(i,L) +
     x             gkc(i)  * mkc(i,L) * gamma(i) +
     x             gkahp(i)* mkahp(i,L)
       gca_tot(i) = gcat(i) * (mcat(i,L)**2) * hcat(i,L) +
     x              gcaL(i) * (mcaL(i,L)**2)
       gca_high(i) =
     x              gcaL(i) * (mcaL(i,L)**2)
       gar_tot(i) = gar(i) * mar(i,L)


       FV(I) = FV(I) - ( gna_tot(i) * (v(i,L) - vna)
     X  + gk_tot(i) * (v(i,L) - vK(i))
     X  + gca_tot(i) * (v(i,L) - vCa)
     X  + gar_tot(i) * (v(i,L) - var) ) * cinv(i)
c        endif
88           continue

         do i = 1, numcomp
         do j = 1, numcomp
          if (i.ne.j) then
            dfv_dv(i,j) = jacob(i,j)
          else
            dfv_dv(i,j) = jacob(i,i) - cinv(i) *
     X  (gna_tot(i) + gk_tot(i) + gca_tot(i) + gar_tot(i)
     X   + ggaba_a(i,L) + ggaba_b(i,L) + gampa(i,L)
     X   + open(i) * gnmda(I,L) )
          endif
         end do
         end do

           do i = 1, numcomp
        dfv_dchi(i)  = - cinv(i) * gkc(i) * mkc(i,L) *
     x                     gamma_prime(i) * (v(i,L)-vK(i))
        dfv_dmnaf(i) = -3.d0 * cinv(i) * (mnaf(i,L)**2) *
     X    (gnaf(i) * hnaf(i,L)          ) * (v(i,L) - vna)
        dfv_dmnap(i) = - cinv(i) *
     X    (               gnap(i)) * (v(i,L) - vna)
        dfv_dhnaf(i) = - cinv(i) * gnaf(i) * (mnaf(i,L)**3) *
     X                    (v(i,L) - vna)
        dfv_dmkdr(i) = -4.d0 * cinv(i) * gkdr(i) * (mkdr(i,L)**3)
     X                   * (v(i,L) - vK(i))
        dfv_dmka(i)  = -4.d0 * cinv(i) * gka(i) * (mka(i,L)**3) *
     X                   hka(i,L) * (v(i,L) - vK(i))
        dfv_dhka(i)  = - cinv(i) * gka(i) * (mka(i,L)**4) *
     X                    (v(i,L) - vK(i))
        dfv_dmk2(i)  = - cinv(i)*gk2(i)*hk2(i,L)*(v(i,L)-vK(i))
        dfv_dhk2(i)  = - cinv(i)*gk2(i)*mk2(i,L)*(v(i,L)-vK(i))
        dfv_dmkm(i)  = - cinv(i) * gkm(i) * (v(i,L) - vK(i))
      dfv_dmkc(i) = - cinv(i) * gkc(i) * gamma(i)*(v(i,L)-vK(i))
        dfv_dmkahp(i)= - cinv(i) * gkahp(i) * (v(i,L) - vK(i))
        dfv_dmcat(i)  = -2.d0 * cinv(i) * gcat(i) * mcat(i,L) *
     X                    hcat(i,L) * (v(i,L) - vCa)
        dfv_dhcat(i) = - cinv(i) * gcat(i) * (mcat(i,L)**2) *
     X                  (v(i,L) - vCa)
        dfv_dmcal(i) = -2.d0 * cinv(i) * gcal(i) * mcal(i,L) *
     X                      (v(i,L) - vCa)
        dfv_dmar(i) = - cinv(i) * gar(i) * (v(i,L) - var)
            end do

         do i = 1, numcomp
          fchi(i) = - cafor(i) * gca_high(i) * (v(i,L) - vca)
     x       - betchi(i) * chi(i,L)
          dfchi_dv(i) = - cafor(i) * gca_high(i)
          dfchi_dchi(i) = - betchi(i)
         end do

       do i = 1, numcomp
c Note possible increase in rate at which AHP current develops
c       alpham_ahp(i) = dmin1(0.2d-4 * chi(i),0.01d0)
        alpham_ahp(i) = dmin1(1.0d-4 * chi(i,L),0.01d0)
        if (chi(i,L).le.500.d0) then
c         alpham_ahp_prime(i) = 0.2d-4
          alpham_ahp_prime(i) = 1.0d-4
        else
          alpham_ahp_prime(i) = 0.d0
        endif
       end do

       do i = 1, numcomp
        fmkahp(i) = alpham_ahp(i) * (1.d0 - mkahp(i,L))
     x                  -.001d0 * mkahp(i,L)
c    x                  -.010d0 * mkahp(i,L)
        dfmkahp_dmkahp(i) = - alpham_ahp(i) - .001d0
c       dfmkahp_dmkahp(i) = - alpham_ahp(i) - .010d0
        dfmkahp_dchi(i) = alpham_ahp_prime(i) *
     x                     (1.d0 - mkahp(i,L))
       end do

          do i = 1, numcomp

       K1 = IDNINT ( 4.d0 * (V(I,L) + 120.d0) )
       IF (K1.GT.640) K1 = 640
       IF (K1.LT.  0) K1 =   0

c      persistentNa_shift =  0.d0
c      persistentNa_shift =  8.d0
       persistentNa_shift = 10.d0
       K2 = IDNINT ( 4.d0 * (V(I,L)+persistentNa_shift+ 120.d0) )
       IF (K2.GT.640) K2 = 640
       IF (K2.LT.  0) K2 =   0

c            fastNa_shift = -2.0d0
c            fastNa_shift = -2.5d0
             fastNa_shift_SD = -3.5d0
             fastNa_shift_axon = fastNa_shift_SD + 7.d0
       K0 = IDNINT ( 4.d0 * (V(I,L)+  fastNa_shift_SD+ 120.d0) )
       K3 = IDNINT ( 4.d0 * (V(I,L)+  fastNa_shift_axon+ 120.d0) )
       IF (K0.GT.640) K0 = 640
       IF (K0.LT.  0) K0 =   0
       IF (K3.GT.640) K3 = 640
       IF (K3.LT.  0) K3 =   0

         if (i.le.44) then
        fmnaf(i) = alpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k0) * mnaf(i,L)
        fhnaf(i) = alphah_naf(k0) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k0) * hnaf(i,L)
         else
        fmnaf(i) = alpham_naf(k3) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k3) * mnaf(i,L)
        fhnaf(i) = alphah_naf(k3) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k3) * hnaf(i,L)
         endif


        fmnap(i) = alpham_naf(k2) * (1.d0 - mnap(i,L)) -
     X              betam_naf(k2) * mnap(i,L)
        fmkdr(i) = alpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X              betam_kdr(k1) * mkdr(i,L)
        fmka(i)  = alpham_ka (k1) * (1.d0 - mka(i,L)) -
     X              betam_ka (k1) * mka(i,L)
        fhka(i)  = alphah_ka (k1) * (1.d0 - hka(i,L)) -
     X              betah_ka (k1) * hka(i,L)
        fmk2(i)  = alpham_k2 (k1) * (1.d0 - mk2(i,L)) -
     X              betam_k2 (k1) * mk2(i,L)
        fhk2(i)  = alphah_k2 (k1) * (1.d0 - hk2(i,L)) -
     X              betah_k2 (k1) * hk2(i,L)
        fmkm(i)  = alpham_km (k1) * (1.d0 - mkm(i,L)) -
     X              betam_km (k1) * mkm(i,L)
        fmkc(i)  = alpham_kc (k1) * (1.d0 - mkc(i,L)) -
     X              betam_kc (k1) * mkc(i,L)
        fmcat(i) = alpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X              betam_cat(k1) * mcat(i,L)
        fhcat(i) = alphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X              betah_cat(k1) * hcat(i,L)
        fmcaL(i) = alpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X              betam_caL(k1) * mcaL(i,L)
        fmar(i)  = alpham_ar (k1) * (1.d0 - mar(i,L)) -
     X              betam_ar (k1) * mar(i,L)

       dfmnaf_dv(i) = dalpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X                  dbetam_naf(k0) * mnaf(i,L)
       dfmnap_dv(i) = dalpham_naf(k2) * (1.d0 - mnap(i,L)) -
     X                  dbetam_naf(k2) * mnap(i,L)
       dfhnaf_dv(i) = dalphah_naf(k1) * (1.d0 - hnaf(i,L)) -
     X                  dbetah_naf(k1) * hnaf(i,L)
       dfmkdr_dv(i) = dalpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X                  dbetam_kdr(k1) * mkdr(i,L)
       dfmka_dv(i)  = dalpham_ka(k1) * (1.d0 - mka(i,L)) -
     X                  dbetam_ka(k1) * mka(i,L)
       dfhka_dv(i)  = dalphah_ka(k1) * (1.d0 - hka(i,L)) -
     X                  dbetah_ka(k1) * hka(i,L)
       dfmk2_dv(i)  = dalpham_k2(k1) * (1.d0 - mk2(i,L)) -
     X                  dbetam_k2(k1) * mk2(i,L)
       dfhk2_dv(i)  = dalphah_k2(k1) * (1.d0 - hk2(i,L)) -
     X                  dbetah_k2(k1) * hk2(i,L)
       dfmkm_dv(i)  = dalpham_km(k1) * (1.d0 - mkm(i,L)) -
     X                  dbetam_km(k1) * mkm(i,L)
       dfmkc_dv(i)  = dalpham_kc(k1) * (1.d0 - mkc(i,L)) -
     X                  dbetam_kc(k1) * mkc(i,L)
       dfmcat_dv(i) = dalpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X                  dbetam_cat(k1) * mcat(i,L)
       dfhcat_dv(i) = dalphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X                  dbetah_cat(k1) * hcat(i,L)
       dfmcaL_dv(i) = dalpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X                  dbetam_caL(k1) * mcaL(i,L)
       dfmar_dv(i)  = dalpham_ar(k1) * (1.d0 - mar(i,L)) -
     X                  dbetam_ar(k1) * mar(i,L)

       dfmnaf_dmnaf(i) =  - alpham_naf(k0) - betam_naf(k0)
       dfmnap_dmnap(i) =  - alpham_naf(k2) - betam_naf(k2)
       dfhnaf_dhnaf(i) =  - alphah_naf(k1) - betah_naf(k1)
       dfmkdr_dmkdr(i) =  - alpham_kdr(k1) - betam_kdr(k1)
       dfmka_dmka(i)  =   - alpham_ka (k1) - betam_ka (k1)
       dfhka_dhka(i)  =   - alphah_ka (k1) - betah_ka (k1)
       dfmk2_dmk2(i)  =   - alpham_k2 (k1) - betam_k2 (k1)
       dfhk2_dhk2(i)  =   - alphah_k2 (k1) - betah_k2 (k1)
       dfmkm_dmkm(i)  =   - alpham_km (k1) - betam_km (k1)
       dfmkc_dmkc(i)  =   - alpham_kc (k1) - betam_kc (k1)
       dfmcat_dmcat(i) =  - alpham_cat(k1) - betam_cat(k1)
       dfhcat_dhcat(i) =  - alphah_cat(k1) - betah_cat(k1)
       dfmcaL_dmcaL(i) =  - alpham_caL(k1) - betam_caL(k1)
       dfmar_dmar(i)  =   - alpham_ar (k1) - betam_ar (k1)

          end do

       dt2 = 0.5d0 * dt * dt

        do i = 1, numcomp
          v(i,L) = v(i,L) + dt * fv(i)
           do j = 1, numcomp
        v(i,L) = v(i,L) + dt2 * dfv_dv(i,j) * fv(j)
           end do
        v(i,L) = v(i,L) + dt2 * ( dfv_dchi(i) * fchi(i)
     X          + dfv_dmnaf(i) * fmnaf(i)
     X          + dfv_dmnap(i) * fmnap(i)
     X          + dfv_dhnaf(i) * fhnaf(i)
     X          + dfv_dmkdr(i) * fmkdr(i)
     X          + dfv_dmka(i)  * fmka(i)
     X          + dfv_dhka(i)  * fhka(i)
     X          + dfv_dmk2(i)  * fmk2(i)
     X          + dfv_dhk2(i)  * fhk2(i)
     X          + dfv_dmkm(i)  * fmkm(i)
     X          + dfv_dmkc(i)  * fmkc(i)
     X          + dfv_dmkahp(i)* fmkahp(i)
     X          + dfv_dmcat(i)  * fmcat(i)
     X          + dfv_dhcat(i) * fhcat(i)
     X          + dfv_dmcaL(i) * fmcaL(i)
     X          + dfv_dmar(i)  * fmar(i) )

        chi(i,L) = chi(i,L) + dt * fchi(i) + dt2 *
     X   (dfchi_dchi(i) * fchi(i) + dfchi_dv(i) * fv(i))
        mnaf(i,L) = mnaf(i,L) + dt * fmnaf(i) + dt2 *
     X   (dfmnaf_dmnaf(i) * fmnaf(i) + dfmnaf_dv(i)*fv(i))
        mnap(i,L) = mnap(i,L) + dt * fmnap(i) + dt2 *
     X   (dfmnap_dmnap(i) * fmnap(i) + dfmnap_dv(i)*fv(i))
        hnaf(i,L) = hnaf(i,L) + dt * fhnaf(i) + dt2 *
     X   (dfhnaf_dhnaf(i) * fhnaf(i) + dfhnaf_dv(i)*fv(i))
        mkdr(i,L) = mkdr(i,L) + dt * fmkdr(i) + dt2 *
     X   (dfmkdr_dmkdr(i) * fmkdr(i) + dfmkdr_dv(i)*fv(i))
        mka(i,L) =  mka(i,L) + dt * fmka(i) + dt2 *
     X   (dfmka_dmka(i) * fmka(i) + dfmka_dv(i) * fv(i))
        hka(i,L) =  hka(i,L) + dt * fhka(i) + dt2 *
     X   (dfhka_dhka(i) * fhka(i) + dfhka_dv(i) * fv(i))
        mk2(i,L) =  mk2(i,L) + dt * fmk2(i) + dt2 *
     X   (dfmk2_dmk2(i) * fmk2(i) + dfmk2_dv(i) * fv(i))
        hk2(i,L) =  hk2(i,L) + dt * fhk2(i) + dt2 *
     X   (dfhk2_dhk2(i) * fhk2(i) + dfhk2_dv(i) * fv(i))
        mkm(i,L) =  mkm(i,L) + dt * fmkm(i) + dt2 *
     X   (dfmkm_dmkm(i) * fmkm(i) + dfmkm_dv(i) * fv(i))
        mkc(i,L) =  mkc(i,L) + dt * fmkc(i) + dt2 *
     X   (dfmkc_dmkc(i) * fmkc(i) + dfmkc_dv(i) * fv(i))
        mkahp(i,L) = mkahp(i,L) + dt * fmkahp(i) + dt2 *
     X (dfmkahp_dmkahp(i)*fmkahp(i) + dfmkahp_dchi(i)*fchi(i))
        mcat(i,L) =  mcat(i,L) + dt * fmcat(i) + dt2 *
     X   (dfmcat_dmcat(i) * fmcat(i) + dfmcat_dv(i) * fv(i))
        hcat(i,L) =  hcat(i,L) + dt * fhcat(i) + dt2 *
     X   (dfhcat_dhcat(i) * fhcat(i) + dfhcat_dv(i) * fv(i))
        mcaL(i,L) =  mcaL(i,L) + dt * fmcaL(i) + dt2 *
     X   (dfmcaL_dmcaL(i) * fmcaL(i) + dfmcaL_dv(i) * fv(i))
        mar(i,L) =   mar(i,L) + dt * fmar(i) + dt2 *
     X   (dfmar_dmar(i) * fmar(i) + dfmar_dv(i) * fv(i))
c            endif
         end do

! Add membrane currents into membcurr for appropriate compartments
          do i = 1, 6
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 13, 17
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 24, 28
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 35, 44
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do

           end do ! do L

         field_1mm = 0.d0
         field_2mm = 0.d0

         do i = 1, 14
          field_1mm = field_1mm + membcurr(i) / dabs(1000.d0 - depth(i))
          field_2mm = field_2mm + membcurr(i) / dabs(2000.d0 - depth(i))
         end do

4000          END

        SUBROUTINE nontuftRSMAJ
C BRANCHED ACTIVE DENDRITES
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM,depth,level)
c Conductances: leak gL, coupling g, delayed rectifier gKDR, A gKA,
c C gKC, AHP gKAHP, K2 gK2, M gKM, low thresh Ca gCAT, high thresh
c gCAL, fast Na gNAF, persistent Na gNAP, h or anom. rectif. gAR.
c Note VAR = equil. potential for anomalous rectifier.
c Soma = comp. 1; 10 dendrites each with 13 compartments, 6-comp. axon
c Drop "glc"-like terms, just using "gl"-like
c CAFOR corresponds to "phi" in Traub et al., 1994
c Consistent set of units: nF, mV, ms, nA, microS

        integer, parameter:: numcomp = 50

        REAL*8 C(numcomp),GL(numcomp),GAM(0:numcomp,0:numcomp)
        REAL*8 GNAF(numcomp),GCAT(numcomp)
        REAL*8 GKDR(numcomp),GKA(numcomp),GKC(numcomp)
        REAL*8 GKAHP(numcomp),GCAL(numcomp)
        REAL*8 GK2(numcomp),GKM(numcomp),GNAP(numcomp),GAR(numcomp)
        REAL*8 JACOB(numcomp,numcomp),RI_SD,RI_AXON,RM_SD,RM_AXON
        INTEGER LEVEL(numcomp)
        REAL*8 GNAF_DENS(0:14), GCAT_DENS(0:14), GKDR_DENS(0:14)
        REAL*8 GKA_DENS(0:14), GKC_DENS(0:14), GKAHP_DENS(0:14)
        REAL*8 GCAL_DENS(0:14), GK2_DENS(0:14), GKM_DENS(0:14)
        REAL*8 GNAP_DENS(0:14), GAR_DENS(0:14)
        REAL*8 RES, RINPUT, CDENS
        REAL*8 RSOMA, PI, BETCHI(numcomp), CAFOR(numcomp)
        REAL*8 RAD(numcomp),LEN(numcomp),GAM1,GAM2,ELEN(numcomp)
        REAL*8 RIN, D(numcomp), AREA(numcomp), RI, Z
        INTEGER NEIGH(numcomp, 7), NNUM(numcomp)
C FOR ESTABLISHING TOPOLOGY OF COMPARTMENTS
        real*8 depth(14) ! in microns, subscript refers to level

        depth(1) = 2200.d0
        depth(2) = 2245.d0
        depth(3) = 2290.d0
        depth(4) = 2335.d0
        depth(5) = 2175.d0
        depth(6) = 2125.d0
        depth(7) = 2075.d0
        depth(8) = 2025.d0
        depth(9) = 1975.d0
        depth(10) = 1925.d0
        depth(11) = 1875.d0
        depth(12) = 1825.d0
        depth(13) = 1775.d0
        depth(14) = 1725.d0

        RI_SD = 250.d0
        RM_SD = 50000.d0
        RI_AXON = 100.d0
        RM_AXON = 1000.d0
        CDENS = 0.9d0

        PI = 3.14159d0

        gnaf_dens =  5.d0
c       gnaf_dens = 10.d0
c       gnaf_dens(0) = 450.d0
        gnaf_dens(0) = 175.d0
c       gnaf_dens(1) = 200.d0
        gnaf_dens(1) = 175.d0
        gnaf_dens(2) =  75.d0
        gnaf_dens(5) = 150.d0
        gnaf_dens(6) =  75.d0

        gkdr_dens = 0.d0
        gkdr_dens(0) = 450.d0
        gkdr_dens(1) = 170.d0
        gkdr_dens(2) =  75.d0
        gkdr_dens(5) = 120.d0
        gkdr_dens(6) =  75.d0

        do i = 1, 14
          gnap_dens(i) = 0.0040d0 * gnaf_dens(i)
        end do

        do i = 1, 14
          gcat_dens(i) = 0.1d0
        end do

        do i = 1, 9
          gcaL_dens(i) = 0.20d0
        end do
        do i = 10, 14
          gcaL_dens(i) = 2.0d0
        end do

        gka_dens    = 4.d0
        gka_dens(1) = 35.d0
        gka_dens(5) = 35.d0
        do i = 1, 14
          gka_dens(i) = 3.4d0 * gka_dens(i)
        end do

        gkc_dens = 0.00d0
         gkc_dens(1) =  7.50d0
         gkc_dens(2) =  7.50d0
         gkc_dens(5) =  7.50d0
         gkc_dens(6) =  7.50d0

        do i = 1, 14
c        gkm_dens(i) = 1.4d0 * 1.50d0
         gkm_dens(i) = 2.8d0 * 1.50d0
        end do

        gk2_dens    = 0.1d0

        do i = 1, 14
c        gkahp_dens(i) = 0.100d0
         gkahp_dens(i) = 0.200d0
        end do

        do i = 1, 14
         gar_dens(i) = 0.25d0
        end do

c        if (thisno.eq.0) then
c       WRITE   (6,9988)
9988    FORMAT(2X,'I',4X,'NADENS',' CADENS(L)',' KDRDEN',' KAHPDE',
     X     ' KCDENS',' KADENS')
        DO 9989, I = 0, 14
c         WRITE (6,9990) I, gnaf_dens(i), gcaL_dens(i), gkdr_dens(i),
c    X  gkahp_dens(i), gkc_dens(i), gka_dens(i)
9990    FORMAT(2X,I2,2X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2)
9989    CONTINUE
c         endif


        level(1) = 1
        do i = 2, 12
         level(i) = 2
        end do
        do i = 13, 23
           level(i) = 3
        end do
        do i = 24, 34
           level(i) = 4
        end do
        level(35) = 5
        level(36) = 6
        level(37) = 7
        level(38) = 8
        level(39) = 9
        level(40) = 10
        level(41) = 11
        level(42) = 12
        level(43) = 13
        level(44) = 14

        do i =  45, 50
         level(i) = 0
        end do

c connectivity of axon
        nnum( 45) = 2
        nnum( 46) = 3
        nnum( 47) = 3
        nnum( 49) = 3
        nnum( 48) = 1
        nnum( 50) = 1
         neigh(45,1) =  1
         neigh(45,2) = 46
         neigh(46,1) = 45
         neigh(46,2) = 47
         neigh(46,3) = 49
         neigh(47,1) = 46
         neigh(47,2) = 48
         neigh(47,3) = 49
         neigh(49,1) = 46
         neigh(49,2) = 47
         neigh(49,3) = 50
         neigh(48,1) = 47
         neigh(50,1) = 49

c connectivity of SD part
          nnum(1) = 7
          neigh(1,1) = 45
          neigh(1,2) =  2
          neigh(1,3) =  3
          neigh(1,4) =  4
          neigh(1,5) =  5
          neigh(1,6) =  6
          neigh(1,7) = 35

          do i = 2, 6
           nnum(i) = 2
           neigh(i,1) = 1
           neigh(i,2) = i + 11
          end do

          do i = 13, 17
            nnum(i) = 2
            neigh(i,1) = i - 11
            neigh(i,2) = i + 11
          end do

          do i = 24, 28
            nnum(i) = 1
            neigh(i,1) = i - 11
          end do

          do i =  7, 12
            nnum(i) = 2
      if ((i.eq.7).or.(i.eq.12)) neigh(i,1) = 35
      if ((i.eq.8).or.(i.eq.11)) neigh(i,1) = 36
      if ((i.eq.9).or.(i.eq.10)) neigh(i,1) = 37
            neigh(i,2) = i + 11
          end do

          do i = 18, 23
            nnum(i) = 2
            neigh(i,1) = i - 11
            neigh(i,2) = i + 11
          end do

          do i = 29, 34
            nnum(i) = 1
            neigh(i,1) = i - 11
          end do

          nnum(35) = 4
          neigh(35,1) = 1
          neigh(35,2) = 36
          neigh(35,3) =  7
          neigh(35,4) = 12

          nnum(36) = 4
          neigh(36,1) = 35
          neigh(36,2) = 37
          neigh(36,3) =  8
          neigh(36,4) = 11

          nnum(37) = 4
          neigh(37,1) = 36
          neigh(37,2) = 38
          neigh(37,3) =  9
          neigh(37,4) = 10

          nnum(38) = 2
          neigh(38,1) = 37
          neigh(38,2) = 39

          nnum(39) = 2
          neigh(39,1) = 38
          neigh(39,2) = 40

          nnum(40) = 2
          neigh(40,1) = 39
          neigh(40,2) = 41

          nnum(41) = 2
          neigh(41,1) = 40
          neigh(41,2) = 42

          nnum(42) = 2
          neigh(42,1) = 41
          neigh(42,2) = 43

          nnum(43) = 2
          neigh(43,1) = 42
          neigh(43,2) = 44

          nnum(44) = 1
          neigh(44,1) = 43

c           if (thisno.eq.0) then
         DO 332, I = 1, numcomp
c          WRITE(6,3330) I, NEIGH(I,1),NEIGH(I,2),NEIGH(I,3),NEIGH(I,4),
c    X NEIGH(I,5),NEIGH(I,6),NEIGH(I,7)
3330     FORMAT(2X, 8I5)
332      CONTINUE
c            endif

          DO 858, I = 1, numcomp
           DO 858, J = 1, NNUM(I)
            K = NEIGH(I,J)
            IT = 0
            DO 859, L = 1, NNUM(K)
             IF (NEIGH(K,L).EQ.I) IT = 1
859         CONTINUE
             IF (IT.EQ.0) THEN
c             WRITE(6,8591) I, K
8591          FORMAT(' ASYMMETRY IN NEIGH MATRIX ',I4,I4)
              STOP
             ENDIF
858       CONTINUE

c length and radius of axonal compartments
c Note shortened "initial segment"
          len(45) = 25.d0
          do i = 46, 50
            len(i) = 50.d0
          end do
          rad( 45) = 0.90d0
          rad( 46) = 0.7d0
          do i = 47, 50
           rad(i) = 0.5d0
          end do

c  length and radius of SD compartments
          len(1) = 20.d0
          rad(1) =  8.d0

          do i = 2, 34
           len(i) = 60.d0
          end do

          do i = 35, 44
           len(i) = 50.d0
          end do

          do i = 2, 6
            rad(i) = 0.85d0
          end do
          do i = 13, 17
            rad(i) = 0.85d0
          end do
          do i = 24, 28
            rad(i) = 0.85d0
          end do

          do i = 7, 12
            rad(i) = 0.62d0
          end do
          do i = 18, 23
            rad(i) = 0.62d0
          end do
          do i = 29, 34
            rad(i) = 0.62d0
          end do

          rad(35) = 1.5d0
          rad(36) = 1.4d0
          rad(37) = 1.3d0
          rad(38) = 1.2d0
          rad(39) = 1.1d0
          rad(40) = 1.0d0
          rad(41) = 0.9d0
          rad(42) = 0.8d0
          rad(43) = 0.7d0
          rad(44) = 0.6d0


c            if (thisno.eq.0) then
        WRITE(6,919)
919     FORMAT('COMPART.',' LEVEL ',' RADIUS ',' LENGTH(MU)')
c       DO 920, I = 1, numcomp
c920      WRITE(6,921) I, LEVEL(I), RAD(I), LEN(I)
921     FORMAT(I3,5X,I2,3X,F6.2,1X,F6.1,2X,F4.3)
c            endif

        DO 120, I = 1, numcomp
          AREA(I) = 2.d0 * PI * RAD(I) * LEN(I)
      if((i.gt.1).and.(i.le.44)) area(i) = 2.d0 * area(i)
C    CORRECTION FOR CONTRIBUTION OF SPINES TO AREA
          K = LEVEL(I)
          C(I) = CDENS * AREA(I) * (1.D-8)

           if (k.ge.1) then
          GL(I) = (1.D-2) * AREA(I) / RM_SD
           else
          GL(I) = (1.D-2) * AREA(I) / RM_AXON
           endif

          GNAF(I) = GNAF_DENS(K) * AREA(I) * (1.D-5)
          GNAP(I) = GNAP_DENS(K) * AREA(I) * (1.D-5)
          GCAT(I) = GCAT_DENS(K) * AREA(I) * (1.D-5)
          GKDR(I) = GKDR_DENS(K) * AREA(I) * (1.D-5)
          GKA(I) = GKA_DENS(K) * AREA(I) * (1.D-5)
          GKC(I) = GKC_DENS(K) * AREA(I) * (1.D-5)
          GKAHP(I) = GKAHP_DENS(K) * AREA(I) * (1.D-5)
          GCAL(I) = GCAL_DENS(K) * AREA(I) * (1.D-5)
          GK2(I) = GK2_DENS(K) * AREA(I) * (1.D-5)
          GKM(I) = GKM_DENS(K) * AREA(I) * (1.D-5)
          GAR(I) = GAR_DENS(K) * AREA(I) * (1.D-5)
c above conductances should be in microS
120           continue

         Z = 0.d0
         DO 1019, I = 2, 44
           Z = Z + AREA(I)
1019     CONTINUE
c             if (thisno.eq.1) then
c        WRITE(6,1020) Z
c              endif
1020     FORMAT(2X,' TOTAL DENDRITIC AREA ',F7.0)

        DO 140, I = 1, numcomp
        DO 140, K = 1, NNUM(I)
         J = NEIGH(I,K)
           if (level(i).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM1 =100.d0 * PI * RAD(I) * RAD(I) / ( RI * LEN(I) )

           if (level(j).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM2 =100.d0 * PI * RAD(J) * RAD(J) / ( RI * LEN(J) )

         GAM(I,J) = 2.d0/( (1.d0/GAM1) + (1.d0/GAM2) )
140     CONTINUE
c gam computed in microS

        DO 299, I = 1, numcomp
299       BETCHI(I) = .05d0
        BETCHI( 1) =  .01d0

        DO 300, I = 1, numcomp
300     D(I) = 4.D-4
        DO 301, I = 1, numcomp
         IF (LEVEL(I).EQ.1) D(I) = 4.D-3
301     CONTINUE


       DO 160, I = 1, numcomp
160     CAFOR(I) = 5200.d0 / (AREA(I) * D(I))
C     NOTE CORRECTION

        do 200, i = 1, numcomp
200     C(I) = 1000.d0 * C(I)
C     TO GO FROM MICROF TO NF.

      DO 909, I = 1, numcomp
       JACOB(I,I) = - GL(I)
      DO 909, J = 1, NNUM(I)
         K = NEIGH(I,J)
         IF (I.EQ.K) THEN
c            WRITE(6,510) I
510          FORMAT(' UNEXPECTED SYMMETRY IN NEIGH ',I4)
         ENDIF
         JACOB(I,K) = GAM(I,K)
         JACOB(I,I) = JACOB(I,I) - GAM(I,K)
909   CONTINUE

c 15 Jan. 2001: make correction for c(i)
          do i = 1, numcomp
          do j = 1, numcomp
             jacob(i,j) = jacob(i,j) / c(i)
          end do
          end do

c          if (thisno.eq.1) then
       DO 500, I = 1, numcomp
c       WRITE (6,501) I,C(I)
501     FORMAT(1X,I3,' C(I) = ',F7.4)
500     CONTINUE
c               endif

        END

C  SETS UP TABLES FOR RATE FUNCTIONS
       SUBROUTINE nontuftRS_SETUP
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)
      INTEGER I,J,K
      real*8 minf, hinf, taum, tauh, V, Z, shift_hnaf,
     X  shift_mkdr,
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
C FOR VOLTAGE, RANGE IS -120 TO +40 MV (absol.), 0.25 MV RESOLUTION


       DO 1, I = 0, 640
          V = dble(I)
          V = (V / 4.d0) - 120.d0

c gNa
           minf = 1.d0/(1.d0 + dexp((-V-38.d0)/10.d0))
           if (v.le.-30.d0) then
            taum = .025d0 + .14d0*dexp((v+30.d0)/10.d0)
           else
            taum = .02d0 + .145d0*dexp((-v-30.d0)/10.d0)
           endif
c from principal c. data, Martina & Jonas 1997, tau x 0.5
c Note that minf about the same for interneuron & princ. cell.
           alpham_naf(i) = minf / taum
           betam_naf(i) = 1.d0/taum - alpham_naf(i)

            shift_hnaf =  0.d0
        hinf = 1.d0/(1.d0 +
     x     dexp((v + shift_hnaf + 62.9d0)/10.7d0))
        tauh = 0.15d0 + 1.15d0/(1.d0+dexp((v+37.d0)/15.d0))
c from princ. cell data, Martina & Jonas 1997, tau x 0.5
            alphah_naf(i) = hinf / tauh
            betah_naf(i) = 1.d0/tauh - alphah_naf(i)

          shift_mkdr = 0.d0
c delayed rectifier, non-inactivating
       minf = 1.d0/(1.d0+dexp((-v-shift_mkdr-29.5d0)/10.0d0))
            if (v.le.-10.d0) then
             taum = .25d0 + 4.35d0*dexp((v+10.d0)/10.d0)
            else
             taum = .25d0 + 4.35d0*dexp((-v-10.d0)/10.d0)
            endif
              alpham_kdr(i) = minf / taum
              betam_kdr(i) = 1.d0 /taum - alpham_kdr(i)
c from Martina, Schultz et al., 1998. See espec. Table 1.

c A current: Huguenard & McCormick 1992, J Neurophysiol (TCR)
            minf = 1.d0/(1.d0 + dexp((-v-60.d0)/8.5d0))
            hinf = 1.d0/(1.d0 + dexp((v+78.d0)/6.d0))
        taum = .185d0 + .5d0/(dexp((v+35.8d0)/19.7d0) +
     x                            dexp((-v-79.7d0)/12.7d0))
        if (v.le.-63.d0) then
         tauh = .5d0/(dexp((v+46.d0)/5.d0) +
     x                  dexp((-v-238.d0)/37.5d0))
        else
         tauh = 9.5d0
        endif
           alpham_ka(i) = minf/taum
           betam_ka(i) = 1.d0 / taum - alpham_ka(i)
           alphah_ka(i) = hinf / tauh
           betah_ka(i) = 1.d0 / tauh - alphah_ka(i)

c h-current (anomalous rectifier), Huguenard & McCormick, 1992
           minf = 1.d0/(1.d0 + dexp((v+75.d0)/5.5d0))
           taum = 1.d0/(dexp(-14.6d0 -0.086d0*v) +
     x                   dexp(-1.87 + 0.07d0*v))
           alpham_ar(i) = minf / taum
           betam_ar(i) = 1.d0 / taum - alpham_ar(i)

c K2 K-current, McCormick & Huguenard
             minf = 1.d0/(1.d0 + dexp((-v-10.d0)/17.d0))
             hinf = 1.d0/(1.d0 + dexp((v+58.d0)/10.6d0))
            taum = 4.95d0 + 0.5d0/(dexp((v-81.d0)/25.6d0) +
     x                  dexp((-v-132.d0)/18.d0))
            tauh = 60.d0 + 0.5d0/(dexp((v-1.33d0)/200.d0) +
     x                  dexp((-v-130.d0)/7.1d0))
             alpham_k2(i) = minf / taum
             betam_k2(i) = 1.d0/taum - alpham_k2(i)
             alphah_k2(i) = hinf / tauh
             betah_k2(i) = 1.d0 / tauh - alphah_k2(i)

c voltage part of C-current, using 1994 kinetics, shift 60 mV
              if (v.le.-10.d0) then
       alpham_kc(i) = (2.d0/37.95d0)*dexp((v+50.d0)/11.d0 -
     x                                     (v+53.5)/27.d0)
       betam_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)-alpham_kc(i)
               else
       alpham_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)
       betam_kc(i) = 0.d0
               endif

c high-threshold gCa, from 1994, with 60 mV shift & no inactivn.
            alpham_cal(i) = 1.6d0/(1.d0+dexp(-.072d0*(v-5.d0)))
            betam_cal(i) = 0.1d0 * ((v+8.9d0)/5.d0) /
     x          (dexp((v+8.9d0)/5.d0) - 1.d0)

c M-current, from plast.f, with 60 mV shift
        alpham_km(i) = .02d0/(1.d0+dexp((-v-20.d0)/5.d0))
        betam_km(i) = .01d0 * dexp((-v-43.d0)/18.d0)

c T-current, from Destexhe, Neubig et al., 1998
         minf = 1.d0/(1.d0 + dexp((-v-56.d0)/6.2d0))
         hinf = 1.d0/(1.d0 + dexp((v+80.d0)/4.d0))
         taum = 0.204d0 + .333d0/(dexp((v+15.8d0)/18.2d0) +
     x                  dexp((-v-131.d0)/16.7d0))
          if (v.le.-81.d0) then
         tauh = 0.333 * dexp((v+466.d0)/66.6d0)
          else
         tauh = 9.32d0 + 0.333d0*dexp((-v-21.d0)/10.5d0)
          endif
              alpham_cat(i) = minf / taum
              betam_cat(i) = 1.d0/taum - alpham_cat(i)
              alphah_cat(i) = hinf / tauh
              betah_cat(i) = 1.d0 / tauh - alphah_cat(i)

1        CONTINUE

         do 2, i = 0, 639

      dalpham_naf(i) = (alpham_naf(i+1)-alpham_naf(i))/.25d0
      dbetam_naf(i) = (betam_naf(i+1)-betam_naf(i))/.25d0
      dalphah_naf(i) = (alphah_naf(i+1)-alphah_naf(i))/.25d0
      dbetah_naf(i) = (betah_naf(i+1)-betah_naf(i))/.25d0
      dalpham_kdr(i) = (alpham_kdr(i+1)-alpham_kdr(i))/.25d0
      dbetam_kdr(i) = (betam_kdr(i+1)-betam_kdr(i))/.25d0
      dalpham_ka(i) = (alpham_ka(i+1)-alpham_ka(i))/.25d0
      dbetam_ka(i) = (betam_ka(i+1)-betam_ka(i))/.25d0
      dalphah_ka(i) = (alphah_ka(i+1)-alphah_ka(i))/.25d0
      dbetah_ka(i) = (betah_ka(i+1)-betah_ka(i))/.25d0
      dalpham_k2(i) = (alpham_k2(i+1)-alpham_k2(i))/.25d0
      dbetam_k2(i) = (betam_k2(i+1)-betam_k2(i))/.25d0
      dalphah_k2(i) = (alphah_k2(i+1)-alphah_k2(i))/.25d0
      dbetah_k2(i) = (betah_k2(i+1)-betah_k2(i))/.25d0
      dalpham_km(i) = (alpham_km(i+1)-alpham_km(i))/.25d0
      dbetam_km(i) = (betam_km(i+1)-betam_km(i))/.25d0
      dalpham_kc(i) = (alpham_kc(i+1)-alpham_kc(i))/.25d0
      dbetam_kc(i) = (betam_kc(i+1)-betam_kc(i))/.25d0
      dalpham_cat(i) = (alpham_cat(i+1)-alpham_cat(i))/.25d0
      dbetam_cat(i) = (betam_cat(i+1)-betam_cat(i))/.25d0
      dalphah_cat(i) = (alphah_cat(i+1)-alphah_cat(i))/.25d0
      dbetah_cat(i) = (betah_cat(i+1)-betah_cat(i))/.25d0
      dalpham_caL(i) = (alpham_cal(i+1)-alpham_cal(i))/.25d0
      dbetam_caL(i) = (betam_cal(i+1)-betam_cal(i))/.25d0
      dalpham_ar(i) = (alpham_ar(i+1)-alpham_ar(i))/.25d0
      dbetam_ar(i) = (betam_ar(i+1)-betam_ar(i))/.25d0
2      CONTINUE
       END
