/* Incorporates C preprocessor as a way of compiling the program with the appropriate physics */

/* SYSTEM: determines the system type compiled in the program
*  0= brush-melt channel, no PBC in z
 * 1= droplet(no brush on top wall) 
 * 2= charged system.No brush, 
 * 3= charged channel+brush  
 * 4= chain plus rings (RINGS set the number of rings) need defined PARTICLE_4 */

#define SYSTEM 4 

/* 0 = channel geometry: no PBC in Z; 
 * 1 = Bulk geometry PBC in 3D */

#define SYMMETRY 1   

/* Potentials and particles */

/* Wall type 
 * 1: explicit wall 
 * 2: implicit wall 9-3, 
 * 3: top: hard wall, bottom: 9-3 4=top and bottom, hard walls   (droplets)
 * 4: hard walls in top and bottom walls */

#define WALL 2

/*      Potential cut-offs for LJ
 * 3= non-additive+poor 
 * 2= non-additive+good 
 * 1= good solvent 
 * 0= poor solvent      */

#define SOLVENT 1    

/*Bending Forces */ 
#define BENDING  /* if def, the chains are assumed to be semiflexible: bending potential */
#define ORIENTATION  /* if def the chains will be oriented through an harmonic potential*//

#define PARTICLE_4 /* If defined the program runs with four different particle type */

#undef STARS /*whether you want to simulate with or without stars, sigma is fixed to 1. As well as sigma for walls */

/* Thermostat */

#define THERMOSTAT 1 /*  1=LGV 0=DPD       */

#define DPD_WEIGHT 0 /*  0=usual choice of DPD weight: Wd=(1-r/rc)^2 ; 1= constant: Wr=Wd=1 ; 2 "quartic" */ 
                     /*  wd=(1-r/rc)^4                                                                    */  
#undef DPD_CUT_OFF  2.24 /* if defined, takes this cutoff for DPD forces */
/*#define DPD_EMBEDDED*/ /* embeds DPD in fluid_fluid calculation   */
#undef DPD_VV       /* Adds a re-calculation of the Fd at the end of the iteration cycle. 
                     * Improves T=const ? Vattullainen JCP 116 2967 
                     * Good for soft potentials (not implemmented) */ 

/*
* Relaxation mechanisms when thermalizing new configurations 
*/

#define FORCE_SWITCH_ON /* switchs on the real force progressively to avoid overflows when there are strong
                       * particles overlaps */
                    /*NOTE: Use with caution because the program will not be doing real MD */    
#undef RELAX       /* When defined, velocities  are set to zero in each time step. 
                    * This is not MD, but force relaxation */


#undef POISEUILLE    /* Adds external constant force to simulate Poiseuille flow      */
#undef SHEARED      /* if defined, the shear protocols are applied, mfa_input is different!! */
                    /* NOTE: if it is not defined, wall velocities can anyway been used */

#define  STORE 1 /* 0=Writes out folded coordinates,1=writes out unfolded (to follow diffusion)   */
#undef  PROFILES      /* this uses mide for calculating profiles in the MD run.       */
#undef DIFF /*ifdef, the program will calculate diffusion coefficients. Not implemmented for all the systems  */
#undef BRUSH_SEP /* ifdef separates the density profiles of top and bottom brushes */                    
#undef FLUKT /*write out observations of the forces between Particle 4 and 3 in file  fort.555*/
#undef PINNED 2 /*applicable only with particle 4 and particle 3 and brushes,*/
                 /* 1 - fixed position, thermostat is not applied */
                 /* 2 - with the spring, thermostat is applied*/
                 /* the distance between them - from system input*/
/* Misc variables */
#define NO_WARNS     /* Does not print the warning messages of the "beads too close"                      */     
/*#define FLUID_ROUTINE 0*/  /* controls if it uses the normal fluid routine (0) or the experimental HPC-tuned (1) */
#undef FREE_HEADS        /* if defined, the heads of brushes are not fixed, epsilon between head and surface is 250, potential - attractive*/
#define BIN_TYPE  0 /*0: uses binning.f90; single counting of each interaction; 1: uses my_binning.f90, from S. Plimpton and Cem Servantie versions; double counts each interaction.  */

#define RINGS 0 /* 1 or 0 (exist or NOT_exist rings in system for SYSTEM == 4 ) for number of rings look for "n_ring" variable */

#define FIXCM 0 /* Used in SYSTEM == 4 :
                   0 - Free CM chain
                   1 - Fix CM  chain  
                */
#define CHAIN_TYPE 0 /*0: random positions in Y and Z axis ; 1: harmonic positions in Z axis */
#define CHAIN_BC 2 /* Boundary condition of the chains: 1 periodic boundary conditions; 2: Fix ends */
