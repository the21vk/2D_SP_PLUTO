#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  CHARACTERISTIC_TRACING
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            4

/* -- physics dependent declarations -- */

#define  EOS                            ISOTHERMAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   DIV_CLEANING
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    EXPLICIT
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      EXPLICIT
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  ETA                            0
#define  WIDTH                          1
#define  PSI0                           2
#define  BETA                           3

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        MC_LIM
#define  CT_EMF_AVERAGE                 CT_CONTACT
#define  USE_RANDOM_PERTURBATION        NO

/* [End] user-defined constants (do not change this line) */
