#FLAGS =   -w
#FLAGS = -c  -w
#FLAGS = -g -d1 -C -w -save
#FLAGS = -O3 -qextchk -qarch=auto -bmaxdata:0x80000000
#FLAGS = -g -d1 -C -w  -bmaxdata:0x80000000


STRUCT = groucho_gapbld.o groucho_gapbld_mix.o  gettime.o otis_table_setup.o dexptablesmall_setup.o dexptablebig_setup.o synaptic_map_construct.o sy        naptic_compmap_construct.o durand.o

INTEGRATE = integrate_tuftIBVx3B.o fnmda.o integrate_deepLTSx.o integrate_deepaxaxx.o integrate_deepbaskx.o integrate_deepng.o  integrate_nrtxB.o in        tegrate_supLTSX.o integrate_supaxaxx.o integrate_supbaskx.o integrate_supng.o integrate_suppyrFRBxPB.o  integrate_tcrxB.o  integrate_tuftRSXXB.o inte        grate_spinstelldiegoxB.o


alphaY : alphaY.f makefile
         mpxlf $(FLAGS) alphaY.f $(STRUCT) $(INTEGRATE) $(TRACE) -o alphaY

clean :
        rm -f alphaY
