
#
#       ---------------     Definitions:   -------------------------
#
# TREES simulation model
#
# This runs consistently with with and without optimization using LLVM
#
# This runs correctly when compiled under GCC with no optimization, but
# behaves differently under GCC with optimization
#

.SUFFIXES: .c .cpp


CXX    = g++
CFLAGS = -O2 -std=gnu++11

CSRC =  \
 bgc.cpp \
 calcRhat.cpp \
 calcVmat.cpp \
 data.cpp \
 data_store.cpp \
 dictionary.cpp \
 getPIL.cpp \
 getPost.cpp \
 hydraulicModel.cpp \
 parameter.cpp \
 process_functions.cpp \
 randomize.cpp \
 simulation_functions.cpp \
 simulator.cpp \
 state_store.cpp \
 time_keeper.cpp \
 trees_main.cpp \
 util.cpp

OBJ = $(CSRC:.cpp=.o)

EXE = calcRhat calcVmat getPIL getPost trees3


#      ----------------------   TOP-LEVEL TARGETS:   ------------------

all:  $(EXE)

clean:
	rm $(OBJ) $(EXE)


#      -----------------------   RULES:   -------------------------

%.o : %.mod        #  Disable "gmake"s obnoxious implicit Modula-2 rule !!

.cpp.o:
	$(CXX) $(CFLAGS) -c $<


#  ---------------------------  $(EXE) Program builds:  -----------------

calcRhat : calcRhat.o util.o
	$(CXX) $(CFLAGS)} $^ -o $@

calcVmat : calcVmat.o util.o
	$(CXX) $(CFLAGS)} $^ -o $@

getPIL : getPIL.o simulator.o simulation_functions.o process_functions.o hydraulicModel.o data_store.o state_store.o dictionary.o parameter.o time_keeper.o data.o randomize.o util.o
	$(CXX) $(CFLAGS)} $^ -o $@

getPost : getPost.o util.o
	$(CXX) $(CFLAGS)} $^ -o $@

trees3 : trees_main.o simulator.o simulation_functions.o process_functions.o hydraulicModel.o bgc.o data_store.o state_store.o dictionary.o parameter.o time_keeper.o data.o randomize.o util.o
	$(CXX) $(CFLAGS) $^ -o $@


#  ---------------------------  Dependencies:  -------------------------
# DO NOT DELETE

bgc.o			: simulator2.h
calcRhat.o              : util.h
calcVmat.o              : util.h
data.o                  : data.h
data_store.o            : data_store.h time_keeper.h dictionary.h
dictionary.o            : dictionary.h
getPIL.o                : simulator2.h state_store.h dictionary.h
getPIL.o                : data.h constants.h data_store.h time_keeper.h
getPost.o               : util.h
hydraulicModel.o        : randomize.h util.h state_store.h dictionary.h
hydraulicModel.o        : data.h constants.h data_store.h time_keeper.h
parameter.o             : parameter.h util.h
process_functions.o     : parameter.h randomize.h util.h data.h data_store.h time_keeper.h
randomize.o             : randomize.h
simulation_functions.o  : constants.h simulator2.h parameter.h randomize.h util.h
simulation_functions.o  : state_store.h dictionary.h data.h data_store.h time_keeper.h
simulator.o             : constants.h simulator2.h parameter.h randomize.h util.h
simulator.o             : state_store.h dictionary.h data.h data_store.h time_keeper.h
state_store.o           : state_store.h dictionary.h
time_keeper.o           : time_keeper.h  constants.h util.h
trees_main.o            : trees_main.h simulator2.h 
trees_main.o            : parameter.h randomize.h util.h state_store.h dictionary.h
trees_main.o            : data.h constants.h data_store.h time_keeper.h
util.o                  : util.h 


##  ---------------------------  Additional related programs:  -----------------
#
#calcpostV : calcpostV.o util.o
#	$(CXX) $(CFLAGS)} $^ -o $@
#
#dic : mcmc_nee_DIC.o simulator.o hydraulicModel.o data_store.o state_store.o dictionary.o parameter.o time_keeper.o data.o randomize.o util.o
#	$(CXX) $(CFLAGS)} $^ -o $@
#
#trees_mix : mcmc_mix.o simulator_mix.o hydraulicModel.o data_store.o state_store.o dictionary.o parameter.o time_keeper.o randomize.o util.o data.o
#	$(CXX) $(CFLAGS)} $^ -o $@
#trees2e : trees_main.o simulator_ecrit.o hydraulicModel.o data_store.o state_store.o dictionary.o parameter.o time_keeper.o randomize.o util.o data.o
#	$(CXX) $(CFLAGS)} $^ -o $@
#
#minDev : minDev.o util.o
#	$(CXX) $(CFLAGS)} $^ -o $@
#
##simulation evaluation programs
#
#eval_sims:eval_sims.o time_keeper.o data.o util.o
#	g++ -O3 eval_sims.o time_keeper.o data.o util.o -o eval_sims
#
#eval_sims.o:eval_sims.cpp
#	g++ -c -O3 eval_sims.cpp
#
##U uncertainty calcs
#
#uunc: uunc.o util.o
#	g++ -O3 uunc.o util.o -o uunc
#
#uunc.o: uunc.cpp
#	g++ -c -O3 uunc.cpp
#
##calculate support reject 0 mgs
#
#ev2mg: ev2mg.o util.o
#	g++ -O3 ev2mg.o util.o -o ev2mg
#
#ev2mg.o: ev2mg.cpp
#	g++ -c -O3 ev2mg.cpp
#
##calculate the u-uncertaity from file from ev2mg
#
#ucalc: ucalc.o util.o
#	g++ -O3 ucalc.o util.o -o ucalc
#
#ucalc.o: ucalc.cpp
#	g++ -c -O3 ucalc.cpp
#
##making restricted sets
#
#makeRset: makeRset.o util.o
#	g++ -O3 makeRset.o util.o -o makeRset
#
#makeRset.o: makeRset.cpp
#	g++ -c -O3 makeRset.cpp
#
##finding optima
#
#showMax: showMax.o util.o
#	g++ -O3 showMax.o util.o -o showMax
#
#showMax.o: showMax.cpp
#	g++ -c -O3 showMax.cpp
##constructing uncertainty interval from list
#construct_ui: construct_ui.o data_store.o data.o time_keeper.o dictionary.o util.o
#	g++ -O3 construct_ui.o data_store.o data.o time_keeper.o dictionary.o util.o -o construct_ui
#construct_ui.o: construct_ui.cpp
#	g++ -c -O3 construct_ui.cpp
#
##MCMC sims for above
#mcet_gda:mcet_gda.o simfunc.o data_store.o dictionary.o parameter.o time_keeper.o data.o randomize.o util.o
#	g++ -O3 mcet_gda.o simfunc.o data_store.o dictionary.o parameter.o time_keeper.o data.o randomize.o util.o -o mcet_gda
#
#mcet_gda.o:mcet_gda.cpp
#	g++ -c -O3 mcet_gda.cpp
#
##prediction for above mcmc posterior
#predict:predict.o simfunc.o data_store.o dictionary.o parameter.o time_keeper.o data.o randomize.o util.o
#	g++ -O3 predict.o simfunc.o data_store.o dictionary.o parameter.o time_keeper.o data.o randomize.o util.o -o predict
#predict.o:predict.cpp
#	g++ -c -O3 predict.cpp
##prediction using all params in file
#pred4:pred4.o simfunc.o data_store.o dictionary.o parameter.o time_keeper.o data.o randomize.o util.o
#	g++ -O3 pred4.o simfunc.o data_store.o dictionary.o parameter.o time_keeper.o data.o randomize.o util.o -o pred4
#pred4.o:pred4.cpp
#	g++ -c -O3 pred4.cpp
#
##grid sim for loess fit
#gridgdalo:gridgdalo.o simfunc.o data_store.o dictionary.o parameter.o time_keeper.o randomize.o util.o
#	g++ -O3 gridgdalo.o simfunc.o data_store.o dictionary.o parameter.o time_keeper.o randomize.o util.o -o gridgdalo
#
#gridgdalo.o:gridgdalo.cpp
#	g++ -c -O3 gridgdalo.cpp
#
#
