#!/bin/bash

curdir="$(pwd)"

#---- online install ----#

#-> 1. compile gsamples
echo ""
echo "==========================="
echo "Section I: compile gsamples"
echo "==========================="
echo ""
cd $curdir
cd gsamples_package
#	./install_gsamples.sh -v 0
#	mv gsamples ../../bin
cd ../

#-> 2. compile pinpoint
echo ""
echo "============================"
echo "Section II: compile pinpoint"
echo "============================"
echo ""
cd $curdir
cd pinpoint_package
#	./install_pinpoint.sh -v 0
#	mv pinpoint ../../bin
cd ../

#-> 3. compile DeepCNF_Pred
echo ""
echo "========================="
echo "Section III: DeepCNF_Pred"
echo "========================="
echo ""
cd $curdir
cd DeepCNF_AUC/source_code/DeepCNF_Pred_src
	make
	mv DeepCNF_Pred ../../../../bin
cd ../../../

#-> 4. compile others
echo ""
echo "=========================="
echo "Section IV: compile others"
echo "=========================="
echo ""
cd $curdir
make

