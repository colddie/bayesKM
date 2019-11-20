#!/bin/bash
#
# this script should not be run directly,
# instead you need to source it from your .bashrc,
# by adding this line:
#   . ~/bin/batchmake.sh
#  
# make sure run sed -i 's/\r//' batchmake.sh first
function batchmake() {
cd /home/tsun/bin/fsl/install/src/fabber_core
make cleana
make install -j4
cd /home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c2
make clean
cp /home/tsun/bin/fsl/install/src/fabber_core/fabber_main.o /home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c2
make all -j4
cd /home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1
make clean
cp /home/tsun/bin/fsl/install/src/fabber_core/fabber_main.o /home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c1
make all -j4
cp /home/tsun/bin/fsl/install/src/fabber_core/fabber_main.o /home/tsun/bin/fsl/install/src/fabber_core/fabber_linear
make all -j4
}