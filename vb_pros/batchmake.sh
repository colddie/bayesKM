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
make install
cd /home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c2
make clean
cp /home/tsun/bin/fsl/install/src/fabber_core/fabber_main.o /home/tsun/bin/fsl/install/src/fabber_core/fabber_pet_c2
make all

}