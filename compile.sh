export PETSC_DIR=$HOME/Downloads/petsc-3.25.5
export PETSC_ARCH=arch-darwin-exodus
cmake -S ~/projects/inc_ns -B ~/projects/inc_ns/build \
      -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=$PETSC_ARCH \
      -DCMAKE_C_COMPILER=$PETSC_DIR/$PETSC_ARCH/bin/mpicc \
      -DCMAKE_CXX_COMPILER=$PETSC_DIR/$PETSC_ARCH/bin/mpicxx
cmake --build ~/projects/inc_ns/build
