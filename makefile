include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

helloworld: helloworld.o chkopts
	-${CLINKER} -o helloworld helloworld.o ${PETSC_LIB}
		${RM} helloworld.o