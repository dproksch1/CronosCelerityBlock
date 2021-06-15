orbitTest: orbitTest.C ../Orbit.H
	gcc $< -lm -lstdc++ -o $@
