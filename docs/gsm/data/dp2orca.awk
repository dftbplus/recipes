#!/usr/bin/awk -f
# SPDX-Identifier: Unlicense

# Produce expected output for GSM driver by writing Orca standard output.
END { orca_output(iat, mermin_energy, forces) }

# Write energy and gradient of the system in expected output format.
# This routine emulates the necessary parts of an Orca standard output.
#
# Input
# -----
#   nat: Number of atoms of the whole system
#   energy: Mermin free energy
#   forces: Atomic forces for each atom
#
# Output
# ------
#   Writes Orca formatted information to standard output
function orca_output(nat, energy, forces) {
  printf "ORCA-Dummy output for GSM TS Optimizer. Not a real ORCA-Output\n"
  printf "Total Energy       :  %20.14f\n", energy
  printf "------------------\n"
  printf "CARTESIAN GRADIENT\n"
  printf "------------------\n"
  printf "\n"
  for (jat = 1; jat <= nat; jat++) {
    printf "%4d%4s   :%15.9f%15.9f%15.9f\n", jat - 1, "X",
           -forces[jat][1], -forces[jat][2], -forces[jat][3]
  }
}

# Collect the energy and reset the flag.
read_energy > 0 {
  mermin_energy = $1
  read_energy--
}

# Collect a line from the force output.
read_forces > 0 {
  iat++
  for (ic = 1; ic <= NF; ic++) {
    forces[iat][ic] = $ic
    read_forces--
  }
}

# Label for total energy found.
/^mermin_energy/ { read_energy = 1 }

# Label for gradient found.
/^forces/ { gsub(/[:,]/, " ") }
/^forces/ {
  read_forces = 1
  for (dim = 1; dim <= $3; dim++) {
    read_forces *= $(3+dim)
  }
}
