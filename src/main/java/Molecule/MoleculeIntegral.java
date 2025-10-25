package Molecule;

import Integrals.BasisSet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MoleculeIntegral {
    private final List<Atom.AtomIntegrals> atoms;

    public MoleculeIntegral(BasisSet basisSet, Molecule mol) {

        List<Atom.AtomIntegrals> atomsList = new ArrayList<>();

        for (Atom atom : mol.getAtoms()) {
            String atSymb = atom.getSymbol();
            for (BasisSet.Orbital orbit : basisSet.getOrbitals(atSymb)) {
                String orbType = orbit.getType();

                List<Double> exponents = new ArrayList<>();
                List<Double> coeffS = new ArrayList<>();
                List<Double> coeffP = new ArrayList<>();

                for (BasisSet.PrimitiveGaussian g : orbit.getGaussians()) {
                    exponents.add(g.exponent());
                    coeffS.add(g.coeffS());
                    coeffP.add(g.coeffP());
                }

                Atom.coordinates coord = new Atom.coordinates(atom.getX(), atom.getY(), atom.getZ());
                Atom.angCoordinates cors = new Atom.angCoordinates(0, 0, 0);
                Atom.angCoordinates corpx = new Atom.angCoordinates(1, 0, 0);
                Atom.angCoordinates corpy = new Atom.angCoordinates(0, 1, 0);
                Atom.angCoordinates corpz = new Atom.angCoordinates(0, 0, 1);

                // Build integrals depending on the orbital type
                switch (orbType) {
                    case "S" -> atomsList.add(new Atom.AtomIntegrals(
                            atSymb, coord, "S", exponents, coeffS, cors));

                    case "P" -> {
                        atomsList.add(new Atom.AtomIntegrals(atSymb, coord,
                                "Px", exponents, coeffP, corpx));
                        atomsList.add(new Atom.AtomIntegrals(atSymb, coord,
                                "Py", exponents, coeffP, corpy));
                        atomsList.add(new Atom.AtomIntegrals(atSymb, coord,
                                "Pz", exponents, coeffP, corpz));
                    }

                    case "SP" -> {
                        // Hybrid orbital: S and P share exponents but have distinct coefficients
                        atomsList.add(new Atom.AtomIntegrals(atSymb, coord,
                                "S", exponents, coeffS, cors));
                        atomsList.add(new Atom.AtomIntegrals(atSymb, coord,
                                "Px", exponents, coeffP, corpx));
                        atomsList.add(new Atom.AtomIntegrals(atSymb, coord,
                                "Py", exponents, coeffP, corpy));
                        atomsList.add(new Atom.AtomIntegrals(atSymb, coord,
                                "Pz", exponents, coeffP, corpz));
                    }

                    default -> {
                        // Handle D/F/G orbitals later if needed
                        System.err.println("Unsupported orbital type: " + orbType);
                    }
                }
            }
        }
        this.atoms = atomsList;
    }

    public List<Atom.AtomIntegrals> getAtoms() {
        return atoms;
    }
}
