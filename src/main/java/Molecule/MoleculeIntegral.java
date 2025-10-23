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

                // Build integrals depending on the orbital type
                switch (orbType) {
                    case "S" -> atomsList.add(new Atom.AtomIntegrals(
                            atSymb, atom.getX(), atom.getY(), atom.getZ(),
                            "S", exponents, coeffS, Arrays.asList(0, 0, 0)));

                    case "P" -> {
                        atomsList.add(new Atom.AtomIntegrals(atSymb, atom.getX(), atom.getY(), atom.getZ(),
                                "Px", exponents, coeffP, Arrays.asList(1, 0, 0)));
                        atomsList.add(new Atom.AtomIntegrals(atSymb, atom.getX(), atom.getY(), atom.getZ(),
                                "Py", exponents, coeffP, Arrays.asList(0, 1, 0)));
                        atomsList.add(new Atom.AtomIntegrals(atSymb, atom.getX(), atom.getY(), atom.getZ(),
                                "Pz", exponents, coeffP, Arrays.asList(0, 0, 1)));
                    }

                    case "SP" -> {
                        // Hybrid orbital: S and P share exponents but have distinct coefficients
                        atomsList.add(new Atom.AtomIntegrals(atSymb, atom.getX(), atom.getY(), atom.getZ(),
                                "S", exponents, coeffS, Arrays.asList(0, 0, 0)));
                        atomsList.add(new Atom.AtomIntegrals(atSymb, atom.getX(), atom.getY(), atom.getZ(),
                                "Px", exponents, coeffP, Arrays.asList(1, 0, 0)));
                        atomsList.add(new Atom.AtomIntegrals(atSymb, atom.getX(), atom.getY(), atom.getZ(),
                                "Py", exponents, coeffP, Arrays.asList(0, 1, 0)));
                        atomsList.add(new Atom.AtomIntegrals(atSymb, atom.getX(), atom.getY(), atom.getZ(),
                                "Pz", exponents, coeffP, Arrays.asList(0, 0, 1)));
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
