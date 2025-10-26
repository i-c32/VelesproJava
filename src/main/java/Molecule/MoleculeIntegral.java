package Molecule;

import Integrals.BasisSet;

import java.util.ArrayList;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MoleculeIntegral {
    private static final Logger log = LoggerFactory.getLogger(MoleculeIntegral.class);

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

                Atom.Coordinates coord = new Atom.Coordinates(atom.getX(), atom.getY(), atom.getZ());
                Atom.AngCoordinates cors = new Atom.AngCoordinates(0, 0, 0);
                Atom.AngCoordinates corpx = new Atom.AngCoordinates(1, 0, 0);
                Atom.AngCoordinates corpy = new Atom.AngCoordinates(0, 1, 0);
                Atom.AngCoordinates corpz = new Atom.AngCoordinates(0, 0, 1);

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

                    default -> log.error("Unsupported orbital type: {}", orbType);

                }
            }
        }
        this.atoms = atomsList;
    }

    public List<Atom.AtomIntegrals> getAtoms() {
        return atoms;
    }
}
