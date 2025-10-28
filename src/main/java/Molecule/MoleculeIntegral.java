package Molecule;

import Integrals.BasisSet;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MoleculeIntegral {
    private static final Logger log = LoggerFactory.getLogger(MoleculeIntegral.class);

    private final List<Atom.AtomIntegrals> atoms;

    public MoleculeIntegral(BasisSet basisSet, Molecule mol) {

        List<TypeOrbital> orbS = new ArrayList<>();
        List<TypeOrbital> orbP = new ArrayList<>();
        List<TypeOrbital> orbSP = new ArrayList<>();

        // Contructor of the orbitals
        orbS.add(new TypeOrbital("S",new Atom.AngCoordinates(0, 0, 0)));

        orbP.add(new TypeOrbital("Px",new Atom.AngCoordinates(1, 0, 0)));
        orbP.add(new TypeOrbital("Py",new Atom.AngCoordinates(0, 1, 0)));
        orbP.add(new TypeOrbital("Pz",new Atom.AngCoordinates(0, 0, 1)));

        orbSP.add(new TypeOrbital("S",new Atom.AngCoordinates(0, 0, 0)));
        orbSP.add(new TypeOrbital("Px",new Atom.AngCoordinates(1, 0, 0)));
        orbSP.add(new TypeOrbital("Py",new Atom.AngCoordinates(0, 1, 0)));
        orbSP.add(new TypeOrbital("Pz",new Atom.AngCoordinates(0, 0, 1)));

        List<Atom.AtomIntegrals> atomsList = new ArrayList<>();

        for (Atom atom : mol.getAtoms()) {
            String atSymb = atom.getSymbol();
            for (BasisSet.Orbital orbit : basisSet.getOrbitals(atSymb)) {
                String orbType = orbit.getType();

                List<Double> exponents = new ArrayList<>();
                List<Double> coeff1 = new ArrayList<>();
                List<Double> coeff2 = new ArrayList<>();

                for (BasisSet.PrimitiveGaussian g : orbit.getGaussians()) {
                    exponents.add(g.exponent());
                    coeff1.add(g.coeffS());
                    coeff2.add(g.coeffP());
                }

                Atom.Coordinates coord = new Atom.Coordinates(atom.getX(), atom.getY(), atom.getZ());

                // Build integrals depending on the orbital type
                switch (orbType) {
                    case "S" -> {
                        for (TypeOrbital o : orbS) {
                            atomsList.add(new Atom.AtomIntegrals(
                                    atSymb, coord, o.type(), exponents, coeff1, o.corAng()));
                        }
                    }

                    case "P" -> {
                        for (TypeOrbital o : orbP) {
                            atomsList.add(new Atom.AtomIntegrals(
                                    atSymb, coord, o.type(), exponents, coeff1, o.corAng()));
                        }
                    }

                    case "SP" -> {
                        for (TypeOrbital o : orbSP) {
                            if (Objects.equals(o.type(), "S")) {
                                atomsList.add(new Atom.AtomIntegrals(
                                        atSymb, coord, o.type(), exponents, coeff1, o.corAng()));
                            }
                            else {
                                atomsList.add(new Atom.AtomIntegrals(
                                        atSymb, coord, o.type(), exponents, coeff2, o.corAng()));
                            }
                        }
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

    public record TypeOrbital(String type, Atom.AngCoordinates corAng){}
}
