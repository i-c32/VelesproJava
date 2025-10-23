package Velespro;

import Integrals.BasisSet;
import Molecule.Molecule;
import Molecule.MoleculeIntegral;
import Molecule.Atom;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import com.typesafe.config.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;


public class Main {
    private static final Logger log = LoggerFactory.getLogger(Main.class);

    static void main(String[] args) throws IOException {
        log.info("Application started");
        // Load file
        Config config = ConfigFactory.parseFile(new File(args[0]));
        Config velConfig = config.getConfig("Velespro");

        // Get all molecule names under Velespro
        Set<String> moleculeNames = velConfig.root().keySet();

        Molecule mol = new Molecule(moleculeNames.iterator().next(),velConfig);

        File outputFile = new File("Mis_ficheros/molecule_info.txt");

        // ✅ Append mode = true
        try (PrintWriter writer = new PrintWriter(new FileWriter(outputFile, true))) {
            writer.println("Molecule: " + mol.getName());
            writer.println("Charge: " + mol.getCharge());
            writer.println("Multiplicity: " + mol.getMultiplicity());
            writer.println("Atoms: " + mol.getNumAtom());
            writer.println("Mass: " + mol.getMass());
            writer.println("Erep: " + mol.getRepEner());
            writer.println("Coordinates:");
            for (Atom atom : mol.getAtoms()) {
                writer.printf("%-3s %15.6f %15.6f %15.6f%n",
                        atom.getSymbol(), atom.getX(), atom.getY(), atom.getZ());
            }
            writer.println(); // blank line between entries
        } catch (IOException e) {
            log.error("Failed to write molecule coordinates to file", e);        }

        String bas = velConfig.getString("water.basis set");
        BasisSet sto3g = BasisSet.readBasisSet(bas);

        System.out.println("Loaded basis: " + sto3g.getName());
        System.out.println("Elements: " + sto3g.getOrbitals("O").size() + " orbitals for Oxygen");

        MoleculeIntegral molInt = new MoleculeIntegral(sto3g,mol);

        for (Atom.AtomIntegrals atom : molInt.getAtoms()) {
            System.out.printf("%-3s %15.6f %15.6f %15.6f%n",
                    atom.symbol(), atom.x(), atom.y(), atom.z());
            System.out.printf("     Orbital: %-4s%n", atom.orbital());
            System.out.printf("     Exponents:   %s%n", atom.exponent());
            System.out.printf("     Coefficients:%s%n", atom.coefficient());
            System.out.printf("     Cart. Ang.:  %s%n%n", atom.cartAngular());
        }
        log.info("Application ended");
    }
}
