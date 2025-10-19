package org.velespro;

import Molecule.Molecule;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import com.typesafe.config.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


public class Main {
    private static final Logger log = LoggerFactory.getLogger(Main.class);

    static void main() {
        log.info("Application started");
        // Load file
        Config config = ConfigFactory.parseFile(new File("Mis_ficheros/input.vel"));

        Molecule mol = new Molecule("water",config);

        File outputFile = new File("Mis_ficheros/molecule_info.txt");

        // ✅ Append mode = true
        try (PrintWriter writer = new PrintWriter(new FileWriter(outputFile, true))) {
            writer.println("Molecule: " + mol.getName());
            writer.println("Charge: " + mol.getCharge());
            writer.println("Multiplicity: " + mol.getMultiplicity());
            writer.println("Atoms: " + mol.getNumAtom());
            writer.println("Coordinates:");
            for (Molecule.Atom atom : mol.getCoords()) {
                writer.printf("%-3s %15.6f %15.6f %15.6f%n",
                        atom.symbol(), atom.x(), atom.y(), atom.z());
            }
            writer.println(); // blank line between entries
        } catch (IOException e) {
            log.error("Failed to write molecule coordinates to file", e);        }
    }
}
