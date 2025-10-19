package org.velespro;

import Molecule.Molecule;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import com.typesafe.config.*;

import java.io.File;
import java.util.List;


public class Main {
    private static final Logger log = LoggerFactory.getLogger(Main.class);

    static void main() {
        log.info("Application started");
        // Load file
        Config config = ConfigFactory.parseFile(new File("Mis_ficheros/input.vel"));

        Molecule mol = new Molecule("water",config);


        System.out.println("Charge: " + mol.getCharge());
        System.out.println("Multiplicity: " + mol.getMultiplicity());

        System.out.println("Coordinates:");
        for (List<Object> atom : mol.getCoords()) {
            System.out.println("  " + atom.getFirst() + " -> " + atom.subList(1, 4));
        }
    }
}
