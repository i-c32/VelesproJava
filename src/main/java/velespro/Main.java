package velespro;

import integrals.BasisSet;
import integrals.OneIntegrals;
import molecule.Molecule;
import molecule.MoleculeIntegral;
import molecule.Atom;
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
        Config velConfig = config.getConfig("velespro");

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

        MoleculeIntegral molInt = new MoleculeIntegral(sto3g,mol);

        // Build overlap integrals object
        OneIntegrals oneIntegrals = new OneIntegrals(molInt);

        // Retrieve the overlap matrix
        double[][] test = oneIntegrals.getKineticInt();


        for (int i = 0; i < test.length; i++) {
            for (int j = 0; j < test[i].length; j++) {
                System.out.printf("%12.9f ", test[i][j]);
            }
            System.out.println();
        }

        log.info("Application ended");
    }
}
