package Molecule;

import com.typesafe.config.Config;
import com.typesafe.config.ConfigFactory;
import org.junit.jupiter.api.Test;
import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

import static Parameter.Parameter.ANGSTROM_TO_BOHR;

public class MoleculeTest {

    @Test
    public void testMoleculeFromConfig() {
        String configStr = """
            Velespro {
              water {
                coord = [
                  {element:"H", x: 0.8668118, y: 0.6014357,  z:0.0000000}
                  {element:"H", x:-0.8668118, y: 0.6014357,  z:0.0000000}
                  {element:"O", x: 0.0000000, y:-0.0757918,  z:0.0000000}
                ]
                method = "HF"
                basis_set = "STO-3G"
                charge = 0
                multiplicity = 1
                num_atom = 3
                option = { bohr = false }
              }
            }
        """;

        Config config = ConfigFactory.parseString(configStr);
        Config velConfig = config.getConfig("Velespro");

        // Get all molecule names under Velespro
        Set<String> moleculeNames = velConfig.root().keySet();

        Molecule water = new Molecule(moleculeNames.iterator().next(),velConfig);

        // ✅ Basic checks
        assertEquals("water", water.getName());
        assertEquals(3, water.getNumAtom());
        assertEquals(0, water.getCharge());
        assertEquals(1, water.getMultiplicity());

        // ✅ Atom list check
        List<Atom> atoms = water.getAtoms();
        assertEquals(3, atoms.size());
        assertEquals("O", atoms.get(2).getSymbol());

        // ✅ Coordinate conversion check (Angstrom → Bohr)
        // 1 Å = 1.8897259886 Bohr
        double expected = 0.8668118 * ANGSTROM_TO_BOHR;
        assertEquals(expected, atoms.get(0).getX(), 1e-6);

        // ✅ Total mass check (approximate)
        assertTrue(water.getMass() > 18.0 && water.getMass() < 19.0, "Water mass should be ~18 u");

        // ✅ Distance matrix should not be null and symmetric
        double[][] matDist = water.getMatDist();
        assertNotNull(matDist);
        assertEquals(matDist[0][1], matDist[1][0], 1e-12);

        // ✅ Repulsion energy should be positive (repulsive)
        assertTrue(water.getRepEner() > 0.0, "Repulsion energy should be positive");
    }
}