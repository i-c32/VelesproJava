package Molecule;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class AtomTest {

    @Test
    public void testAtomInitialization() {
        // Create a simple atom (Hydrogen)
        Atom h = new Atom("H", 0.0, 0.1, 0.2);

        // Check symbol
        assertEquals("H", h.getSymbol());

        // Check coordinates
        assertEquals(0.0, h.getX(), 1e-8);
        assertEquals(0.1, h.getY(), 1e-8);
        assertEquals(0.2, h.getZ(), 1e-8);

        // Check atomic number (from your ATOMIC_NUMBERS map)
        assertEquals(1, h.getAtNumb());

        // Check mass (from your ATOMIC_MASSES map)
        assertTrue(h.getMass() > 0.0, "Mass should be positive for known atoms");
    }

    @Test
    public void testUnknownAtom() {
        Atom x = new Atom("Xx", 1.0, 2.0, 3.0);

        // Unknown elements should return default 0 or 0.0
        assertEquals(0, x.getAtNumb());
        assertEquals(0.0, x.getMass());
    }
}