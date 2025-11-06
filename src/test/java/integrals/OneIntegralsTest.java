package integrals;

import molecule.Molecule;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.List;
import molecule.Atom;
import molecule.MoleculeIntegral;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

class OneIntegralsTest {

    private static final String TEST_PATH = "src/test/resources/";

    @Test
    void computeOverlap() throws IOException {
        // --- Mock Molecule ---
        Molecule mol = mock(Molecule.class);

        // --- Mock atoms to return ---
        Atom hydrogen1 = new Atom("H", 0.0000000, 1.43233673, -0.96104039);
        Atom hydrogen2  = new Atom("H", 0.0000000, -1.43233673, -0.96104039);

        when(mol.getAtoms()).thenReturn(List.of(hydrogen1, hydrogen2));

        BasisSet hBas = BasisSet.readBasisSet(TEST_PATH, "H_basis");

        //TODO hacer esto de manera que no llame a esta funcion
        MoleculeIntegral molInt = new MoleculeIntegral(hBas,mol);

        OneIntegrals oneIntegrals = new OneIntegrals(molInt);

        // Retrieve the overlap matrix
        double[][] ovInt = oneIntegrals.getOverlapInt();

        assertEquals(ovInt[0][1], ovInt[1][0], 1e-12);
        assertEquals(0.250987, ovInt[1][0], 1e-6);
        assertEquals(1.000000, ovInt[0][0], 1e-6);

    }

    @Test
    void computeKinetic() throws IOException {
        //TODO poner como AFTERALL
        // --- Mock Molecule ---
        Molecule mol = mock(Molecule.class);

        // --- Mock atoms to return ---
        Atom hydrogen1 = new Atom("H", 0.0000000, 1.43233673, -0.96104039);
        Atom hydrogen2  = new Atom("H", 0.0000000, -1.43233673, -0.96104039);

        when(mol.getAtoms()).thenReturn(List.of(hydrogen1, hydrogen2));

        BasisSet hBas = BasisSet.readBasisSet(TEST_PATH, "H_basis");

        //TODO hacer esto de manera que no llame a esta funcion
        MoleculeIntegral molInt = new MoleculeIntegral(hBas,mol);

        OneIntegrals oneIntegrals = new OneIntegrals(molInt);

        // Retrieve the overlap matrix
        double[][] ovInt = oneIntegrals.getKineticInt();

        assertEquals(ovInt[0][1], ovInt[1][0], 1e-12);
        assertEquals(0.0083216, ovInt[1][0], 1e-6);
        assertEquals(0.760032, ovInt[0][0], 1e-6);
    }
}