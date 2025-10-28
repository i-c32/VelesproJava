package Molecule;

import Integrals.BasisSet;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.mockito.Mockito.*;
import static org.junit.jupiter.api.Assertions.*;

class MoleculeIntegralTest {

    @Test
    void testMoleculeIntegralConstruction() {
        // --- Mock Molecule ---
        Molecule mol = mock(Molecule.class);

        // --- Mock atoms to return ---
        Atom hydrogen = new Atom("H", 0.0, 0.0, 0.0);
        Atom oxygen   = new Atom("O", 0.0, 0.0, 1.0);

        when(mol.getAtoms()).thenReturn(List.of(hydrogen, oxygen));

        // --- Mock BasisSet ---
        BasisSet basisSet = mock(BasisSet.class);

        // --- Create mock orbitals ---
        BasisSet.Orbital sOrbital = new BasisSet.Orbital("S", 1, 1.0);
        sOrbital.addPrimitive(1.0, 0.5);

        BasisSet.Orbital pOrbital = new BasisSet.Orbital("P", 1, 1.0);
        pOrbital.addPrimitive(0.5, 0.0, 1.0);

        BasisSet.Orbital spOrbital = new BasisSet.Orbital("SP", 1, 1.0);
        spOrbital.addPrimitive(1.0, 0.6, 0.4);

        // --- Mock getOrbitals to return the orbitals for any atom symbol ---
        when(basisSet.getOrbitals(anyString())).thenReturn(List.of(sOrbital, pOrbital, spOrbital));

        // --- Construct MoleculeIntegral ---
        MoleculeIntegral molIntegral = new MoleculeIntegral(basisSet, mol);

        // --- Assertions ---
        List<Atom.AtomIntegrals> mInt = molIntegral.getAtoms();
        assertNotNull(mInt);
        assertFalse(mInt.isEmpty());
        assertTrue(mInt.stream().anyMatch(i -> i.orbital().equals("S")), "Should contain S orbitals");
        assertTrue(mInt.stream().anyMatch(i -> i.orbital().equals("Px")), "Should contain Px orbitals");
    }
}