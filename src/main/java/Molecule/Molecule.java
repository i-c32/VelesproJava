package Molecule;

import com.typesafe.config.Config;

import java.util.ArrayList;
import java.util.List;

import static Parameter.Parameter.ANGSTROM_TO_BOHR;

public class Molecule {
    private final String name;
    private final int numAtom;
    private final int charge;
    private final int multiplicity;
    public record Atom(String symbol, double x, double y, double z) {}
    private final List<Atom> coords;

    public Molecule(String name, Config config) {
        this.name = name;

        Config molec = config.getConfig("Velespro."+name);
        this.charge = molec.getInt("charge");
        this.multiplicity = molec.getInt("multiplicity");
        this.numAtom = molec.getInt("num_atom");
        List<?> rawList = molec.getAnyRefList("coord");

        List<Atom> atoms = new ArrayList<>();

        for (Object o : rawList) {
            if (o instanceof List<?> inner) {
                String symbol = (String) inner.get(0);
                double x = ((Number) inner.get(1)).doubleValue();
                double y = ((Number) inner.get(2)).doubleValue();
                double z = ((Number) inner.get(3)).doubleValue();

                // ✅ Convert if in Angstroms
                if (!molec.getBoolean("option.bohr")) {
                    x *= ANGSTROM_TO_BOHR;
                    y *= ANGSTROM_TO_BOHR;
                    z *= ANGSTROM_TO_BOHR;
                }

                atoms.add(new Atom(symbol, x, y, z));
            }
        }
        this.coords = atoms;
    }

    // ✅ Getters
    public String getName() { return name; }
    public int getNumAtom() { return numAtom; }
    public int getCharge() { return charge; }
    public int getMultiplicity() { return multiplicity; }
    public List<Atom> getCoords() { return coords; }
}
