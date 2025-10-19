package Molecule;

import com.typesafe.config.Config;

import java.util.List;

public class Molecule {
    private String name;
    private int numAtom;
    private int charge;
    private int multiplicity;
    private List<List<Object>> coords;

    public Molecule(String name, Config config) {
        this.name = name;

        Config molec = config.getConfig("Velespro."+name);
        this.charge = molec.getInt("charge");
        this.multiplicity = molec.getInt("multiplicity");
        this.numAtom = molec.getInt("num_atom");;
        this.coords = (List<List<Object>>) (Object) molec.getAnyRefList("coord");
    }

    // ✅ Getters
    public String getName() { return name; }
    public int getNumAtom() { return numAtom; }
    public int getCharge() { return charge; }
    public int getMultiplicity() { return multiplicity; }
    public List<List<Object>> getCoords() { return coords; }
}
