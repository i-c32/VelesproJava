package Molecule;

import com.typesafe.config.Config;

import java.util.ArrayList;
import java.util.List;

import static Parameter.Parameter.ANGSTROM_TO_BOHR;

public class Molecule {
    private final String name;
    private final double mass;
    private final int numAtom;
    private final int charge;
    private final int multiplicity;
    private final List<Atom> atoms;
    private final double[][] matDist;
    private final double repEnergy;

    public Molecule(String name, Config config) {
        this.name = name;

        Config molec = config.getConfig(name);
        this.charge = molec.getInt("charge");
        this.multiplicity = molec.getInt("multiplicity");
        this.numAtom = molec.getInt("num_atom");

        List<Atom> atomsList = new ArrayList<>();

        for (Config atom : molec.getConfigList("coord")) {
                String element = atom.getString("element");
                double x = atom.getDouble("x");
                double y = atom.getDouble("y");
                double z = atom.getDouble("z");

                // Convert if in Angstroms
                if (!molec.getBoolean("option.bohr")) {
                    x *= ANGSTROM_TO_BOHR;
                    y *= ANGSTROM_TO_BOHR;
                    z *= ANGSTROM_TO_BOHR;
                }

                atomsList.add(new Atom(element, x, y, z));
        }
        this.atoms = atomsList;

        this.mass = this.atoms.stream()
                .mapToDouble(Atom::getMass)
                .sum();

        this.matDist = getDistanceMatrix(this.numAtom);

        this.repEnergy = getRepulsionEnergy(this.numAtom);
    }

    private double[][] getDistanceMatrix(int n) {
        double[][] distances = new double[n][n];

        for (int i = 0; i < n; i++) {
            Atom a1 = atoms.get(i);
            for (int j = i; j < n; j++) { // use symmetry
                Atom a2 = atoms.get(j);

                double dx = a1.getX() - a2.getX();
                double dy = a1.getY() - a2.getY();
                double dz = a1.getZ() - a2.getZ();

                double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                distances[i][j] = dist;
                distances[j][i] = dist; // symmetric
            }
        }
        return distances;
    }

    private double getRepulsionEnergy(int n) {
        double energy = 0.0;
        double[][] dist = getDistanceMatrix(n);

        for (int i = 0; i < n; i++) {
            Atom a1 = atoms.get(i);
            for (int j = i + 1; j < n; j++) {
                Atom a2 = atoms.get(j);
                double r = dist[i][j];
                if (r > 1e-12) { // avoid division by zero
                    // Simple approximation: C_ij proportional to atomic numbers
                    energy += a1.getAtNumb() * a2.getAtNumb() / r; // repulsive term
                }
            }
        }
        return energy;
    }

    // ✅ Getters
    public String getName() { return name; }
    public int getNumAtom() { return numAtom; }
    public int getCharge() { return charge; }
    public int getMultiplicity() { return multiplicity; }
    public List<Atom> getAtoms() { return atoms; }
    public double getMass() { return mass; }
    public double[][] getMatDist() { return matDist; }
    public double getRepEner() { return repEnergy; }
}
