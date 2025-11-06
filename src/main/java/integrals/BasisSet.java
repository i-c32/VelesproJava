package integrals;

import java.io.*;
import java.util.*;

import static parameter.Parameter.PATH;

public class BasisSet {
    private String name;
    private Map<String, List<Orbital>> elementBasis = new HashMap<>();

    public BasisSet(String name) {
        this.name = name;
    }

    public void addOrbital(String element, Orbital orbital) {
        elementBasis.computeIfAbsent(element, k -> new ArrayList<>()).add(orbital);
    }

    public List<Orbital> getOrbitals(String element) {
        return elementBasis.getOrDefault(element, List.of());
    }

    public String getName() { return name; }


    // Original for production use
    public static BasisSet readBasisSet(String name) throws IOException {
        return readBasisSet(PATH, name);
    }

    // Overloaded version for testing or alternative paths
    public static BasisSet readBasisSet(String directory, String name) throws IOException {
        BasisSet basisSet = new BasisSet(name);
        File file = new File(directory + name);
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            String currentElement = null;
            Orbital currentOrbital = null;

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                if (line.equals("****")) {
                    currentElement = null;
                    continue;
                }

                String[] parts = line.split("\\s+");

                // Element definition line
                if (parts.length == 2 && parts[1].equals("0")) {
                    currentElement = parts[0];
                    continue;
                }

                // Orbital definition line
                if (parts.length >= 3 && Character.isLetter(parts[0].charAt(0))) {
                    currentOrbital = new Orbital(parts[0], Integer.parseInt(parts[1]), Double.parseDouble(parts[2]));
                    basisSet.addOrbital(currentElement, currentOrbital);
                    continue;
                }

                // Gaussian primitive line
                if (currentOrbital != null) {
                    double exponent = Double.parseDouble(parts[0]);
                    if (currentOrbital.getType().equalsIgnoreCase("SP")) {
                        // SP: 3 columns — exponent, coeff_S, coeff_P
                        double coeffS = Double.parseDouble(parts[1]);
                        double coeffP = Double.parseDouble(parts[2]);
                        currentOrbital.addPrimitive(exponent, coeffS, coeffP);
                    } else {
                        // S/P/D/F: 2 columns — exponent, coeff
                        double coeff = Double.parseDouble(parts[1]);
                        currentOrbital.addPrimitive(exponent, coeff);
                    }
                }
            }
        }
        return basisSet;
    }

    // Nested helper classes
    public static class Orbital {
        private final String type; // S, P, D, SP, etc.
        private final int primitives;
        private final double scale;
        private final List<PrimitiveGaussian> gaussians = new ArrayList<>();

        public Orbital(String type, int primitives, double scale) {
            this.type = type;
            this.primitives = primitives;
            this.scale = scale;
        }

        public void addPrimitive(double exponent, double coefficient) {
            gaussians.add(new PrimitiveGaussian(exponent, coefficient));
        }

        // Overload for SP orbitals
        public void addPrimitive(double exponent, double coeffS, double coeffP) {
            gaussians.add(new PrimitiveGaussian(exponent, coeffS, coeffP));
        }

        public List<PrimitiveGaussian> getGaussians() { return gaussians; }

        public String getType() { return type; }
        public int getPrimitives() { return primitives; }
        public double getScale() { return scale; }
    }

    public record PrimitiveGaussian(double exponent, double coeffS, double coeffP) {
        // Constructor for S/P/D
        public PrimitiveGaussian(double exponent, double coefficient) {
            this(exponent, coefficient, 0.0);
        }
    }
}
