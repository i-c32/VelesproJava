package Integrals;

import java.io.*;
import java.util.*;

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

    // ✅ Static parser method belongs inside the same class
    public static BasisSet readBasisSet(File file, String name) throws IOException {
        BasisSet basisSet = new BasisSet(name);
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
                if (currentOrbital != null && parts.length >= 3) {
                    double exponent = Double.parseDouble(parts[0]);
                    double coefficient = Double.parseDouble(parts[1]);
                    double coefficient1 = Double.parseDouble(parts[2]);
                    currentOrbital.addPrimitive(exponent, coefficient);
                    currentOrbital.addPrimitive(exponent, coefficient1);
                } else if (currentOrbital != null && parts.length == 2) {
                    double exponent = Double.parseDouble(parts[0]);
                    double coefficient = Double.parseDouble(parts[1]);
                    currentOrbital.addPrimitive(exponent, coefficient);
                }
            }
        }
        return basisSet;
    }

    // Nested helper classes (can also be in their own files)
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

        public List<PrimitiveGaussian> getGaussians() { return gaussians; }

        public String getType() { return type; }
        public int getPrimitives() { return primitives; }
        public double getScale() { return scale; }
    }

    public static class PrimitiveGaussian {
        private final double exponent;
        private final double coefficient;

        public PrimitiveGaussian(double exponent, double coefficient) {
            this.exponent = exponent;
            this.coefficient = coefficient;
        }

        public double getExponent() { return exponent; }
        public double getCoefficient() { return coefficient; }
    }
}
