package Molecule;

import java.util.List;
import java.util.Map;

public class Atom {
    private final String symbol;
    private final double x;
    private final double y;
    private final double z;
    private final double mass;
    private final int atomicNumber;

    // Static lookup tables for atomic numbers and masses
    private static final Map<String, Integer> ATOMIC_NUMBERS = Map.<String, Integer>ofEntries(
        Map.entry("H", 1),
        Map.entry("He", 2),
        Map.entry("Li", 3),
        Map.entry("Be", 4),
        Map.entry("B", 5),
        Map.entry("C", 6),
        Map.entry("N", 7),
        Map.entry("O", 8),
        Map.entry("F", 9),
        Map.entry("Ne", 10),
        Map.entry("Na", 11),
        Map.entry("Mg", 12),
        Map.entry("Al", 13),
        Map.entry("Si", 14),
        Map.entry("P", 15),
        Map.entry("S", 16),
        Map.entry("Cl", 17),
        Map.entry("Ar", 18),
        Map.entry("K", 19),
        Map.entry("Ca", 20),
        Map.entry("Sc", 21),
        Map.entry("Ti", 22),
        Map.entry("V", 23),
        Map.entry("Cr", 24),
        Map.entry("Mn", 25),
        Map.entry("Fe", 26),
        Map.entry("Co", 27),
        Map.entry("Ni", 28),
        Map.entry("Cu", 29),
        Map.entry("Zn", 30),
        Map.entry("Ga", 31),
        Map.entry("Ge", 32),
        Map.entry("As", 33),
        Map.entry("Se", 34),
        Map.entry("Br", 35),
        Map.entry("Kr", 36),
        Map.entry("Rb", 37),
        Map.entry("Sr", 38),
        Map.entry("Y", 39),
        Map.entry("Zr", 40),
        Map.entry("Nb", 41),
        Map.entry("Mo", 42),
        Map.entry("Tc", 43),
        Map.entry("Ru", 44),
        Map.entry("Rh", 45),
        Map.entry("Pd", 46),
        Map.entry("Ag", 47),
        Map.entry("Cd", 48),
        Map.entry("In", 49),
        Map.entry("Sn", 50),
        Map.entry("Sb", 51),
        Map.entry("Te", 52),
        Map.entry("I", 53),
        Map.entry("Xe", 54),
        Map.entry("Cs", 55),
        Map.entry("Ba", 56),
        Map.entry("La", 57),
        Map.entry("Ce", 58),
        Map.entry("Pr", 59),
        Map.entry("Nd", 60),
        Map.entry("Pm", 61),
        Map.entry("Sm", 62),
        Map.entry("Eu", 63),
        Map.entry("Gd", 64),
        Map.entry("Tb", 65),
        Map.entry("Dy", 66),
        Map.entry("Ho", 67),
        Map.entry("Er", 68),
        Map.entry("Tm", 69),
        Map.entry("Yb", 70),
        Map.entry("Lu", 71),
        Map.entry("Hf", 72),
        Map.entry("Ta", 73),
        Map.entry("W", 74),
        Map.entry("Re", 75),
        Map.entry("Os", 76),
        Map.entry("Ir", 77),
        Map.entry("Pt", 78),
        Map.entry("Au", 79),
        Map.entry("Hg", 80),
        Map.entry("Tl", 81),
        Map.entry("Pb", 82),
        Map.entry("Bi", 83),
        Map.entry("Po", 84),
        Map.entry("At", 85),
        Map.entry("Rn", 86),
        Map.entry("Fr", 87),
        Map.entry("Ra", 88),
        Map.entry("Ac", 89),
        Map.entry("Th", 90),
        Map.entry("Pa", 91),
        Map.entry("U", 92),
        Map.entry("Np", 93),
        Map.entry("Pu", 94),
        Map.entry("Am", 95),
        Map.entry("Cm", 96),
        Map.entry("Bk", 97),
        Map.entry("Cf", 98),
        Map.entry("Es", 99),
        Map.entry("Fm", 100),
        Map.entry("Md", 101),
        Map.entry("No", 102),
        Map.entry("Lw", 103)
    );

    private static final Map<String, Double> ATOMIC_MASSES = Map.<String, Double>ofEntries(
            Map.entry("H", 1.00797),
            Map.entry("He", 4.0026),
            Map.entry("Li", 6.939),
            Map.entry("Be", 9.0122),
            Map.entry("B", 10.811),
            Map.entry("C", 12.01115),
            Map.entry("N", 14.0067),
            Map.entry("O", 15.9994),
            Map.entry("F", 18.9984032),
            Map.entry("Ne", 20.183),
            Map.entry("Na", 22.98976928),
            Map.entry("Mg", 24.312),
            Map.entry("Al", 26.9815386),
            Map.entry("Si", 28.0855),
            Map.entry("P", 30.9737620),
            Map.entry("S", 32.064),
            Map.entry("Cl", 35.4527),
            Map.entry("Ar", 39.948),
            Map.entry("K", 39.102),
            Map.entry("Ca", 40.08),
            Map.entry("Sc", 44.956),
            Map.entry("Ti", 47.90),
            Map.entry("V", 50.942),
            Map.entry("Cr", 51.996),
            Map.entry("Mn", 54.938),
            Map.entry("Fe", 55.847),
            Map.entry("Co", 58.933),
            Map.entry("Ni", 58.71),
            Map.entry("Cu", 63.54),
            Map.entry("Zn", 65.37),
            Map.entry("Ga", 69.72),
            Map.entry("Ge", 72.59),
            Map.entry("As", 74.922),
            Map.entry("Se", 78.96),
            Map.entry("Br", 79.909),
            Map.entry("Kr", 83.80),
            Map.entry("Rb", 85.47),
            Map.entry("Sr", 87.62),
            Map.entry("Y", 88.905),
            Map.entry("Zr", 91.22),
            Map.entry("Nb", 92.906),
            Map.entry("Mo", 95.94),
            Map.entry("Tc", 98.00),
            Map.entry("Ru", 101.07),
            Map.entry("Rh", 102.905),
            Map.entry("Pd", 106.4),
            Map.entry("Ag", 107.870),
            Map.entry("Cd", 112.40),
            Map.entry("In", 114.82),
            Map.entry("Sn", 118.69),
            Map.entry("Sb", 121.75),
            Map.entry("Te", 127.60),
            Map.entry("I", 126.904),
            Map.entry("Xe", 131.30),
            Map.entry("Cs", 132.90545),
            Map.entry("Ba", 137.34),
            Map.entry("La", 138.90547),
            Map.entry("Ce", 140.116),
            Map.entry("Pr", 140.90765),
            Map.entry("Nd", 144.242),
            Map.entry("Pm", 145.0),
            Map.entry("Sm", 150.36),
            Map.entry("Eu", 151.964),
            Map.entry("Gd", 157.25),
            Map.entry("Tb", 158.92535),
            Map.entry("Dy", 162.500),
            Map.entry("Ho", 164.93032),
            Map.entry("Er", 167.259),
            Map.entry("Tm", 168.93421),
            Map.entry("Yb", 173.054),
            Map.entry("Lu", 174.9668),
            Map.entry("Hf", 178.49),
            Map.entry("Ta", 180.948),
            Map.entry("W", 183.85),
            Map.entry("Re", 186.2),
            Map.entry("Os", 190.2),
            Map.entry("Ir", 192.2),
            Map.entry("Pt", 195.09),
            Map.entry("Au", 196.967),
            Map.entry("Hg", 200.59),
            Map.entry("Tl", 204.37),
            Map.entry("Pb", 207.19),
            Map.entry("Bi", 208.980),
            Map.entry("Po", 209.0),
            Map.entry("At", 210.0),
            Map.entry("Rn", 222.0),
            Map.entry("Fr", 223.0),
            Map.entry("Ra", 226.0),
            Map.entry("Ac", 227.0),
            Map.entry("Th", 232.038),
            Map.entry("Pa", 231.036),
            Map.entry("U", 238.029),
            Map.entry("Np", 237.0),
            Map.entry("Pu", 244.0),
            Map.entry("Am", 243.0),
            Map.entry("Cm", 247.0),
            Map.entry("Bk", 247.0),
            Map.entry("Cf", 251.0),
            Map.entry("Es", 252.0),
            Map.entry("Fm", 257.0),
            Map.entry("Md", 258.0),
            Map.entry("No", 259.0),
            Map.entry("Lw", 262.0)
    );

    public Atom(String symbol, double x, double y, double z) {
        this.symbol = symbol;
        this.x = x;
        this.y = y;
        this.z = z;
        this.atomicNumber = calculateAtomicNumber(symbol);
        this.mass = calculateMass(symbol);
    }

    // Calculation methods
    private int calculateAtomicNumber(String symbol) {
        return ATOMIC_NUMBERS.getOrDefault(symbol, 0); // 0 if unknown
    }

    private double calculateMass(String symbol) {
        return ATOMIC_MASSES.getOrDefault(symbol, 0.0); // 0.0 if unknown
    }

    // Getters
    public String getSymbol() {
        return symbol;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getZ() {
        return z;
    }

    public double getMass() {
        return mass;
    }

    public int getAtNumb() {
        return atomicNumber;
    }

    public record AtomIntegrals(String symbol, Coordinates coord,
                                String orbital, List<Double> exponent, List<Double> coefficient, AngCoordinates cartAngular) {
    }

    public record Coordinates(double x, double y, double z) {
    }

    public record AngCoordinates(int x, int y, int z) {
    }
}