package integrals;

import molecule.Atom;
import molecule.MoleculeIntegral;

import java.util.List;

public class OneIntegrals {
    private final double[][] overlapInt;
    private final double[][] kineticInt;
    private final double[][] potentialInt;

    // Constructor — builds overlap integral matrix automatically
    public OneIntegrals(MoleculeIntegral molInt) {
        this.overlapInt = computeOverlap(molInt);
        this.kineticInt = computeKinetic(molInt);
        this.potentialInt = computePotential(molInt);
    }

    // Getter integrals matrix
    public double[][] getOverlapInt() {
        return overlapInt;
    }
    public double[][] getKineticInt() {
        return kineticInt;
    }
    public double[][] getPotentialInt() {
        return potentialInt;
    }

    public static double[][] computeOverlap(MoleculeIntegral molInt) {

        List<Atom.AtomIntegrals> atomInt = molInt.getAtoms();
        int nOrb = atomInt.size();

        double[][] intOvelap = new double[nOrb][nOrb];

        for (int i = 0; i < nOrb; i++) {
            for (int j = i; j < nOrb; j++) {
                intOvelap[i][j] = overlap(atomInt.get(i),atomInt.get(j));
                intOvelap[j][i] = intOvelap[i][j];
            }
        }
        return intOvelap;
    }

    public static double[][] computeKinetic(MoleculeIntegral molInt) {

        List<Atom.AtomIntegrals> atomInt = molInt.getAtoms();
        int nOrb = atomInt.size();

        double[][] intKinetic = new double[nOrb][nOrb];

        for (int i = 0; i < nOrb; i++) {
            for (int j = i; j < nOrb; j++) {
                intKinetic[i][j] = kinetic(atomInt.get(i),atomInt.get(j));
                intKinetic[j][i] = intKinetic[i][j];
            }
        }
        return intKinetic;
    }

    public static double[][] computePotential(MoleculeIntegral molInt) {

        List<Atom.AtomIntegrals> atomInt = molInt.getAtoms();
        List<Integer> numAtomic = molInt.getAtNumbers();
        List<Atom.Coordinates> nucCoord = molInt.getNucCoord();
        int nOrb = atomInt.size();

        double[][] intPotential = new double[nOrb][nOrb];

        for (int i = 0; i < nOrb; i++) {
            for (int j = i; j < nOrb; j++) {
                intPotential[i][j] = potential(atomInt.get(i), atomInt.get(j), numAtomic, nucCoord);
                intPotential[j][i] = intPotential[i][j];
            }
        }
        return intPotential;
    }

    private static double overlap(Atom.AtomIntegrals atomInt1,
                                 Atom.AtomIntegrals atomInt2) {

        double resultado = 0.0;

        for (int i = 0; i < atomInt1.exponent().size(); i++) {
            for (int j = 0; j < atomInt2.exponent().size(); j++) {
                resultado += normOrbital(atomInt1.exponent().get(i), atomInt1.cartAngular())
                        * normOrbital(atomInt2.exponent().get(j), atomInt2.cartAngular())
                        * atomInt1.coefficient().get(i) * atomInt2.coefficient().get(j)
                        * ovInt(atomInt1.exponent().get(i), atomInt2.exponent().get(j), atomInt1.coord(),
                        atomInt2.coord(), atomInt1.cartAngular(), atomInt2.cartAngular());
            }
        }

        return resultado;
    }

    private static double kinetic(Atom.AtomIntegrals atomInt1,
                                 Atom.AtomIntegrals atomInt2) {

        double resultado = 0.0;

        for (int i = 0; i < atomInt1.exponent().size(); i++) {
            for (int j = 0; j < atomInt2.exponent().size(); j++) {
                resultado += normOrbital(atomInt1.exponent().get(i), atomInt1.cartAngular())
                        * normOrbital(atomInt2.exponent().get(j), atomInt2.cartAngular())
                        * atomInt1.coefficient().get(i) * atomInt2.coefficient().get(j)
                        * kinInt(atomInt1.exponent().get(i), atomInt2.exponent().get(j), atomInt1.coord(),
                        atomInt2.coord(), atomInt1.cartAngular(), atomInt2.cartAngular());
            }
        }

        return resultado;
    }

    private static double potential(Atom.AtomIntegrals atomInt1,
                                  Atom.AtomIntegrals atomInt2, List<Integer> numAtomic, List<Atom.Coordinates> nucCoord) {

        double resultado = 0.0;
        for (int n = 0; n < nucCoord.size(); n++) {
            for (int i = 0; i < atomInt1.exponent().size(); i++) {
                for (int j = 0; j < atomInt2.exponent().size(); j++) {
                    resultado += normOrbital(atomInt1.exponent().get(i), atomInt1.cartAngular())
                            * normOrbital(atomInt2.exponent().get(j), atomInt2.cartAngular())
                            * atomInt1.coefficient().get(i) * atomInt2.coefficient().get(j)
                            * potInt(atomInt1.exponent().get(i), atomInt2.exponent().get(j),
                            atomInt1.coord(), atomInt2.coord(), nucCoord.get(n),
                            atomInt1.cartAngular(), atomInt2.cartAngular());
                }
            }
            resultado *= -numAtomic.get(n);
        }

        return resultado;
    }

    public static double ovInt(double orbC1, double orbC2,
                               Atom.Coordinates coord1, Atom.Coordinates coord2,
                               Atom.AngCoordinates angCoord1, Atom.AngCoordinates angCoord2) {

        double producto = 1.0;

        // Product over x, y, z components
        producto *= sRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x());
        producto *= sRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y());
        producto *= sRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        // Compute distance squared between coord1 and coord2
        double dist2 = Math.pow(coord1.x() - coord2.x(), 2) +
                Math.pow(coord1.y() - coord2.y(), 2) +
                Math.pow(coord1.z() - coord2.z(), 2);

        double eAB = Math.exp(-(orbC1 * orbC2 / (orbC1 + orbC2)) * dist2);

        return eAB * Math.pow(Math.PI / (orbC1 + orbC2), 1.5) * producto;
    }

    private static double kinInt(double orbC1, double orbC2,
                               Atom.Coordinates coord1, Atom.Coordinates coord2,
                               Atom.AngCoordinates angCoord1, Atom.AngCoordinates angCoord2) {


        double producto1 = kiRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x())
                * sRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y())
                * sRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        double producto2 = sRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x())
                * kiRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y())
                * sRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        double producto3 = sRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x())
                * sRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y())
                * kiRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        double suma = producto1 + producto2 + producto3;

        // Distance squared between centers
        double dist2 = Math.pow(coord1.x() - coord2.x(), 2) +
                Math.pow(coord1.y() - coord2.y(), 2) +
                Math.pow(coord1.z() - coord2.z(), 2);

        double energyAB = Math.exp(-(orbC1 * orbC2 / (orbC1 + orbC2)) * dist2);

        return energyAB * Math.pow(Math.PI / (orbC1 + orbC2), 1.5) * suma;

    }

    private static double potInt(double orbC1, double orbC2,
                                 Atom.Coordinates coord1, Atom.Coordinates coord2, Atom.Coordinates coordNuc,
                                 Atom.AngCoordinates angCoord1, Atom.AngCoordinates angCoord2) {

        // Distance squared between centers
        double dist2 = Math.pow(coord1.x() - coord2.x(), 2) +
                Math.pow(coord1.y() - coord2.y(), 2) +
                Math.pow(coord1.z() - coord2.z(), 2);

        double enerAB = Math.exp(-(orbC1 * orbC2 / (orbC1 + orbC2)) * dist2);

        return enerAB * 2.0 * Math.PI / (orbC1 + orbC2);

    }

    private static double sRecurr(double alpha1, double beta1,
                                 double coord1, double coord2,
                                 int angCoord1coord, int angCoord2coord) {

        // Base cases
        if (angCoord1coord == 0 && angCoord2coord == 0) {
            return 1.0;
        } else if (angCoord1coord == 1 && angCoord2coord == 0) {
            return -(coord1 - ((alpha1 * coord1 + beta1 * coord2) / (alpha1 + beta1)));
        } else if (angCoord1coord > 1 && angCoord2coord == 0) {
            double term1 = -(coord1 - ((alpha1 * coord1 + beta1 * coord2) / (alpha1 + beta1)))
                    * sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord - 1, angCoord2coord);
            double term2 = ((angCoord1coord - 1) / (2.0 * (alpha1 + beta1)))
                    * sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord - 2, angCoord2coord);
            return term1 + term2;
        } else {
            return sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord + 1, angCoord2coord - 1)
                    + (coord1 - coord2)
                    * sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord, angCoord2coord - 1);
        }
    }

    private static double kiRecurr(double alpha, double beta,
                                        double coord1, double coord2,
                                        int angCoord1coord, int angCoord2coord) {

        double resultado;

        // Case 1: both angular momenta are zero
        if (angCoord1coord == 0 && angCoord2coord == 0) {
            resultado = 2.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, 1, 1);

            // Case 2: angCoord1coord > 0, angCoord2coord== 0
        } else if (angCoord1coord > 0 && angCoord2coord == 0) {
            resultado = -angCoord1coord * beta * sRecurr(alpha, beta, coord1, coord2, angCoord1coord - 1, 1)
                    + 2.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, angCoord1coord + 1, 1);

            // Case 3: angCoord1coord == 0, angCoord2coord > 0
        } else if (angCoord1coord == 0 && angCoord2coord > 0) {
            resultado = -angCoord2coord * alpha * sRecurr(alpha, beta, coord1, coord2, 1, angCoord2coord - 1)
                    + 2.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, 1, angCoord2coord + 1);

            // Case 4: both angCoord1coord and CA2 > 0
        } else {
            resultado = (
                    (angCoord1coord * angCoord2coord * sRecurr(alpha, beta, coord1, coord2, angCoord1coord - 1, angCoord2coord - 1))
                            - (2.0 * angCoord1coord * beta * sRecurr(alpha, beta, coord1, coord2, angCoord1coord - 1, angCoord2coord + 1))
                            - (2.0 * angCoord2coord * alpha * sRecurr(alpha, beta, coord1, coord2, angCoord1coord + 1, angCoord2coord - 1))
                            + (4.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, angCoord1coord + 1, angCoord2coord + 1))
            ) / 2.0;
        }

        return resultado;
    }

    private static double potRecurr(double alpha, double beta,
                                    double coord1, double coord2,
                                    int angCoord1coord, int angCoord2coord,
                                    double nucCoord, double xX) {

        // Case 1: both angular momenta are zero
        if (angCoord1coord == 0 && angCoord2coord == 0) {
            return 1.0;
        } else if (angCoord1coord == 1 && angCoord2coord == 0) {
            return -(coord1 - ((alpha * coord1 + beta * coord2) / (alpha + beta))) +
                    Math.pow((xX+1)/2,2)*((alpha * coord1 + beta * coord2) / (alpha + beta)- nucCoord);
        } else if (angCoord1coord > 0 && angCoord2coord == 0) {
            return -(coord1 - ((alpha * coord1 + beta * coord2) / (alpha + beta))) +
                    Math.pow((xX+1)/2,2)*((alpha * coord1 + beta * coord2) / (alpha + beta)- nucCoord)
                    * potRecurr(alpha, beta, coord1, coord2, angCoord1coord - 1, angCoord2coord, nucCoord, xX)
                    +(angCoord1coord - 1)/(2.0*(alpha + beta))*(1 - Math.pow((xX+1)/2,2))
                    * potRecurr(alpha, beta, coord1, coord2, angCoord1coord - 2, angCoord2coord, nucCoord, xX);
        } else {
            return potRecurr(alpha, beta, coord1, coord2, angCoord1coord + 1, angCoord2coord - 1, nucCoord, xX)
                    + (coord1 - coord2)
                    * potRecurr(alpha, beta, coord1, coord2, angCoord1coord, angCoord2coord - 1, nucCoord, xX);
        }

    }

    public static double normOrbital(double orbC, Atom.AngCoordinates vecA) {

        double pi = Math.PI;
        double b = Math.pow((2.0 * orbC / pi), 0.75);
        b *= Math.pow((4 * orbC), (vecA.x() + vecA.y() + vecA.z()) / 2.0);

        b /= Math.sqrt(
                doubleFactorial(2 * vecA.x() - 1)
                        * doubleFactorial(2 * vecA.y() - 1)
                        * doubleFactorial(2 * vecA.z() - 1)
        );

        return b;
    }

    public static double doubleFactorial(int n) {
        if (n <= 0) return 1.0;
        double result = 1.0;
        for (int i = n; i > 0; i -= 2) {
            result *= i;
        }
        return result;
    }

    public static double abscissa(int n, int i) {
        double theta = i * Math.PI / (n + 1.0);

        double term1 = (n + 1.0 - 2.0 * i) / (n + 1.0);
        double term2 = (2.0 / Math.PI) * (1.0 + (2.0 / 3.0) * Math.pow(Math.sin(theta), 2))
                * Math.sin(theta) * Math.cos(theta);
        return term1 + term2;
    }

    public static double omega(int n, int i) {
        return (16.0 / (3.0*(n + 1)) * Math.pow(Math.sin((i * Math.PI) / (n + 1)), 4));
    }

//    public static double IntChebyshev(double eps, int M) {
//        final double PI = Math.PI;
//
//        int n = 3;
//        int j = 0;
//
//        double C0 = Math.cos(PI / 6.0);
//        double S0 = Math.sin(PI / 6.0);
//
//        double C1 = S0;
//        double S1 = C0;
//
//        double q = (func1(abscissa(2, 1)) + func1(-abscissa(2, 1))) * omega(2, 1);
//        double p = func1(0.0);
//        double CHP = p + q;
//
//        double error = 10.0;
//
//        while (error > eps && ((2 * n * (1 - j) + j * 4 * n / 3 - 1) <= M)) {
//            j = 1 - j;
//
//            C1 = j * C1 + (1 - j) * C0;
//            S1 = j * S1 + (1 - j) * S0;
//            C0 = j * C0 + (1 - j) * Math.sqrt((1.0 + C0) * 0.5);
//            S0 = j * S0 + (1 - j) * (S0 / (C0 + C0));
//
//            double C = C0;
//            double S = S0;
//
//            // Inner loop
//            for (int i = 1; i <= n - 1; i += 2) {
//                double x = 1 + 2.0 / (3.0 * PI) * S * C * (3 + 2 * S * S) - ((double) i / (double) n);
//
//                if (3 * Math.rint((i + j + j) / 3.0) > (i + j)) {
//                    CHP += (func1(-x) + func1(x)) * Math.pow(S, 4);
//                }
//
//                double temp = S;
//                S = S * C1 + C * S1;
//                C = C * C1 - temp * S1;
//            }
//
//            n = (1 + j) * n;
//            p = p + (1 - j) * (CHP - q);
//
//            // Error estimate
//            error = 16 * Math.abs((1 - j) * (q - 3 * p / 2.0) + j * (CHP - 2 * q)) / (3.0 * n);
//            q = (1 - j) * q + j * CHP;
//        }
//
//        CHP = 16 * q / (3.0 * n);
//        return CHP;
//    }
//
//    private static double func1(double x, double alpha, double beta,
//                                double[] coord1, double[] coord2, double[] coordNuc) {
//        double[] cons1 = new double[3];
//        for (int i = 0; i < 3; i++) {
//            cons1[i] = ((alpha * coord1[i] + beta * coord2[i]) / (alpha + beta)) - coordNuc[i];
//        }
//        return 0.5*Math.exp(-((alpha+beta)*(x+1)/2*dotProduct(cons1,cons1))/4.0);
//    }
//
//    /**
//     * Computes the dot product of two vectors of the same length.
//     *
//     * @param a first vector
//     * @param b second vector
//     * @return dot product a · b
//     * @throws IllegalArgumentException if vectors have different lengths
//     */
//    public static double dotProduct(double[] a, double[] b) {
//        if (a.length != b.length) {
//            throw new IllegalArgumentException("Vectors must have the same length");
//        }
//
//        double sum = 0.0;
//        for (int i = 0; i < a.length; i++) {
//            sum += a[i] * b[i];
//        }
//        return sum;
//    }

}
